# SPDX-License-Identifier: LGPL-3.0-or-later
"""Command-line validation and semantic comparison for timeline artifacts."""

import argparse
import json
import os
import sys
import tempfile
from dataclasses import asdict
from pathlib import Path

from .timedoutput import (
    TimedOutputValidationError,
    compare_semantic_manifests,
    semantic_manifest,
    validate_timed_output,
)


def _parser():
    parser = argparse.ArgumentParser(
        description="Validate a ReacNetGenerator timeline and compare its meaning."
    )
    parser.add_argument("timeline", help="timeline.h5 artifact to validate")
    parser.add_argument(
        "--write-manifest",
        metavar="PATH",
        help="atomically write the validated semantic manifest as JSON",
    )
    parser.add_argument(
        "--compare-manifest",
        metavar="PATH",
        help="return exit code 1 if the artifact differs from this manifest",
    )
    parser.add_argument(
        "--include-provenance",
        action="store_true",
        help="include paths, source metadata, build version, and creation time",
    )
    parser.add_argument(
        "--block-rows",
        type=int,
        default=8192,
        help="maximum numeric rows read per validation/hash block (default: 8192)",
    )
    return parser


def _write_json_atomic(path, value):
    destination = Path(path).expanduser()
    handle, temporary = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        suffix=".incomplete",
        dir=destination.parent,
    )
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as file:
            json.dump(value, file, ensure_ascii=False, sort_keys=True, indent=2)
            file.write("\n")
            file.flush()
            os.fsync(file.fileno())
        os.replace(temporary, destination)
    except BaseException:
        try:
            os.close(handle)
        except OSError:
            pass
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def main(argv=None):
    """Run validation; return 0 for success and 1 for invalid or different data."""
    args = _parser().parse_args(argv)
    try:
        timeline = Path(args.timeline).expanduser().resolve()
        manifest_path = (
            Path(args.write_manifest).expanduser().resolve()
            if args.write_manifest
            else None
        )
        comparison_path = (
            Path(args.compare_manifest).expanduser().resolve()
            if args.compare_manifest
            else None
        )
        if manifest_path == timeline:
            raise ValueError("--write-manifest must not overwrite the input timeline")
        if manifest_path is not None and manifest_path == comparison_path:
            raise ValueError(
                "--write-manifest must not overwrite the comparison manifest"
            )
        manifest = None
        if args.write_manifest or args.compare_manifest:
            manifest = semantic_manifest(
                timeline,
                include_provenance=args.include_provenance,
                block_rows=args.block_rows,
            )
            summary = manifest["summary"]
        else:
            summary = asdict(
                validate_timed_output(timeline, block_rows=args.block_rows)
            )

        if args.compare_manifest:
            assert comparison_path is not None
            with comparison_path.open(encoding="utf-8") as file:
                reference = json.load(file)
            differences = compare_semantic_manifests(reference, manifest)
            if differences:
                for path in differences:
                    print(f"manifest differs at {path}", file=sys.stderr)
                return 1
        if args.write_manifest:
            _write_json_atomic(manifest_path, manifest)
        status = "equivalent" if args.compare_manifest else "valid"
        print(json.dumps({"status": status, "summary": summary}, sort_keys=True))
        return 0
    except (OSError, TimedOutputValidationError, ValueError) as exc:
        print(f"timed-output validation failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
