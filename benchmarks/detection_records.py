# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compare complete DETECT runs and isolated accumulation/save replays."""

import argparse
import hashlib
import importlib.metadata
import itertools
import json
import multiprocessing
import os
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

# Stage, pattern, frames, atoms/definitions, workers. Replays retain real save
# workers in the baseline while replacing only the upstream detection producer.
CASES = {
    "replay-small": ("replay", "dense", 8, 4, 1),
    "replay-dense-1": ("replay", "dense", 65536, 128, 1),
    "replay-dense-2": ("replay", "dense", 65536, 128, 2),
    "replay-dense-4": ("replay", "dense", 65536, 128, 4),
    "replay-wide": ("replay", "dense", 131072, 64, 2),
    "replay-sparse": ("replay", "sparse", 32768, 256, 2),
    "replay-unique": ("replay", "unique", 8192, 32768, 2),
    "replay-unique-threshold": ("replay", "unique", 6144, 24576, 2),
    "detect-small-1": ("detect", "dense", 8, 3, 1),
    "detect-small-2": ("detect", "dense", 8, 3, 2),
    "detect-dense-1": ("detect", "dense", 512, 128, 1),
    "detect-dense-2": ("detect", "dense", 512, 128, 2),
    "detect-dense-4": ("detect", "dense", 512, 128, 4),
    "detect-changing-2": ("detect", "changing", 512, 256, 2),
}


def _usage():
    """Include parent and terminated-worker CPU; keep RSS high-water marks separate."""
    parent = resource.getrusage(resource.RUSAGE_SELF)
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    scale = 1 if sys.platform == "darwin" else 1024
    return {
        "cpu_seconds": parent.ru_utime
        + parent.ru_stime
        + child.ru_utime
        + child.ru_stime,
        "parent_peak_mib": parent.ru_maxrss * scale / 2**20,
        "child_peak_mib": child.ru_maxrss * scale / 2**20,
    }


def _sample(case):
    """Time detection or replay, then validate complete records and provenance."""
    import numpy as np

    from reacnetgenerator import ReacNetGenerator
    from reacnetgenerator import _detect as detection
    from reacnetgenerator.utils import bytestolist, listtobytes, read_compressed_block

    stage, pattern, frames, atoms, nproc = CASES[case]
    with tempfile.TemporaryDirectory(prefix="rng-detection-sample-") as folder:
        root = Path(folder)
        source = root / "input.bond"
        # Isolated H atoms exercise real parsing, connected components and
        # structure compression without chemistry-dependent graph ambiguity.
        body = "".join(f"{atom + 1} 1 0 0 0 0 0\n" for atom in range(atoms))
        random = np.random.default_rng(20261030)
        with source.open("w") as output:
            for frame in range(frames if stage == "detect" else 1):
                if pattern == "changing":
                    pairs = random.permutation(atoms).reshape((-1, 2))
                    partners = np.empty(atoms, dtype=int)
                    partners[pairs[:, 0]] = pairs[:, 1]
                    partners[pairs[:, 1]] = pairs[:, 0]
                    body = "".join(
                        f"{atom + 1} 1 1 {partner + 1} 0 1.0 1.0 0 0\n"
                        for atom, partner in enumerate(partners)
                    )
                output.write(
                    f"# Timestep {frame * 10}\n# Number of particles {atoms}\n{body}"
                )
        del body
        rng = ReacNetGenerator(
            inputfilename=source,
            inputfiletype="bond",
            atomname=["H"],
            nproc=nproc,
            output_dir=root,
            timed_output=root / "timeline.h5" if stage == "detect" else None,
        )
        original = detection.run_mp
        if stage == "replay":
            definitions = [
                b"".join(listtobytes(value) for value in ([atom], [], []))
                for atom in range(atoms)
            ]

            def producer():
                """Replay fixed occurrences without timing frame detection."""
                for frame in range(frames):
                    if pattern == "dense":
                        molecules = definitions
                    elif pattern == "unique":
                        molecules = definitions[4 * frame : 4 * frame + 4]
                    else:
                        molecules = [
                            definitions[(2 * frame + offset) % atoms]
                            for offset in (0, 1)
                        ]
                    yield molecules, (frame, frame * 10)

            def replay(nproc, **kwargs):
                """Replace detection while preserving the real save-stage scheduler."""
                if kwargs["desc"] == "Save molecules":
                    return original(nproc, **kwargs)
                return producer()

            detection.run_mp = replay
        try:
            before = _usage()
            started = time.perf_counter()
            detection._Detect.gettype(rng).detect()
            wall = time.perf_counter() - started
            after = _usage()
        finally:
            detection.run_mp = original
        path = Path(rng.moleculetempfilename)
        try:
            # Internal unsigned widths and worker arrival order may differ.
            # Canonicalize complete definitions and sorted frame multisets;
            # retain a separate order-sensitive digest for deterministic replay.
            semantic, ordered = [], hashlib.sha256()
            payload, widths = 0, Counter()
            with path.open("rb") as stream:
                blocks = iter(read_compressed_block(stream))
                for atom_blob, bond_blob, level_blob, frame_blob in zip(
                    *[blocks] * 4, strict=True
                ):
                    definition = [
                        bytestolist(blob) for blob in (atom_blob, bond_blob, level_blob)
                    ]
                    visits = bytestolist(frame_blob)
                    assert visits.dtype.kind in "iu"
                    assert np.all(visits >= 0) and np.all(visits < frames)
                    record = json.dumps(
                        [definition, visits.tolist()], separators=(",", ":")
                    ).encode()
                    ordered.update(record + b"\n")
                    record = json.dumps(
                        [definition, sorted(visits.tolist())], separators=(",", ":")
                    ).encode()
                    semantic.append(hashlib.sha256(record).hexdigest())
                    payload += visits.nbytes
                    widths[str(visits.dtype)] += 1
            metadata = [
                rng.N,
                rng.atomtype.tolist(),
                rng.step,
                sorted(rng.timestep.items()),
                rng.source_frame_counts,
            ]
            digest = hashlib.sha256(
                json.dumps([metadata, sorted(semantic)], separators=(",", ":")).encode()
            ).hexdigest()
            return {
                "wall_seconds": wall,
                "cpu_seconds": after["cpu_seconds"] - before["cpu_seconds"],
                "parent_peak_mib": after["parent_peak_mib"],
                "child_peak_mib": after["child_peak_mib"],
                "output_bytes": path.stat().st_size,
                "frame_payload_bytes": payload,
                "records": len(semantic),
                "frame_dtypes": dict(widths),
                "semantic_sha256": digest,
                "ordered_sha256": ordered.hexdigest(),
                "environment": {
                    "python": sys.version.split()[0],
                    "platform": platform.platform(),
                    "numpy": importlib.metadata.version("numpy"),
                    "lz4": importlib.metadata.version("lz4"),
                    "source_sha256": hashlib.sha256(
                        Path(detection.__file__).read_bytes()
                    ).hexdigest(),
                    "start_method": multiprocessing.get_start_method(),
                    "source": detection.__file__,
                },
            }
        finally:
            path.unlink()


def main():
    """Alternate independent installations and fail on semantic output drift."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-python", type=Path)
    parser.add_argument("--candidate-python", type=Path, default=Path(sys.executable))
    parser.add_argument(
        "--previous-python",
        type=Path,
        help="Optional fixed pre-fix PR installation for a three-way comparison",
    )
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument(
        "--start-method",
        choices=multiprocessing.get_all_start_methods(),
        default="spawn",
    )
    parser.add_argument("--cases", nargs="+", choices=CASES, default=list(CASES))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--worker", choices=CASES, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        multiprocessing.set_start_method(args.start_method, force=True)
        args.output.write_text(json.dumps(_sample(args.worker), indent=2))
        return
    if args.baseline_python is None or args.repeats < 1:
        parser.error("Separate baseline interpreter and positive repeat count required")
    if max(CASES[case][-1] for case in args.cases) > (os.cpu_count() or 1):
        parser.error("Selected worker count exceeds the available CPU count")
    environment = dict(
        os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1"
    )
    environment.pop("PYTHONPATH", None)
    revisions = (
        ("baseline", "previous", "candidate")
        if args.previous_python
        else ("baseline", "candidate")
    )
    orders = list(itertools.permutations(revisions))
    rows = []
    with tempfile.TemporaryDirectory(prefix="rng-detection-comparison-") as folder:
        for case in args.cases:
            samples = {revision: [] for revision in revisions}
            expected = None
            for repeat in range(args.repeats):
                # Six permutations balance each three-way revision in every
                # execution position; two-way runs retain alternating order.
                order = orders[repeat % len(orders)]
                for revision in order:
                    output = Path(folder) / "sample.json"
                    result = subprocess.run(
                        [
                            str(getattr(args, f"{revision}_python").absolute()),
                            str(Path(__file__).resolve()),
                            "--worker",
                            case,
                            "--start-method",
                            args.start_method,
                            "--output",
                            str(output),
                        ],
                        cwd=folder,
                        env=environment,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.PIPE,
                        text=True,
                    )
                    if result.returncode:
                        raise RuntimeError(
                            f"{case}, {revision}: {result.stderr[-4000:]}"
                        )
                    sample = json.loads(output.read_text())
                    digest = (
                        sample["semantic_sha256"],
                        sample["ordered_sha256"] if case.startswith("replay") else None,
                    )
                    if expected is None:
                        expected = digest
                    if digest != expected:
                        raise RuntimeError(f"Output mismatch: {case}, {revision}")
                    samples[revision].append(sample)
            envs = [values[0]["environment"] for values in samples.values()]
            if len({env["source"] for env in envs}) != len(envs) or any(
                envs[0][key] != env[key]
                for env in envs[1:]
                for key in ("python", "platform", "numpy", "lz4", "start_method")
            ):
                raise RuntimeError(
                    "Comparison requires separate sources and matched environments"
                )
            medians = {
                revision: {
                    key: statistics.median(sample[key] for sample in values)
                    for key in (
                        "wall_seconds",
                        "cpu_seconds",
                        "parent_peak_mib",
                        "child_peak_mib",
                    )
                }
                for revision, values in samples.items()
            }
            rows.append(
                {
                    "case": case,
                    "configuration": CASES[case],
                    "execution_orders": [
                        orders[i % len(orders)] for i in range(args.repeats)
                    ],
                    "medians": medians,
                    "samples": samples,
                }
            )
            args.output.write_text(json.dumps(rows, indent=2))
            print(case, medians, flush=True)


if __name__ == "__main__":
    main()
