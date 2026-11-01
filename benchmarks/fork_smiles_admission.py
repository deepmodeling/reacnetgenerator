# SPDX-License-Identifier: LGPL-3.0-or-later
"""Replay fork admission controls with explicit input warming in each process.

Real states use isomer_candidates.py --prepare. Synthetic controls reuse
smiles_batches.py, including oversized returned frames. Audits serialize tasks
again and must be run separately from timing samples. Use identical arguments
and this script with both installed source versions.
"""

import argparse
import hashlib
import importlib.metadata
import json
import multiprocessing
import platform
import tempfile
from pathlib import Path
from types import SimpleNamespace

import isomer_candidates
import smiles_batches

from reacnetgenerator import _smilesworker


def warm(path):
    """Pre-read compressed records outside the measured production entry point."""
    with open(path, "rb") as stream:
        for _ in iter(lambda: stream.read(1024 * 1024), b""):
            pass


def main():
    """Measure one fresh process and retain exact ordered output fingerprints."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        choices=["real", "repeated", "low-repeat", "mixed", "large-output"],
        default="real",
    )
    parser.add_argument("--state", type=Path)
    parser.add_argument(
        "--records",
        type=int,
        default=1024,
        help="Synthetic record count; real states replay all records",
    )
    parser.add_argument("--nproc", type=int, default=32)
    parser.add_argument("--miso", type=int, choices=[0, 1, 2], default=1)
    parser.add_argument(
        "--method",
        choices=multiprocessing.get_all_start_methods(),
        default=multiprocessing.get_start_method(),
    )
    parser.add_argument("--audit", action="store_true")
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    if args.nproc < 1 or args.records < 1:
        parser.error("--nproc and --records must be positive")
    if args.case == "real" and args.state is None:
        parser.error("Real replays require --state")
    multiprocessing.set_start_method(args.method, force=True)
    with tempfile.TemporaryDirectory() as directory:
        folder = Path(directory)
        if args.case == "real":
            audit = smiles_batches._audit() if args.audit else None
            result = isomer_candidates.replay(
                SimpleNamespace(
                    case="real",
                    state=args.state,
                    nproc=args.nproc,
                    miso=args.miso,
                    audit=False,
                ),
                folder,
            )
            result["audit"] = audit
        else:
            original_fixture = smiles_batches._fixture

            def prepared(options, root):
                """Warm the same prepared stream for either installed version."""
                rng, source = original_fixture(options, root)
                warm(rng.moleculetemp2filename)
                return rng, source

            smiles_batches._fixture = prepared
            try:
                result = smiles_batches.sample(
                    SimpleNamespace(
                        case=args.case,
                        records=args.records,
                        state=None,
                        nproc=args.nproc,
                        miso=args.miso,
                        hmm=False,
                        stage="results" if args.case == "large-output" else "names",
                        audit=args.audit,
                    ),
                    folder,
                )
            finally:
                smiles_batches._fixture = original_fixture
            result["seconds"] = result["wall"]
        result.update(
            method=args.method,
            batching=_smilesworker._should_batch_smiles(result["records"], args.nproc),
            worker_source_sha256=hashlib.sha256(
                Path(_smilesworker.__file__).read_bytes()
            ).hexdigest(),
            python=platform.python_version(),
            dependencies={
                name: importlib.metadata.version(name)
                for name in ["numpy", "rdkit", "networkx", "lz4"]
            },
        )
        args.json.write_text(json.dumps(result, indent=2))
        print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
