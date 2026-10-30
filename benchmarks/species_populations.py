# SPDX-License-Identifier: LGPL-3.0-or-later
"""Fresh-process baseline/candidate comparisons of complete species writing."""

import argparse
import gc
import hashlib
import importlib.metadata
import json
import os
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

CASES = {
    "small": (32, 4, 12),
    "sparse": (200000, 64, 128),
    "sparse-large": (200000, 128, 256),
    "dense": (32768, 64, 128),
    "mixed": (32768, 128, 256),
    "many-species": (8, 4096, 8192),
}


def _sample(case):
    """Measure one complete write in a fresh process, outside fixture and hash work."""
    import numpy as np

    from reacnetgenerator._matrix import _GenerateMatrix
    from reacnetgenerator._step3state import _MoleculeNameTable
    from reacnetgenerator.utils import listtobytes

    frames, species, molecules = CASES[case]
    random = np.random.default_rng(20261030)
    names = _MoleculeNameTable.from_names(f"S{i % species}" for i in range(molecules))
    with tempfile.TemporaryDirectory(prefix="rng-species-benchmark-") as root:
        root = Path(root)
        source = root / "records.bin"
        with source.open("wb") as handle:
            for molecule in range(molecules):
                if case == "dense" or (case == "mixed" and molecule < 32):
                    observed = np.arange(frames)
                else:
                    visits = (
                        128 if case == "sparse-large" else (8 if case == "small" else 2)
                    )
                    observed = random.integers(0, frames, visits)
                for value in ([molecule], [], [], observed):
                    handle.write(listtobytes(value))
        del observed, value
        gc.collect()
        matrix = object.__new__(_GenerateMatrix)
        matrix.mname, matrix.timestep = names, np.arange(frames) * 10 + 100
        matrix.moleculetemp2filename = str(source)
        matrix.speciesfilename = str(root / "output.species")
        initial = resource.getrusage(resource.RUSAGE_SELF)
        started = time.perf_counter()
        diagnostics = matrix._printspecies()
        wall = time.perf_counter() - started
        final = resource.getrusage(resource.RUSAGE_SELF)
        scale = 1 if sys.platform == "darwin" else 1024
        digest = hashlib.sha256()
        with Path(matrix.speciesfilename).open("rb") as output:
            while chunk := output.read(1024 * 1024):
                digest.update(chunk)
        return {
            "wall_seconds": wall,
            "cpu_seconds": final.ru_utime
            + final.ru_stime
            - initial.ru_utime
            - initial.ru_stime,
            "peak_mib": final.ru_maxrss * scale / 2**20,
            "sha256": digest.hexdigest(),
            "output_bytes": Path(matrix.speciesfilename).stat().st_size,
            "diagnostics": diagnostics,
            "environment": {
                "python": sys.version.split()[0],
                "platform": platform.platform(),
                "numpy": importlib.metadata.version("numpy"),
                "source": sys.modules[_GenerateMatrix.__module__].__file__,
            },
        }


def main():
    """Compare alternating fresh interpreters and require byte-identical output."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-python", type=Path)
    parser.add_argument("--candidate-python", type=Path, default=Path(sys.executable))
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=list(CASES))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--worker", choices=CASES, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        args.output.write_text(json.dumps(_sample(args.worker), indent=2))
        return
    if args.baseline_python is None or args.repeats < 1:
        parser.error("A baseline interpreter and positive repeat count are required")
    environment = dict(
        os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1"
    )
    environment.pop("PYTHONPATH", None)
    rows = []
    with tempfile.TemporaryDirectory(prefix="rng-species-samples-") as root:
        for case in args.cases:
            samples = {"baseline": [], "candidate": []}
            expected = None
            for repeat in range(args.repeats):
                order = (
                    ("baseline", "candidate")
                    if repeat % 2 == 0
                    else ("candidate", "baseline")
                )
                for revision in order:
                    output = Path(root) / "sample.json"
                    command = [
                        str(getattr(args, f"{revision}_python").absolute()),
                        str(Path(__file__).resolve()),
                        "--worker",
                        case,
                        "--output",
                        str(output),
                    ]
                    subprocess.run(command, cwd=root, env=environment, check=True)
                    sample = json.loads(output.read_text())
                    if expected is None:
                        expected = sample["sha256"]
                    if sample["sha256"] != expected:
                        raise RuntimeError(f"Output mismatch: {case}, {revision}")
                    samples[revision].append(sample)
            envs = [values[0]["environment"] for values in samples.values()]
            if envs[0]["source"] == envs[1]["source"] or any(
                envs[0][key] != envs[1][key] for key in ("python", "platform", "numpy")
            ):
                raise RuntimeError(
                    "Comparison requires separate sources and matched environments"
                )
            medians = {
                revision: {
                    key: statistics.median(sample[key] for sample in values)
                    for key in ("wall_seconds", "cpu_seconds", "peak_mib")
                }
                for revision, values in samples.items()
            }
            rows.append({"case": case, "medians": medians, "samples": samples})
            args.output.write_text(json.dumps(rows, indent=2))
            print(case, medians, flush=True)


if __name__ == "__main__":
    main()
