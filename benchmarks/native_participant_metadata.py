# SPDX-License-Identifier: LGPL-3.0-or-later
"""Measure the full native-output pipeline on a bounded LAMMPS dump window."""

import argparse
import hashlib
import importlib.metadata
import json
import multiprocessing
import os
import platform
import time
from pathlib import Path


def extract_window(source, target, frames):
    """Stream an initial frame window without reading the rest of a trajectory."""
    count = 0
    digest = hashlib.sha256()
    with source.open("rb") as reader, target.open("xb") as writer:
        for line in reader:
            if line.startswith(b"ITEM: TIMESTEP"):
                count += 1
                if count > frames:
                    break
            writer.write(line)
            digest.update(line)
    if count < frames:
        raise ValueError(f"Requested {frames} frames, found only {count}")
    return digest.hexdigest()


def main():
    """Write fresh outputs, timing and a validated semantic manifest."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--atom-names", nargs="+", required=True)
    parser.add_argument("--frames", type=int, default=201)
    parser.add_argument("--nproc", type=int, default=32)
    parser.add_argument("--miso", type=int, choices=(0, 1, 2), default=0)
    parser.add_argument("--hmm", action="store_true")
    parser.add_argument(
        "--start-method",
        choices=multiprocessing.get_all_start_methods(),
        default="fork",
    )
    args = parser.parse_args()
    if args.frames < 1 or args.nproc < 1:
        parser.error("--frames and --nproc must be positive")
    source = args.trajectory.resolve(strict=True)
    folder = args.output_dir.absolute()
    folder.mkdir(parents=True, exist_ok=False)
    window_sha256 = extract_window(source, folder / "trajectory.lammpstrj", args.frames)
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"
    multiprocessing.set_start_method(args.start_method, force=True)
    import reacnetgenerator
    from reacnetgenerator import ReacNetGenerator
    from reacnetgenerator.timedoutput import semantic_manifest

    os.chdir(folder)
    # Input preparation, initialization, imports and validation are untimed.
    parameters = {
        "inputfilename": "trajectory.lammpstrj",
        "inputfiletype": "dump",
        "atomname": args.atom_names,
        "nproc": args.nproc,
        "miso": args.miso,
        "runHMM": args.hmm,
        "stepinterval": 1,
        "pbc": True,
        "SMILES": True,
        "split": 1,
        "needprintspecies": True,
        "getoriginfile": False,
        "printfiltersignal": False,
        "p": [0.5, 0.5],
        "a": [[0.999, 0.001], [0.001, 0.999]],
        "b": [[0.6, 0.4], [0.4, 0.6]],
        "timed_output": "trajectory.lammpstrj.timeline.h5",
    }
    rng = ReacNetGenerator(**parameters)
    started = time.perf_counter()
    rng.run()
    elapsed = time.perf_counter() - started
    package = Path(reacnetgenerator.__file__).parent
    measurement = {
        "pipeline_seconds": elapsed,
        "frames": int(rng.step),
        "atoms": int(rng.N),
        "parameters": parameters,
        "input_window_sha256": window_sha256,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "affinity": sorted(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity")
        else [],
        "start_method": args.start_method,
        "native_math_threads": 1,
        "dependencies": {
            name: importlib.metadata.version(name)
            for name in (
                "numpy",
                "scipy",
                "openbabel",
                "h5py",
                "rdkit",
                "ase",
                "hmmlearn",
                "lz4",
                "networkx",
            )
        },
        "source_sha256": {
            name: hashlib.sha256((package / name).read_bytes()).hexdigest()
            for name in (
                "_timedoutput.py",
                "_timedoutputrange.py",
                "_timedoutputvalidate.py",
            )
        },
        "native_extension_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in package.glob("dps.*")
            if path.suffix in {".so", ".pyd"}
        },
    }
    Path("measurement.json").write_text(json.dumps(measurement, indent=2) + "\n")
    print(json.dumps(measurement), flush=True)
    manifest = semantic_manifest("trajectory.lammpstrj.timeline.h5")
    Path("semantic-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
