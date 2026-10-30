# SPDX-License-Identifier: LGPL-3.0-or-later
"""Measure uncached single-process names/HMM and complete DETECT-to-MATRIX runs.

Invoke the same script with separately installed base and head interpreters.
Each invocation emits one JSON sample; imports, input preparation and validation
are outside the timer. Use --stage full to include detection and all core stages.
"""

import argparse
import hashlib
import importlib.metadata
import json
import multiprocessing
import platform
import resource
import sys
import tempfile
import time
from pathlib import Path


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def _usage():
    parent = resource.getrusage(resource.RUSAGE_SELF)
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    scale = 2**20 if sys.platform == "darwin" else 1024
    return {
        "parent_cpu_s": parent.ru_utime + parent.ru_stime,
        "child_cpu_s": child.ru_utime + child.ru_stime,
        "parent_peak_mib": parent.ru_maxrss / scale,
        "child_peak_mib": child.ru_maxrss / scale,
    }


def sample(args, root):
    """Exercise installed production entry points and hash their complete outputs."""
    from reacnetgenerator import ReacNetGenerator, _hmmfilter, _path
    from reacnetgenerator._detect import _Detect
    from reacnetgenerator.timedoutput import semantic_manifest, validate_timed_output
    from reacnetgenerator.utils import read_compressed_block

    source = args.input
    if source is None:
        source = root / "input.bond"
        # A small changing H/H2 trajectory keeps the benchmark self-contained.
        with source.open("w") as file:
            for frame in range(32):
                file.write(f"# Timestep {frame * 10}\n# Number of particles 8\n")
                for atom in range(8):
                    if 8 <= frame < 24:
                        partner = atom ^ 1
                        file.write(f"{atom + 1} 1 1 {partner + 1} 0 1.0 1.0 0 0\n")
                    else:
                        file.write(f"{atom + 1} 1 0 0 0 0 0\n")
    output = root / "output"
    output.mkdir()
    rng = ReacNetGenerator(
        inputfilename=str(source),
        inputfiletype=args.input_type if args.input else "bond",
        atomname=args.atomname if args.input else ["H"],
        nproc=args.nproc,
        miso=args.miso,
        runHMM=args.hmm,
        getoriginfile=args.origin,
        printfiltersignal=args.filtered,
        printmoleculetime=True,
        output_dir=output,
        timed_output=output / "timeline.h5",
    )
    if args.stage != "full":
        _Detect.gettype(rng).detect()
    if args.stage == "names":
        _hmmfilter._HMMFilter(rng).filter()
        collector = _path._CollectPaths.getstype(rng)
        collector.atomnames = collector.atomname[collector.atomtype]
        action = collector._printmoleculename
    elif args.stage == "hmm":
        action = _hmmfilter._HMMFilter(rng).filter
    else:
        action = rng.run
    before = _usage()
    started = time.perf_counter()
    action()
    wall = time.perf_counter() - started
    after = _usage()
    # Capture peaks before validation; these still include imports/preparation.
    metrics = {"wall_s": wall, **after}
    for key in ("parent_cpu_s", "child_cpu_s"):
        metrics[key] -= before[key]
    metrics["cpu_s"] = metrics["parent_cpu_s"] + metrics["child_cpu_s"]
    text = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in output.iterdir()
        if path.is_file() and path.suffix != ".h5"
    }
    manifest = None
    pairs = None
    if args.stage == "full":
        validate_timed_output(rng.timed_output)
        manifest = _digest(semantic_manifest(rng.timed_output))
    elif args.stage == "hmm" and args.hmm:
        # Compare entire records paired with BOTH optional signal streams;
        # parallel HMM completion order is allowed to differ.
        with open(rng.moleculetemp2filename, "rb") as file:
            blocks = read_compressed_block(file)
            records = list(zip(*[blocks] * 4, strict=True))
        with open(rng.hmmfilename, "rb") as file:
            signals = list(read_compressed_block(file))
        if rng.originfilename:
            with open(rng.originfilename, "rb") as file:
                origins = list(read_compressed_block(file))
        else:
            origins = [b""] * len(records)
        pairs = _digest(
            sorted(
                hashlib.sha256(b"".join(record) + origin + signal).hexdigest()
                for record, origin, signal in zip(
                    records, origins, signals, strict=True
                )
            )
        )
    elif args.stage == "hmm":
        # Dense/sparse no-HMM streams are ordered and reuse the same definitions.
        pairs = _digest(
            [
                hashlib.sha256(Path(name).read_bytes()).hexdigest()
                for name in (rng.moleculetemp2filename, rng.originfilename)
            ]
        )
    return {
        **metrics,
        "stage": args.stage,
        "hmm": args.hmm,
        "origin": args.origin,
        "filtered": args.filtered,
        "nproc": args.nproc,
        "miso": args.miso,
        "start_method": multiprocessing.get_start_method(),
        "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "text_sha256": text,
        "manifest_sha256": manifest,
        "record_signal_sha256": pairs,
        "records": rng.hmmit,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "path_source": _path.__file__,
            "hmm_source": _hmmfilter.__file__,
            "source_sha256": {
                Path(module.__file__).name: hashlib.sha256(
                    Path(module.__file__).read_bytes()
                ).hexdigest()
                for module in (_path, _hmmfilter)
            },
            "dependencies": {
                name: importlib.metadata.version(name)
                for name in ("numpy", "rdkit", "hmmlearn", "lz4", "h5py")
            },
        },
    }


def main():
    """Run one fresh sample, keeping internal temporary files inside its directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("names", "hmm", "full"), required=True)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--input-type", default="dump")
    parser.add_argument("--atomname", nargs="+", default=["C", "H", "O"])
    parser.add_argument("--hmm", action="store_true")
    parser.add_argument("--origin", action="store_true")
    parser.add_argument("--filtered", action="store_true")
    parser.add_argument("--nproc", type=int, default=1)
    parser.add_argument("--miso", type=int, choices=(0, 1, 2), default=0)
    parser.add_argument(
        "--start-method",
        choices=multiprocessing.get_all_start_methods(),
        default="spawn",
    )
    args = parser.parse_args()
    multiprocessing.set_start_method(args.start_method, force=True)
    with tempfile.TemporaryDirectory(prefix="rng-single-process-") as folder:
        previous = tempfile.tempdir
        tempfile.tempdir = folder
        try:
            result = sample(args, Path(folder))
        finally:
            tempfile.tempdir = previous
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
