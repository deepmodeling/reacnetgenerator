# SPDX-License-Identifier: LGPL-3.0-or-later
"""Measure installed SMILES scheduling with reproducible inputs and optional IPC audit.

Use --prepare --input trajectory --state directory once for real-record replays.
Then invoke the same script with separate base/head interpreters, alternating
order. Synthetic cases need no external input. Audit runs measure payloads and
must be kept separate from performance samples because they serialize twice.
"""

import argparse
import hashlib
import importlib.metadata
import json
import multiprocessing
import pickle
import platform
import resource
import shutil
import sys
import tempfile
import time
from contextlib import closing
from multiprocessing.reduction import ForkingPickler
from pathlib import Path


def _hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def _usage():
    parent = resource.getrusage(resource.RUSAGE_SELF)
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    scale = 2**20 if sys.platform == "darwin" else 1024
    return {
        "parent_cpu": parent.ru_utime + parent.ru_stime,
        "child_cpu": child.ru_utime + child.ru_stime,
        "parent_peak_mib": parent.ru_maxrss / scale,
        "child_peak_mib": child.ru_maxrss / scale,
    }


def _outputs(rng, folder):
    for key in (
        "moleculefilename",
        "moleculetimelinefilename",
        "atomroutefilename",
        "reactionfilename",
        "reactionabcdfilename",
        "speciesfilename",
        "reactioneventfilename",
    ):
        setattr(rng, key, str(folder / key))
    rng.timed_output = str(folder / "timeline.h5")


def prepare(args):
    """Save one portable, fixed detection stream for both versions to replay."""
    from reacnetgenerator import ReacNetGenerator
    from reacnetgenerator._detect import _Detect
    from reacnetgenerator._hmmfilter import _HMMFilter

    args.state.mkdir(parents=True, exist_ok=False)
    source = args.state / "input.lammpstrj"
    shutil.copyfile(args.input, source)
    rng = ReacNetGenerator(
        inputfilename=source,
        inputfiletype="dump",
        atomname=["C", "H", "O"],
        nproc=1,
        runHMM=False,
        timed_output=args.state / "prepare.h5",
    )
    _Detect.gettype(rng).detect()
    original = Path(rng.moleculetempfilename)
    target = args.state / "records.bin"
    shutil.copyfile(original, target)
    original.unlink()
    rng.moleculetempfilename = str(target)
    _HMMFilter(rng).filter()
    original = Path(rng.originfilename)
    shutil.copyfile(original, args.state / "origin.bin")
    original.unlink()
    (args.state / "state.pkl").write_bytes(pickle.dumps(rng))
    print(
        json.dumps(
            {
                "frames": rng.step,
                "atoms": rng.N,
                "records": rng.hmmit,
                "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            }
        )
    )


def _fixture(args, root):
    """Create bounded synthetic records, or repeat a fixed real detection stream."""
    import numpy as np

    from reacnetgenerator import ReacNetGenerator
    from reacnetgenerator._hmmfilter import _HMMFilter
    from reacnetgenerator._moleculepresence import _iter_molecule_records
    from reacnetgenerator.utils import listtobytes

    if args.case == "real":
        rng = pickle.loads((args.state / "state.pkl").read_bytes())
        source = args.state / "input.lammpstrj"
        rng.inputfilename = [str(source)]
        rng.moleculetempfilename = str(args.state / "records.bin")
        rng.moleculetemp2filename = rng.moleculetempfilename
        rng.originfilename = str(args.state / "origin.bin")
        if args.records:
            # Repeat/truncate in original order; duplicates are stress inputs,
            # not additional physical observations of the original trajectory.
            target = root / "records.bin"
            remaining = args.records
            with target.open("wb") as output:
                while remaining:
                    before = remaining
                    with open(rng.moleculetempfilename, "rb") as file:
                        for record in _iter_molecule_records(file):
                            output.write(b"".join(record))
                            remaining -= 1
                            if not remaining:
                                break
                    if remaining == before:
                        raise ValueError("Cannot repeat an empty detection stream")
            rng.moleculetempfilename = str(target)
            rng.moleculetemp2filename = str(target)
            rng.hmmit = rng.temp1it = args.records
    else:
        count = args.records or 16384
        source = root / "input.bond"
        source.touch()
        rng = ReacNetGenerator(
            inputfilename=source,
            inputfiletype="bond",
            atomname=["C", "O"],
            nproc=args.nproc,
            runHMM=False,
        )
        rng.N = 256
        rng.atomtype = np.arange(rng.N) % 2
        rng.step = 262144 if args.case == "large-output" else 64
        rng.timestep = (
            {} if args.case == "large-output" else dict(enumerate(range(rng.step)))
        )
        rng.source_frame_counts = [rng.step]
        rng.hmmit = rng.temp1it = count
        target = root / "records.bin"
        with target.open("wb") as file:
            for index in range(count):
                width = 4 if args.case == "repeated" else 14
                if args.case == "mixed":
                    width = (2, 8, 32, 128)[index % 4]
                atoms = np.asarray(
                    [2 * j + ((index >> (j % 14)) & 1) for j in range(width)]
                )
                frames = (
                    np.arange(rng.step)
                    if args.case == "large-output" and index in (0, count // 2)
                    else np.asarray([index % 64])
                )
                for value in (
                    atoms,
                    np.column_stack((atoms[:-1], atoms[1:])),
                    np.ones(width - 1, dtype=int),
                    frames,
                ):
                    file.write(listtobytes(value))
        rng.moleculetempfilename = rng.moleculetemp2filename = str(target)
        rng.originfilename = str(root / "origin.bin")
        Path(rng.originfilename).touch()
    rng.miso = args.miso
    rng.runHMM = args.hmm
    if args.hmm:
        # PATH/name replays use a fixed ordered filter preparation, outside the
        # timer. This is not a full multi-worker HMM pipeline measurement.
        rng.nproc = 1
        _HMMFilter(rng).filter()
    rng.nproc = args.nproc
    rng.printmoleculetime = args.case == "large-output"
    rng.moleculeframes = rng.moleculetimesteps = None
    return rng, source


def _audit():
    """Count application payload bytes separately from timed performance runs."""
    from reacnetgenerator import utils

    original_submit = utils._submit_pool_item
    original_result = utils._record_serialized_pool_result
    original_decode = utils._decode_pool_result
    metrics = dict.fromkeys(
        (
            "tasks",
            "input_bytes",
            "max_input_bytes",
            "result_bytes",
            "max_result_bytes",
            "max_pending",
            "pending",
            "spilled_bytes",
            "max_spilled_bytes",
        ),
        0,
    )

    def submit(pool, completed, pending, func, index, item):
        size = len(ForkingPickler.dumps((func, item)))
        metrics["tasks"] += 1
        metrics["input_bytes"] += size
        metrics["max_input_bytes"] = max(metrics["max_input_bytes"], size)
        metrics["pending"] += 1
        metrics["max_pending"] = max(metrics["max_pending"], metrics["pending"])
        return original_submit(pool, completed, pending, func, index, item)

    def result(completed, index, value):
        size = len(value[1])
        metrics["result_bytes"] += size
        metrics["max_result_bytes"] = max(metrics["max_result_bytes"], size)
        return original_result(completed, index, value)

    def decode(*args):
        metrics["pending"] -= 1
        return original_decode(*args)

    utils._submit_pool_item = submit
    utils._record_serialized_pool_result = result
    utils._decode_pool_result = decode
    try:
        from reacnetgenerator import _smilesworker
    except ImportError:
        return metrics
    original_run = _smilesworker.run_mp

    def run(nproc, **kwargs):
        folder = Path(kwargs["initargs"][-2])
        with closing(original_run(nproc, **kwargs)) as batches:
            for item in batches:
                index, _, payload = item
                if payload is None:
                    size = (folder / f"{index}.pickle").stat().st_size
                    metrics["spilled_bytes"] += size
                    metrics["max_spilled_bytes"] = max(
                        metrics["max_spilled_bytes"], size
                    )
                yield item

    _smilesworker.run_mp = run
    return metrics


def sample(args, root):
    """Time real entry points, then validate their complete outputs."""
    from reacnetgenerator import _path
    from reacnetgenerator.reacnetgen import _json_compatible
    from reacnetgenerator.timedoutput import semantic_manifest, validate_timed_output

    rng, source = _fixture(args, root)
    output = root / "output"
    output.mkdir()
    _outputs(rng, output)
    # Replays change runtime options after restoring DETECT state. Keep the
    # writer's provenance snapshot consistent with the configuration measured.
    rng.parameters = {
        key: _json_compatible(getattr(rng, key)) for key in rng.parameters
    }
    rng.explicit_parameters = sorted(
        set(rng.explicit_parameters) | {"nproc", "miso", "runHMM", "timed_output"}
    )
    collector = _path._CollectPaths.getstype(rng)
    collector.atomnames = collector.atomname[collector.atomtype]
    metrics = _audit() if args.audit else None
    raw_digest = hashlib.sha256()

    def results():
        # Large-return stress measures conversion plus incremental consumption;
        # hashing is included identically, without retaining all decoded arrays.
        with open(collector.moleculetemp2filename, "rb") as file:
            with closing(
                collector._getSMILESresults(file, collector._calmoleculeSMILESname)
            ) as stream:
                for name, atoms, bonds, frames in stream:
                    raw_digest.update(
                        pickle.dumps(
                            (name, atoms, bonds, frames),
                            protocol=pickle.HIGHEST_PROTOCOL,
                        )
                    )

    action = {
        "names": collector._printmoleculename,
        "results": results,
        "path": collector.collect,
    }[args.stage]
    before = _usage()
    started = time.perf_counter()
    action()
    wall = time.perf_counter() - started
    after = _usage()
    for key in ("parent_cpu", "child_cpu"):
        after[key] -= before[key]
    text = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in output.iterdir()
        if path.is_file() and path.suffix != ".h5"
    }
    manifest = None
    if args.stage == "path":
        validate_timed_output(rng.timed_output)
        manifest = _hash(semantic_manifest(rng.timed_output))
    package = Path(_path.__file__).parent
    return {
        "case": args.case,
        "stage": args.stage,
        "records": rng.hmmit,
        "nproc": args.nproc,
        "miso": args.miso,
        "hmm": args.hmm,
        "method": multiprocessing.get_start_method(),
        "wall": wall,
        **after,
        "audit": metrics,
        "semantic": {
            "text": text,
            "manifest": manifest,
            "results": raw_digest.hexdigest() if args.stage == "results" else None,
        },
        "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "package": str(package),
            "dependencies": {
                name: importlib.metadata.version(name)
                for name in ("numpy", "rdkit", "hmmlearn", "lz4", "h5py")
            },
            "sources": {
                name: hashlib.sha256((package / name).read_bytes()).hexdigest()
                for name in ("_path.py", "_smilesworker.py", "utils.py")
                if (package / name).exists()
            },
        },
    }


def main():
    """Use one fresh process per sample; clean internal files after validation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--input", type=Path)
    parser.add_argument("--state", type=Path)
    parser.add_argument(
        "--case",
        choices=("real", "repeated", "low-repeat", "mixed", "large-output"),
        default="real",
    )
    parser.add_argument(
        "--stage", choices=("names", "results", "path"), default="names"
    )
    parser.add_argument("--records", type=int, default=0)
    parser.add_argument("--nproc", type=int, default=2)
    parser.add_argument("--miso", type=int, choices=(0, 1, 2), default=0)
    parser.add_argument("--hmm", action="store_true")
    parser.add_argument("--audit", action="store_true")
    parser.add_argument(
        "--method", choices=multiprocessing.get_all_start_methods(), default="spawn"
    )
    args = parser.parse_args()
    if args.records < 0:
        parser.error("--records must be non-negative")
    if args.prepare:
        if not args.input or not args.state:
            parser.error("--prepare requires --input and --state")
    elif args.case == "real" and not args.state:
        parser.error("real replays require --state")
    elif args.case == "large-output" and args.stage != "results":
        parser.error("large-output requires --stage results")
    elif args.stage == "path" and args.case != "real":
        parser.error("PATH replays require real detection and origin records")
    elif args.audit and args.stage not in ("names", "results"):
        parser.error("audit only names/results; keep it separate from timing samples")
    multiprocessing.set_start_method(args.method, force=True)
    if args.prepare:
        prepare(args)
        return
    with tempfile.TemporaryDirectory(prefix="rng-smiles-sample-") as folder:
        previous = tempfile.tempdir
        tempfile.tempdir = folder
        try:
            result = sample(args, Path(folder))
        finally:
            tempfile.tempdir = previous
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
