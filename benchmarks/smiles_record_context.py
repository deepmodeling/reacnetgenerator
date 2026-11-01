# SPDX-License-Identifier: LGPL-3.0-or-later
"""Replay per-record SMILES with a real RNG owner and frame-sized metadata.

Run the same script with base/head PYTHONPATH, alternating fresh processes.
Synthetic records isolate context size from molecule graph size. Optional real
states use isomer_candidates.py --prepare; --records limits an ordered prefix.
Audits measure actual serialized calls and run separately from timing samples.
"""

import argparse
import hashlib
import importlib.metadata
import json
import multiprocessing
import pickle
import platform
import tempfile
import time
from contextlib import closing
from pathlib import Path

import numpy as np
from smiles_batches import _usage

from reacnetgenerator import ReacNetGenerator, _path, utils
from reacnetgenerator import _smilesworker as workers


def main():
    """Measure conversion or full naming, including worker startup and cleanup."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frames", type=int, default=100000)
    parser.add_argument("--records", type=int, default=400)
    parser.add_argument("--nproc", type=int, default=4)
    parser.add_argument("--miso", type=int, choices=[0, 1, 2], default=1)
    parser.add_argument("--state", type=Path)
    parser.add_argument("--large-output", action="store_true")
    parser.add_argument("--audit", action="store_true")
    parser.add_argument(
        "--method", choices=multiprocessing.get_all_start_methods(), default="spawn"
    )
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    if min(args.frames, args.records, args.nproc) < 1:
        parser.error("Frame, record and worker counts must be positive")
    multiprocessing.set_start_method(args.method, force=True)
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        source = root / "input.bond"
        source.touch()
        rng = ReacNetGenerator(
            inputfilename=source,
            inputfiletype="bond",
            atomname=["C"],
            nproc=args.nproc,
            miso=args.miso,
            runHMM=False,
            output_dir=root / "output",
            printmoleculetime=args.large_output,
        )
        rng.N = 256
        rng.atomtype = np.zeros(rng.N, dtype=int)
        rng.step = args.frames
        rng.timestep = dict(enumerate(range(args.frames)))
        rng.source_frame_counts = [args.frames]
        records = root / "records.bin"
        if args.state:
            metadata = json.loads((args.state / "metadata.json").read_text())
            rng.atomname = np.asarray(metadata["atomname"])
            names = np.load(args.state / "atomnames.npy", allow_pickle=False)
            rng.N = len(names)
            lookup = {name: index for index, name in enumerate(rng.atomname)}
            rng.atomtype = np.asarray([lookup[name] for name in names])
            count = min(args.records, metadata["hmmit"])
            from reacnetgenerator._moleculepresence import _iter_molecule_records

            with (args.state / "molecules.bin").open("rb") as stream:
                with records.open("wb") as output:
                    for index, record in enumerate(_iter_molecule_records(stream)):
                        if index == count:
                            break
                        output.write(b"".join(record))
        else:
            count = args.records
            with records.open("wb") as output:
                for index in range(count):
                    atoms = np.asarray([index % 128 * 2, index % 128 * 2 + 1])
                    frames = (
                        np.arange(args.frames)
                        if args.large_output and index in (0, count // 2)
                        else np.asarray([index % args.frames])
                    )
                    for value in (
                        atoms,
                        atoms.reshape(1, 2),
                        np.asarray([1 + index % 3]),
                        frames,
                    ):
                        output.write(utils.listtobytes(value))
        rng.hmmit = count
        rng.moleculetemp2filename = str(records)
        collector = _path._CollectSMILESPaths(rng)
        collector.atomnames = collector.atomname[collector.atomtype]
        with records.open("rb") as stream:
            for _ in iter(lambda: stream.read(1024 * 1024), b""):
                pass
        audit = {
            "tasks": 0,
            "input_bytes": 0,
            "max_input_bytes": 0,
            "converter_initializer_bytes": 0,
        }
        if args.audit:
            original = utils._forking_pickle_dumps
            original_run = utils.run_mp

            def measured_run(nproc, **kwargs):
                if kwargs.get("initializer") is not None:
                    audit["converter_initializer_bytes"] = len(
                        original((kwargs["initializer"], kwargs["initargs"]))
                    )
                return original_run(nproc, **kwargs)

            _path.run_mp = workers.run_mp = measured_run

            def measured(value):
                payload = original(value)
                audit["tasks"] += 1
                audit["input_bytes"] += len(payload)
                audit["max_input_bytes"] = max(audit["max_input_bytes"], len(payload))
                return payload

            utils._forking_pickle_dumps = measured
        digest = hashlib.sha256()
        before = _usage()
        started = time.perf_counter()
        if args.large_output:
            with records.open("rb") as stream:
                with closing(
                    collector._getSMILESresults(
                        stream, collector._calmoleculeSMILESname
                    )
                ) as results:
                    for result in results:
                        digest.update(
                            pickle.dumps(result, protocol=pickle.HIGHEST_PROTOCOL)
                        )
        else:
            collector._printmoleculename()
        elapsed = time.perf_counter() - started
        after = _usage()
        for key in ("parent_cpu", "child_cpu"):
            after[key] -= before[key]
        outputs = {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (root / "output").glob("*")
            if path.is_file()
        }
        result = {
            "seconds": elapsed,
            "frames": args.frames,
            "records": count,
            "nproc": args.nproc,
            "miso": args.miso,
            "method": args.method,
            "batched": workers._should_batch_smiles(count, args.nproc),
            "real_state": args.state is not None,
            "large_output": args.large_output,
            "audit": audit if args.audit else None,
            "outputs": outputs,
            "results": digest.hexdigest() if args.large_output else None,
            **after,
            "python": platform.python_version(),
            "dependencies": {
                name: importlib.metadata.version(name)
                for name in ("numpy", "rdkit", "lz4", "networkx")
            },
            "source_sha256": {
                name: hashlib.sha256(
                    (Path(_path.__file__).parent / name).read_bytes()
                ).hexdigest()
                for name in ("_path.py", "_smilesworker.py", "utils.py")
            },
        }
        args.json.write_text(json.dumps(result, indent=2))
        print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
