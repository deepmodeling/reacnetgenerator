# SPDX-License-Identifier: LGPL-3.0-or-later
"""Measure installed numeric timeline writing and its integration costs.

Run the same script with separately installed base/head interpreters in AB/BA
order. Each invocation is a fresh process. Audit runs instrument buffer usage
and must be kept separate from timing samples.
"""

import argparse
import hashlib
import importlib.metadata
import json
import multiprocessing
import platform
import sys
import tempfile
import time
from pathlib import Path

from single_process import _digest, _usage


def _small_source(path):
    """Write a small changing H/H2 trajectory for complete parser-to-matrix runs."""
    with path.open("w") as file:
        for frame in range(32):
            file.write(f"# Timestep {frame * 10}\n# Number of particles 8\n")
            for atom in range(8):
                if 8 <= frame < 24:
                    file.write(f"{atom + 1} 1 1 {(atom ^ 1) + 1} 0 1.0 1.0 0 0\n")
                else:
                    file.write(f"{atom + 1} 1 0 0 0 0 0\n")


def _synthetic_records(rng, root, args):
    """Create table stress inputs, not additional chemical observations.

    Fragmented pairs of identical monomer definitions alternate presence: each
    atom is present in exactly one record per frame, with no chemical reaction.
    The array case is one large chain and exercises the writer alone.
    """
    import numpy as np

    from reacnetgenerator.utils import listtobytes

    rng.N = args.atoms
    rng.step = args.frames if args.case == "fragmented" else 2
    rng.atomtype = np.zeros(rng.N, dtype=np.int64)
    rng.timestep = dict(enumerate(range(rng.step)))
    rng.source_frame_counts = [rng.step]
    rng.temp1it = 2 * rng.N if args.case == "fragmented" else 1
    rng.moleculetempfilename = str(root / "records.bin")
    with open(rng.moleculetempfilename, "wb") as file:
        for number in range(rng.temp1it):
            atoms = (
                np.asarray([number // 2])
                if args.case == "fragmented"
                else np.arange(rng.N)
            )
            frames = (
                np.arange(number % 2, rng.step, 2, dtype=np.min_scalar_type(rng.step))
                if args.case == "fragmented"
                else np.arange(rng.step)
            )
            for values in (
                atoms,
                np.column_stack((atoms[:-1], atoms[1:])),
                np.ones(len(atoms) - 1, dtype=np.int64),
                frames,
            ):
                file.write(listtobytes(values))


def _audit():
    """Observe persistent column buffers without retaining their historical rows."""
    from reacnetgenerator import _timedoutput as timed

    metrics = {
        "peak_total_buffer_bytes": 0,
        "peak_column_rows": 0,
        "peak_column_bytes": 0,
        "largest_numeric_input_rows": 0,
        "flushes": 0,
    }
    current = {}
    total = 0

    def observe(column):
        nonlocal total
        total += column.bytes - current.get(id(column), 0)
        current[id(column)] = column.bytes
        metrics["peak_total_buffer_bytes"] = max(
            metrics["peak_total_buffer_bytes"], total
        )
        metrics["peak_column_rows"] = max(metrics["peak_column_rows"], len(column.rows))
        metrics["peak_column_bytes"] = max(metrics["peak_column_bytes"], column.bytes)

    original_append, original_flush = timed._Column.append, timed._Column.flush

    def append(column, value):
        original_append(column, value)
        observe(column)

    def flush(column):
        observe(column)
        if column.rows:
            metrics["flushes"] += 1
        original_flush(column)
        observe(column)

    timed._Column.append, timed._Column.flush = append, flush
    if hasattr(timed._Column, "extend_numeric"):
        original_extend = timed._Column.extend_numeric

        def extend(column, values):
            metrics["largest_numeric_input_rows"] = max(
                metrics["largest_numeric_input_rows"], len(values)
            )
            original_extend(column, values)
            observe(column)

        timed._Column.extend_numeric = extend
    return metrics


def sample(args, root):
    """Measure a production entry point, then compare text and official manifests."""
    from smiles_batches import _fixture, _outputs

    from reacnetgenerator import ReacNetGenerator, _timedoutput
    from reacnetgenerator._detect import _Detect
    from reacnetgenerator._hmmfilter import _HMMFilter
    from reacnetgenerator._path import _CollectPaths
    from reacnetgenerator.reacnetgen import _json_compatible
    from reacnetgenerator.timedoutput import semantic_manifest

    source = args.input or (
        args.state / "input.lammpstrj" if args.state else root / "input.bond"
    )
    if args.case == "small":
        _small_source(source)
    elif args.case in ("fragmented", "array"):
        source.write_text("Synthetic ordered molecule-table replay metadata.\n")
    output = root / "output"
    output.mkdir()
    rng = ReacNetGenerator(
        inputfilename=str(source),
        inputfiletype="dump" if args.case == "real" else "bond",
        atomname=["C", "H", "O"] if args.case == "real" else ["H"],
        nproc=1,
        runHMM=args.hmm,
        output_dir=output,
        timed_output=output / "timeline.h5",
    )
    cached = args.state is not None and args.stage != "full"
    if cached:
        # Reuse a baseline-prepared detection stream so stage comparisons also
        # retain IDs/order. HMM preparation in _fixture remains outside timing.
        rng, source = _fixture(args, root)
        _outputs(rng, output)
        rng.parameters = {
            key: _json_compatible(getattr(rng, key)) for key in rng.parameters
        }
    collector = None
    if args.stage != "full":
        if not cached:
            if args.case in ("fragmented", "array"):
                _synthetic_records(rng, root, args)
            else:
                _Detect.gettype(rng).detect()
            _HMMFilter(rng).filter()
        collector = _CollectPaths.getstype(rng)
        if args.stage == "writer":
            if args.case in ("fragmented", "array"):
                rng.mname = [
                    "[H]" if args.case == "fragmented" else "chain"
                ] * rng.hmmit
            else:
                collector.atomnames = collector.atomname[collector.atomtype]
                collector._printmoleculename()
                collector.returnkeys()
    metrics = _audit() if args.audit else None
    before = _usage()
    started = time.perf_counter()
    if args.stage == "full":
        rng.run()
    elif args.stage == "path":
        collector.collect()
    else:
        with _timedoutput._TimedOutputWriter(rng) as writer:
            writer.write_molecules(rng.mname)
    wall = time.perf_counter() - started
    after = _usage()
    for key in ("parent_cpu_s", "child_cpu_s"):
        after[key] -= before[key]
    # The official manifest performs full validation before hashing every field.
    manifest = semantic_manifest(rng.timed_output)
    text = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in output.iterdir()
        if p.is_file() and p.suffix != ".h5"
    }
    return {
        "case": args.case,
        "stage": args.stage,
        "hmm": args.hmm,
        "atoms": rng.N,
        "frames": rng.step,
        "records": rng.hmmit,
        "method": multiprocessing.get_start_method(),
        "wall_s": wall,
        **after,
        "cpu_s": after["parent_cpu_s"] + after["child_cpu_s"],
        "audit": metrics,
        "semantic": {"manifest": _digest(manifest), "text": text},
        "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "summary": manifest["summary"],
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "source": _timedoutput.__file__,
            "source_sha256": hashlib.sha256(
                Path(_timedoutput.__file__).read_bytes()
            ).hexdigest(),
            "dependencies": {
                name: importlib.metadata.version(name)
                for name in ("numpy", "h5py", "rdkit", "hmmlearn", "lz4")
            },
        },
    }


def main():
    """Run one isolated sample, keeping temporary intermediates private."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case", choices=("small", "real", "fragmented", "array"), default="small"
    )
    parser.add_argument("--stage", choices=("writer", "path", "full"), default="writer")
    parser.add_argument("--input", type=Path)
    parser.add_argument("--state", type=Path)
    parser.add_argument("--atoms", type=int, default=16)
    parser.add_argument("--frames", type=int, default=8192)
    parser.add_argument("--hmm", action="store_true")
    parser.add_argument("--audit", action="store_true")
    parser.add_argument(
        "--method", choices=multiprocessing.get_all_start_methods(), default="spawn"
    )
    parser.set_defaults(nproc=1, miso=0, records=0)
    args = parser.parse_args()
    if args.case == "real":
        if (args.input is None) == (args.state is None):
            parser.error("real input requires exactly one of --input or --state")
    elif args.input or args.state:
        parser.error("--input and --state apply only to real input")
    if args.atoms < 1 or args.frames < 1:
        parser.error("--atoms and --frames must be positive")
    if args.case in ("fragmented", "array") and (args.hmm or args.stage == "full"):
        parser.error("synthetic table replays do not use the parser or HMM")
    if args.case == "array" and args.stage != "writer":
        parser.error("the large-array case measures the writer only")
    multiprocessing.set_start_method(args.method, force=True)
    with tempfile.TemporaryDirectory(prefix="rng-numeric-sample-") as folder:
        previous = tempfile.tempdir
        tempfile.tempdir = folder
        try:
            result = sample(args, Path(folder))
        finally:
            tempfile.tempdir = previous
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
