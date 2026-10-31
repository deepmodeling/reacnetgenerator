# SPDX-License-Identifier: LGPL-3.0-or-later
"""Measure production name grouping on fixed real or synthetic molecule records.

Prepare bounded real windows once with --prepare, --input and --atom-types.
Then replay the same state with base/head environments. Synthetic controls use
explicit candidate labels to isolate grouping, not to model chemical SMILES.
Preparation, imports and correctness fingerprints are outside the timer.
"""

import argparse
import hashlib
import importlib.metadata
import json
import resource
import shutil
import sys
import tempfile
import time
from pathlib import Path

import networkx as nx
import numpy as np

from reacnetgenerator import ReacNetGenerator, _path
from reacnetgenerator._detect import _Detect
from reacnetgenerator._hmmfilter import _HMMFilter
from reacnetgenerator.utils import listtobytes


def prepare(args):
    """Save one immutable no-HMM detection stream for both versions to replay."""
    if args.input is None or not args.atom_types or args.state is None:
        raise ValueError("Preparation requires --input, --atom-types and --state")
    args.state.mkdir(parents=True, exist_ok=False)
    with tempfile.TemporaryDirectory() as scratch:
        previous = tempfile.tempdir
        tempfile.tempdir = scratch
        try:
            rng = ReacNetGenerator(
                inputfilename=args.input.resolve(),
                inputfiletype="dump",
                atomname=args.atom_types,
                nproc=args.nproc,
                runHMM=False,
                output_dir=Path(scratch) / "output",
            )
            _Detect.gettype(rng).detect()
            _HMMFilter(rng).filter()
            shutil.copyfile(rng.moleculetemp2filename, args.state / "molecules.bin")
            np.save(args.state / "atomnames.npy", rng.atomname[rng.atomtype])
            metadata = {"hmmit": rng.hmmit, "atomname": rng.atomname.tolist()}
            (args.state / "metadata.json").write_text(json.dumps(metadata, indent=2))
        finally:
            tempfile.tempdir = previous
    return metadata


def synthetic(folder, case):
    """Exercise cheap admissions, rejected candidates and unavoidable collisions."""
    random = np.random.default_rng(2522)
    if case == "singletons":
        graphs = [nx.path_graph(size) for size in range(2, 66)]
    elif case == "distinct":
        graphs = [
            nx.from_prufer_sequence(random.integers(0, 12, 10)) for _ in range(64)
        ]
    elif case == "collisions":
        graphs = [nx.circular_ladder_graph(3), nx.complete_bipartite_graph(3, 3)] * 16
    elif case == "repeated":
        graphs = [nx.path_graph(12) for _ in range(64)]
    else:
        graphs = [nx.path_graph(4), nx.star_graph(3), nx.path_graph(4)]
    names = {}
    offset = 0
    with (folder / "molecules.bin").open("wb") as file:
        for index, graph in enumerate(graphs):
            atoms = np.arange(len(graph), dtype=np.int64) + offset
            names[offset] = f"candidate-{index:04d}"
            pairs = np.asarray(list(graph.edges), dtype=np.int64) + offset
            levels = 1 + ((index >> np.arange(len(pairs))) & 1)
            frames = np.arange(1 + index % 3, dtype=np.int64)
            for value in (atoms, pairs, levels, frames):
                file.write(listtobytes(value))
            offset += len(graph)
    np.save(folder / "atomnames.npy", np.full(offset, "C"))
    (folder / "metadata.json").write_text(
        json.dumps({"hmmit": len(graphs), "atomname": ["C"]})
    )
    return names


def replay(args, folder):
    """Run the unchanged production naming entry point and hash exact outputs."""
    labels = None
    if args.case == "real":
        if args.state is None:
            raise ValueError("Real replays require --state")
        state = args.state
    else:
        if args.nproc != 1:
            raise ValueError("Synthetic grouping controls require --nproc 1")
        state = folder
        labels = synthetic(folder, args.case)
    metadata = json.loads((state / "metadata.json").read_text())
    collector = object.__new__(_path._CollectSMILESPaths)
    collector.atomname = np.asarray(metadata["atomname"])
    collector.atomnames = np.load(state / "atomnames.npy", allow_pickle=False)
    lookup = {name: index for index, name in enumerate(collector.atomname)}
    collector.atomtype = np.asarray([lookup[name] for name in collector.atomnames])
    collector.hmmit, collector.nproc, collector.miso = (
        metadata["hmmit"],
        args.nproc,
        args.miso,
    )
    collector.printmoleculetime = False
    collector._moleculeframefilter = collector._moleculetimestepfilter = None
    collector.moleculetemp2filename = str(state / "molecules.bin")
    collector.moleculefilename = str(folder / "names.txt")
    if labels is not None:
        collector.convertSMILES = lambda atoms, bonds: labels[int(atoms[0])]
    # Use the same explicit pre-read for both variants, outside the timer.
    with (state / "molecules.bin").open("rb") as file:
        for _ in iter(lambda: file.read(1024 * 1024), b""):
            pass
    calls = 0
    original = _path.Molecule.isomorphic

    def counted(left, right):
        """Count grouping comparisons only in separate audit runs."""
        nonlocal calls
        calls += 1
        return original(left, right)

    if args.audit:
        _path.Molecule.isomorphic = counted
    try:
        start = time.perf_counter()
        collector._printmoleculename()
        elapsed = time.perf_counter() - start
    finally:
        _path.Molecule.isomorphic = original
    names = list(map(str, collector.mname))
    scale = 1024 * 1024 if sys.platform == "darwin" else 1024
    return {
        "case": args.case,
        "records": collector.hmmit,
        "nproc": args.nproc,
        "miso": args.miso,
        "seconds": elapsed,
        "grouping_vf2_calls": calls if args.audit else None,
        "unique_names": len(set(names)),
        "names_sha256": hashlib.sha256(json.dumps(names).encode()).hexdigest(),
        "definitions_sha256": hashlib.sha256(
            (folder / "names.txt").read_bytes()
        ).hexdigest(),
        "parent_peak_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / scale,
        "largest_child_peak_mib": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        / scale,
    }


def main():
    """Prepare a replay state or measure one fresh-process sample."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--input", type=Path)
    parser.add_argument("--atom-types", nargs="+")
    parser.add_argument("--state", type=Path)
    parser.add_argument(
        "--case",
        choices=["real", "small", "singletons", "repeated", "distinct", "collisions"],
        default="real",
    )
    parser.add_argument("--nproc", type=int, default=1)
    parser.add_argument("--miso", type=int, choices=[0, 1, 2], default=1)
    parser.add_argument("--audit", action="store_true")
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    if args.nproc < 1:
        parser.error("--nproc must be positive")
    if args.prepare:
        result = prepare(args)
    else:
        with tempfile.TemporaryDirectory() as folder:
            result = replay(args, Path(folder))
    result["source_sha256"] = hashlib.sha256(
        Path(_path.__file__).read_bytes()
    ).hexdigest()
    result["dependencies"] = {
        name: importlib.metadata.version(name)
        for name in ["numpy", "networkx", "rdkit", "lz4"]
    }
    args.json.write_text(json.dumps(result, indent=2))
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
