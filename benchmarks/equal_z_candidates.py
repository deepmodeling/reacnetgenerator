# SPDX-License-Identifier: LGPL-3.0-or-later
"""Measure periodic candidate dispatch with fixed, reproducible z ties.

Use the same script in fresh base/head processes, alternating execution order.
Kernel reference validation is excluded from elapsed time. Pipeline mode accepts
an external LAMMPS dump and validates the resulting HDF5 after timing.
"""

import argparse
import hashlib
import importlib.metadata
import json
import resource
import time
from pathlib import Path

import numpy as np
from ase import Atoms
from openbabel import openbabel

from reacnetgenerator import ReacNetGenerator, _detect
from reacnetgenerator._detect import _DetectLAMMPSdump
from reacnetgenerator.timedoutput import semantic_manifest


def fixture(case):
    """Construct independent pairs, safe controls, and a connected ambiguity."""
    count = 64 if case == "small" else 2048
    index = np.arange(count // 2)
    centers = np.column_stack(
        (2 + 4 * (index % 11), 2 + 4 * ((index // 11) % 11), 2 + 4 * (index // 121))
    ).astype(float)
    positions = np.repeat(centers, 2, axis=0)
    positions[1::2, 2] += 0.74
    numbers = np.ones(count, dtype=int)
    cell = np.eye(3) * 60
    if case == "unique":
        positions[:, 2] += np.arange(count) * 1e-6
    elif case == "connected":
        positions[1] = positions[0] + np.array([0.74, 0, 0])
    elif case == "rounded":
        random = np.random.default_rng(20261031)
        positions = random.uniform(1, 79, (count, 3))
        positions[:, 2] = np.round(positions[:, 2], 2)
        numbers = random.choice([1, 6, 7, 8], count)
        cell = np.eye(3) * 80
    return Atoms(numbers=numbers, positions=positions, cell=cell, pbc=True), cell


def record_digest(detector, bonds):
    """Hash exact per-molecule records without their list ordering."""
    digest = hashlib.sha256()
    for record in sorted(detector._connectmolecule(*bonds)):
        digest.update(len(record).to_bytes(8, "little"))
        digest.update(record)
    return digest.hexdigest()


def usage():
    """Return parent/child CPU and separate process RSS high-water marks."""
    own = resource.getrusage(resource.RUSAGE_SELF)
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {
        "cpu_seconds": own.ru_utime + own.ru_stime + child.ru_utime + child.ru_stime,
        "parent_peak_rss_KiB": own.ru_maxrss,
        "largest_child_peak_rss_KiB": child.ru_maxrss,
    }


def kernel(args):
    """Time dispatch after comparing its exact graphs with original Open Babel."""
    rng = ReacNetGenerator(
        inputfilename="synthetic",
        inputfiletype="dump",
        atomname=["H", "C", "N", "O"],
        nproc=1,
        max_component_atoms=0,
    )
    detector = _DetectLAMMPSdump(rng)
    atoms, cell = fixture(args.case)
    reference = record_digest(detector, detector._getbondfromopenbabel(atoms, cell))
    candidate = record_digest(detector, detector._getbondfromcrd(atoms, cell))
    if reference != candidate:
        raise AssertionError("Candidate dispatch changed exact molecule records")
    eligible = detector._getperiodicbondcandidates(atoms, cell) is not None
    before = usage()
    start = time.perf_counter()
    for _ in range(args.iterations):
        detector._getbondfromcrd(atoms, cell)
    elapsed = time.perf_counter() - start
    after = usage()
    return {
        "mode": "kernel",
        "case": args.case,
        "atoms": len(atoms),
        "iterations": args.iterations,
        "seconds": elapsed / args.iterations,
        "cpu_seconds": (after["cpu_seconds"] - before["cpu_seconds"]) / args.iterations,
        "parent_peak_rss_KiB": after["parent_peak_rss_KiB"],
        "largest_child_peak_rss_KiB": after["largest_child_peak_rss_KiB"],
        "candidate_eligible": eligible,
        "accelerated_dispatch": eligible
        and len(atoms) >= _detect._OPENBABEL_PERIODIC_NEIGHBOR_MIN_ATOMS,
        "records_sha256": reference,
    }


def pipeline(args):
    """Measure DETECT through MATRIX, then validate its published timeline."""
    if args.input is None or not args.atom_types:
        raise ValueError("pipeline requires --input and --atom-types")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    rng = ReacNetGenerator(
        inputfilename=str(args.input.resolve()),
        inputfiletype="dump",
        atomname=args.atom_types,
        nproc=args.nproc,
        runHMM=args.hmm,
        miso=args.miso,
        output_dir=output,
        timed_output=str(output / "timeline.h5"),
    )
    before = usage()
    start = time.perf_counter()
    rng.run()
    elapsed = time.perf_counter() - start
    after = usage()
    manifest = semantic_manifest(output / "timeline.h5")
    (output / "semantic_manifest.json").write_text(json.dumps(manifest, indent=2))
    return {
        "mode": "pipeline",
        "atoms": rng.N,
        "frames": rng.step,
        "nproc": args.nproc,
        "hmm": args.hmm,
        "miso": args.miso,
        "seconds": elapsed,
        "cpu_seconds": after["cpu_seconds"] - before["cpu_seconds"],
        "parent_peak_rss_KiB": after["parent_peak_rss_KiB"],
        "largest_child_peak_rss_KiB": after["largest_child_peak_rss_KiB"],
        "native_validation": manifest["summary"],
        "ordered_manifest_sha256": manifest["semantic_sha256"],
    }


def main():
    """Run one fresh-process measurement with source and dependency identity."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["kernel", "pipeline"], default="kernel")
    parser.add_argument(
        "--case",
        choices=["independent", "unique", "connected", "rounded", "small"],
        default="independent",
    )
    parser.add_argument("--iterations", type=int, default=3)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--atom-types", nargs="+")
    parser.add_argument("--output", type=Path, default=Path("output"))
    parser.add_argument("--nproc", type=int, default=1)
    parser.add_argument("--miso", type=int, default=1)
    parser.add_argument("--hmm", action="store_true")
    parser.add_argument("--json", type=Path, default=Path("measurement.json"))
    args = parser.parse_args()
    if args.iterations < 1:
        parser.error("--iterations must be positive")
    result = kernel(args) if args.mode == "kernel" else pipeline(args)
    result["source_sha256"] = hashlib.sha256(
        Path(_detect.__file__).read_bytes()
    ).hexdigest()
    result["native_openbabel"] = openbabel.OBReleaseVersion()
    result["dependencies"] = {
        name: importlib.metadata.version(name)
        for name in ["numpy", "scipy", "openbabel", "h5py", "rdkit"]
    }
    args.json.write_text(json.dumps(result, indent=2))
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
