# SPDX-License-Identifier: LGPL-3.0-or-later
"""Replay presence queries against immutable HDF5 columns, with or without caching.

Use separate fresh processes and alternate --cached and the scalar reference.
Real inputs are bounded benchmark windows, not an entire production trajectory.
The timing excludes query preparation and validation against the scalar oracle.
"""

import argparse
import hashlib
import importlib.metadata
import json
import resource
import sys
import tempfile
import time
from pathlib import Path

import h5py
import numpy as np

from reacnetgenerator import _timedoutputrange as ranges


class _ReadAudit:
    """Count HDF5 reads in separate audit runs, excluded from speed claims."""

    def __init__(self, dataset, counts):
        self.dataset = dataset
        self.counts = counts

    def __len__(self):
        return len(self.dataset)

    def __getitem__(self, key):
        self.counts["slices" if isinstance(key, slice) else "scalars"] += 1
        return self.dataset[key]


def queries_from_timeline(file):
    """Recover exactly the participant checks made by the timeline writer."""
    group = file["transition_evidence"]
    offsets = group["participant_offsets"][:]
    ids = group["participant_molecule_id"][:]
    sides = group["participant_side"][:]
    return [
        (int(ids[i]), int(transition) + int(sides[i]))
        for row, transition in enumerate(group["transition"][:])
        for i in range(int(offsets[row]), int(offsets[row + 1]))
    ]


def synthetic(path, case):
    """Prepare small, repeated and low-locality controls outside the timer."""
    count = 7 if case == "small" else 262147
    with h5py.File(path, "w") as file:
        group = file.create_group("molecule_ranges")
        group.create_dataset("molecule_id", data=np.arange(1, count + 1), dtype="<i8")
        group.create_dataset("start_frame", data=np.zeros(count), dtype="<i8")
        group.create_dataset("end_frame", data=np.full(count, 5), dtype="<i8")
    random = np.random.default_rng(2521)
    ids = random.integers(0, count + 2, size=1000)
    if case == "repeated":
        ids[:] = count // 2
    if case == "small":
        ids = ids[:7]
    return [
        (int(mid), int(frame))
        for mid, frame in zip(ids, random.integers(-1, 7, len(ids)))
    ]


def main():
    """Measure one isolated sample and reject any changed presence result."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path)
    parser.add_argument(
        "--case", choices=["small", "repeated", "random"], default="repeated"
    )
    parser.add_argument("--cached", action="store_true")
    parser.add_argument("--audit", action="store_true")
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as folder:
        path = args.input or Path(folder) / "ranges.h5"
        queries = None if args.input else synthetic(path, args.case)
        with h5py.File(path, "r") as file:
            if queries is None:
                queries = queries_from_timeline(file)
            group = file["molecule_ranges"]
            columns = [
                group[name] for name in ("molecule_id", "start_frame", "end_frame")
            ]
            read_counts = {"scalars": 0, "slices": 0}
            timed_columns = (
                [_ReadAudit(column, read_counts) for column in columns]
                if args.audit
                else columns
            )
            start = time.perf_counter()
            if args.cached:
                lookup = ranges._MoleculeRangeLookup(*timed_columns)
                result = [lookup.contains(mid, frame) for mid, frame in queries]
            else:
                result = [
                    ranges.molecule_present(mid, frame, *timed_columns)
                    for mid, frame in queries
                ]
            elapsed = time.perf_counter() - start
            expected = [
                ranges.molecule_present(mid, frame, *columns) for mid, frame in queries
            ]
            if result != expected:
                raise AssertionError("Cached presence differs from scalar reference")
            peak_rss_kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            if sys.platform == "darwin":
                # Darwin reports bytes; Linux already reports KiB.
                peak_rss_kib /= 1024
            output = {
                "case": args.case if args.input is None else args.input.name,
                "cached": args.cached,
                "queries": len(queries),
                "rows": len(columns[0]),
                "seconds": elapsed,
                "audit": read_counts if args.audit else None,
                "results_sha256": hashlib.sha256(bytes(result)).hexdigest(),
                "peak_rss_KiB": peak_rss_kib,
                "cached_array_bytes": sum(
                    block.nbytes for block in lookup.columns[0].cache.values()
                )
                if args.cached
                else 0,
                "source_sha256": hashlib.sha256(
                    Path(ranges.__file__).read_bytes()
                ).hexdigest(),
                "versions": {
                    name: importlib.metadata.version(name) for name in ("h5py", "numpy")
                },
            }

    args.json.write_text(json.dumps(output, indent=2))
    print(json.dumps(output), flush=True)


if __name__ == "__main__":
    main()
