# SPDX-License-Identifier: LGPL-3.0-or-later
"""Measure filtering and matrix construction using the installed RNG version.

Run the same script in separate baseline and candidate environments, alternating
their order over at least five rounds. These synthetic occupancy records isolate
the changed stages; they are not an end-to-end trajectory benchmark. RSS includes
imports, input preparation, and mapped pages. Timings exclude fixture creation
and semantic hashing. Limit BLAS threads identically in both environments.
"""

import argparse
import hashlib
import json
import resource
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._hmmfilter import _HMMFilter
from reacnetgenerator._path import _CollectSMILESPaths
from reacnetgenerator._timedoutput import _TimedOutputWriter
from reacnetgenerator.timedoutput import iter_molecule_ranges, validate_timed_output
from reacnetgenerator.utils import bytestolist, listtobytes, read_compressed_block


def _records(case):
    if case == "small":
        return (
            8,
            16,
            [
                (np.asarray([atom]), np.arange(start, start + 8))
                for atom in range(8)
                for start in (0, 8)
            ],
        )
    step = 131072
    if case == "dense":
        return 512, step, [(np.arange(512), np.arange(step))]
    return (
        128,
        step,
        [
            (
                np.asarray([atom]),
                np.arange(3 * step // 4)
                if case == "mixed" and atom % 4 == 0
                else np.arange(atom, step, 1000),
            )
            for atom in range(128)
        ],
    )


def _key(atoms, frames):
    digest = hashlib.sha256(np.asarray(atoms, dtype="<i8").tobytes())
    digest.update(np.asarray(frames, dtype="<i8").tobytes())
    return digest.digest()


def _cpu():
    return sum(
        usage.ru_utime + usage.ru_stime
        for usage in (
            resource.getrusage(resource.RUSAGE_SELF),
            resource.getrusage(resource.RUSAGE_CHILDREN),
        )
    )


def benchmark(args):
    """Return stage costs and ID-independent semantic checksums for one run."""
    result = vars(args).copy()
    result["cpu_s"] = 0.0

    def measure(stage, function):
        cpu, wall = _cpu(), time.perf_counter()
        value = function()
        result[f"{stage}_wall_s"] = time.perf_counter() - wall
        result["cpu_s"] += _cpu() - cpu
        return value

    with tempfile.TemporaryDirectory(prefix="rng-presence-benchmark-") as temporary:
        directory = Path(temporary)
        source = directory / "input.bond"
        source.touch()
        rng = ReacNetGenerator(
            inputfilename=str(source),
            inputfiletype="bond",
            atomname=["H"],
            runHMM=args.hmm,
            nproc=args.nproc,
            getoriginfile=args.hmm,
            output_dir=directory,
            timed_output=directory / "timeline.h5" if args.timed else None,
        )
        rng.N, rng.step, records = _records(args.case)
        rng.atomtype = np.zeros(rng.N, dtype=int)
        rng.timestep = dict(enumerate(range(rng.step)))
        rng.source_frame_counts = [rng.step]
        rng.temp1it = len(records)
        definitions = directory / "molecules.bin"
        canonical = {
            _key(atoms, frames): i for i, (atoms, frames) in enumerate(records, 1)
        }
        with definitions.open("wb") as file:
            for atoms, frames in records:
                for value in (atoms, [], [], frames):
                    file.write(listtobytes(value))
        del records
        rng.moleculetempfilename = str(definitions)
        try:
            measure("filter", _HMMFilter(rng).filter)
            files = {
                Path(name)
                for name in (
                    rng.moleculetempfilename,
                    rng.moleculetemp2filename,
                    rng.originfilename,
                    rng.hmmfilename,
                )
                if name is not None
            }
            result["additional_temp_bytes"] = sum(
                p.stat().st_size for p in files if p != definitions
            )
            result["total_temp_bytes"] = sum(p.stat().st_size for p in files)
            ids = [0]
            with open(rng.moleculetemp2filename, "rb") as file:
                blocks = read_compressed_block(file)
                for record in zip(*[blocks] * 4, strict=True):
                    ids.append(
                        canonical[_key(bytestolist(record[0]), bytestolist(record[3]))]
                    )
            ids = np.asarray(ids, dtype="<i8")
            with measure("matrix", _CollectSMILESPaths(rng)._getatomeach) as store:
                matrix_hash, conflict_hash = hashlib.sha256(), hashlib.sha256()
                for start in range(0, rng.N, 16):
                    matrix_hash.update(
                        ids[store.atomeach[start : start + 16]].tobytes()
                    )
                    conflict_hash.update(
                        np.asarray(store.conflict[start : start + 16]).tobytes()
                    )
                result["matrix_digest"] = matrix_hash.hexdigest()
                result["conflict_digest"] = conflict_hash.hexdigest()

            if args.timed:

                def write_timeline():
                    with _TimedOutputWriter(rng) as writer:
                        writer.write_molecules(["[H]"] * rng.hmmit)

                measure("timeline", write_timeline)
                validate_timed_output(rng.timed_output)
                ranges = sorted(
                    (int(ids[row.molecule_id]), row.start_frame, row.end_frame)
                    for row in iter_molecule_ranges(rng.timed_output)
                )
                result["ranges_digest"] = hashlib.sha256(
                    json.dumps(ranges).encode()
                ).hexdigest()
            result["wall_s"] = sum(
                value for key, value in result.items() if key.endswith("_wall_s")
            )
            scale = 1024 * 1024 if sys.platform == "darwin" else 1024
            result["parent_peak_rss_mib"] = (
                resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / scale
            )
            result["child_peak_rss_mib"] = (
                resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / scale
            )
        finally:
            for filename in {
                rng.moleculetemp2filename,
                rng.originfilename,
                rng.hmmfilename,
            }:
                if filename:
                    Path(filename).unlink(missing_ok=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case", choices=("small", "sparse", "dense", "mixed"), required=True
    )
    parser.add_argument("--nproc", type=int, default=1)
    parser.add_argument("--hmm", action="store_true")
    parser.add_argument("--timed", action="store_true")
    print(json.dumps(benchmark(parser.parse_args()), sort_keys=True))
