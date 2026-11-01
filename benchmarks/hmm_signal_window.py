# SPDX-License-Identifier: LGPL-3.0-or-later
"""Measure actual HMM filtering on deterministic long and short signal records.

Use this script with separately installed base/head versions. Preparation,
input warming and normalized stream fingerprints are outside the timer.
--audit measures submitted-but-unconsumed records separately from speed runs.
--markers allows an external Linux process-tree memory sampler to delimit the
timed phase without adding threads to the process that creates worker pools.
"""

import argparse
import hashlib
import importlib.metadata
import json
import multiprocessing
import os
import platform
import resource
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path

import numpy as np

from reacnetgenerator import ReacNetGenerator, _hmmfilter, utils
from reacnetgenerator._moleculepresence import _iter_molecule_records


def warm(path):
    """Read the compressed input in bounded blocks before timing starts."""
    with path.open("rb") as stream:
        for _ in iter(lambda: stream.read(1024 * 1024), b""):
            pass


def prepare(root, args):
    """Prepare varied observation signals, retaining identifiable atom records."""
    source = root / "records.bin"
    with source.open("wb") as output:
        for index in range(args.records):
            if index % 4 == 0:
                frames = np.arange(args.frames, dtype=np.int64)
            elif index % 4 == 1:
                frames = np.arange(
                    args.frames // 4, 3 * args.frames // 4, dtype=np.int64
                )
            else:
                frames = np.arange(
                    index % 5, args.frames, 2 + index % 3, dtype=np.int64
                )
            for value in (
                np.array([index]),
                np.empty((0, 2), dtype=int),
                np.zeros(0, dtype=int),
                frames,
            ):
                output.write(utils.listtobytes(value))
    rng = ReacNetGenerator(
        inputfilename=source,
        inputfiletype="bond",
        atomname=["H"],
        nproc=args.nproc,
        runHMM=True,
        getoriginfile=args.signals == "both",
        printfiltersignal=args.signals == "both",
        output_dir=root / "output",
        p=[0.5, 0.5],
        a=[[0.999, 0.001], [0.001, 0.999]],
        b=[[0.99, 0.01], [0.01, 0.99]],
    )
    rng.N, rng.step, rng.temp1it = args.records, args.frames, args.records
    rng.atomtype = np.zeros(args.records, dtype=int)
    rng.timestep = dict(enumerate(range(args.frames)))
    rng.source_frame_counts = [args.frames]
    rng.moleculetempfilename = str(source)
    filter = _hmmfilter._HMMFilter(rng)
    filter._initHMM()
    return filter, source


@contextmanager
def signal_blocks(filename):
    """Keep one optional signal file open only while its blocks are consumed."""
    if filename:
        with open(filename, "rb") as stream:
            yield utils.read_compressed_block(stream)
    else:
        yield None


def fingerprint(filter):
    """Verify aligned definition/origin/HMM streams independent of worker order."""
    signatures = []
    with (
        open(filter.moleculetemp2filename, "rb") as definitions_file,
        signal_blocks(filter.originfilename) as origin,
        signal_blocks(filter.hmmfilename) as filtered,
    ):
        definitions = _iter_molecule_records(definitions_file)
        streams = (origin, filtered)
        for record in definitions:
            key = int(utils.bytestolist(record[0])[0])
            row = [key, hashlib.sha256(b"".join(record)).hexdigest()]
            for stream in streams:
                block = next(stream) if stream is not None else None
                row.append(
                    hashlib.sha256(utils.bytestolist(block).tobytes()).hexdigest()
                    if block
                    else None
                )
            signatures.append(row)
        for stream in streams:
            if stream is not None and next(stream, None) is not None:
                raise AssertionError("Extra signal record")
    assert len(signatures) == filter.hmmit
    assert len({row[0] for row in signatures}) == len(signatures)
    return hashlib.sha256(json.dumps(sorted(signatures)).encode()).hexdigest()


def audit():
    """Count queue admissions and consumption without modifying worker payloads."""
    counts = {"submitted": 0, "consumed": 0, "max_outstanding": 0}
    original_produce, original_submit, original_run = (
        utils.produce,
        utils._submit_pool_item,
        _hmmfilter.run_mp,
    )

    def admitted():
        """Record one input admitted by either scheduler."""
        counts["submitted"] += 1
        counts["max_outstanding"] = max(
            counts["max_outstanding"], counts["submitted"] - counts["consumed"]
        )

    def produce(*args, **kwargs):
        """Observe the legacy record producer without changing its semaphore."""
        for item in original_produce(*args, **kwargs):
            admitted()
            yield item

    def submit(*args, **kwargs):
        """Observe the bounded runner's task admission."""
        admitted()
        return original_submit(*args, **kwargs)

    def run(*args, **kwargs):
        """Count delivered records before the HMM stage consumes them."""
        for result in original_run(*args, **kwargs):
            counts["consumed"] += 1
            yield result

    utils.produce, utils._submit_pool_item, _hmmfilter.run_mp = produce, submit, run
    return counts


def main():
    """Measure one fresh process, including worker startup and output writing."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frames", type=int, default=1000000)
    parser.add_argument("--records", type=int, default=40)
    parser.add_argument("--nproc", type=int, default=4)
    parser.add_argument(
        "--method",
        choices=multiprocessing.get_all_start_methods(),
        default=multiprocessing.get_start_method(),
    )
    parser.add_argument("--signals", choices=["both", "filtered"], default="both")
    parser.add_argument("--audit", action="store_true")
    parser.add_argument("--markers", type=Path)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    if min(args.frames, args.records, args.nproc) < 1:
        parser.error("frames, records and nproc must be positive")
    multiprocessing.set_start_method(args.method, force=True)
    with tempfile.TemporaryDirectory(prefix="rng-hmm-window-") as directory:
        root = Path(directory)
        previous = tempfile.tempdir
        tempfile.tempdir = directory
        try:
            filter, source = prepare(root, args)
            warm(source)
            counts = audit() if args.audit else None
            if args.markers:
                args.markers.mkdir(parents=True, exist_ok=True)
                (args.markers / "ready").write_text(str(os.getpid()))
            start = time.perf_counter()
            filter._calhmm()
            elapsed = time.perf_counter() - start
            if args.markers:
                (args.markers / "done").touch()
            scale = 1024 * 1024 if platform.system() == "Darwin" else 1024
            result = {
                "frames": args.frames,
                "records": args.records,
                "nproc": args.nproc,
                "method": args.method,
                "signals": args.signals,
                "seconds": elapsed,
                "parent_peak_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                / scale,
                "largest_child_peak_mib": resource.getrusage(
                    resource.RUSAGE_CHILDREN
                ).ru_maxrss
                / scale,
                "retained": filter.hmmit,
                "fingerprint": fingerprint(filter),
                "audit": counts,
                "python": platform.python_version(),
                "source_sha256": hashlib.sha256(
                    Path(_hmmfilter.__file__).read_bytes()
                ).hexdigest(),
                "dependencies": {
                    n: importlib.metadata.version(n)
                    for n in ["numpy", "hmmlearn", "lz4"]
                },
            }
        finally:
            tempfile.tempdir = previous
    args.json.write_text(json.dumps(result, indent=2))
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
