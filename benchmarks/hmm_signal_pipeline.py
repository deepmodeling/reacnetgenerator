# SPDX-License-Identifier: LGPL-3.0-or-later
"""Check HMM through MATRIX on prepared long synthetic molecule observations.

Independent hydrogen pairs change once from dimers to individual atoms. This
tests data plumbing and evidence consistency, not a chemical trajectory or
DETECT performance. Preparation and native output validation are outside timing.
"""

import argparse
import json
import multiprocessing
import time
from pathlib import Path

import numpy as np

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator.timedoutput import validate_timed_output
from reacnetgenerator.utils import listtobytes


def main():
    """Run a deterministic prepared observation stream through production stages."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frames", type=int, default=100000)
    parser.add_argument("--pairs", type=int, default=20)
    parser.add_argument("--nproc", type=int, default=4)
    parser.add_argument(
        "--method",
        choices=multiprocessing.get_all_start_methods(),
        default=multiprocessing.get_start_method(),
    )
    parser.add_argument("--miso", type=int, choices=[0, 1, 2], default=1)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    if args.frames < 16 or min(args.pairs, args.nproc) < 1:
        parser.error("Use at least 16 frames, one pair and one worker")
    multiprocessing.set_start_method(args.method, force=True)
    args.output.mkdir(parents=True, exist_ok=False)
    source = args.output / "synthetic.bond"
    source.touch()
    records = args.output / "prepared.bin"
    with records.open("wb") as stream:
        for pair in range(args.pairs):
            left, right = 2 * pair, 2 * pair + 1
            boundary = args.frames // 3 + pair % max(1, args.frames // 3)
            definitions = [
                ([left, right], [(left, right)], [1], np.arange(boundary)),
                ([left], [], [], np.arange(boundary, args.frames)),
                ([right], [], [], np.arange(boundary, args.frames)),
            ]
            for definition in definitions:
                for value in definition:
                    stream.write(listtobytes(np.asarray(value, dtype=np.int64)))
    rng = ReacNetGenerator(
        inputfilename=source.resolve(),
        inputfiletype="bond",
        atomname=["H"],
        nproc=args.nproc,
        runHMM=True,
        miso=args.miso,
        output_dir=args.output.resolve(),
        timed_output=args.output.resolve() / "timeline.h5",
    )
    rng.N, rng.step, rng.temp1it = 2 * args.pairs, args.frames, 3 * args.pairs
    rng.atomtype = np.zeros(rng.N, dtype=int)
    rng.timestep = dict(enumerate(range(args.frames)))
    rng.source_frame_counts = [args.frames]
    rng.moleculetempfilename = str(records.resolve())
    with records.open("rb") as file:
        for _ in iter(lambda: file.read(1024 * 1024), b""):
            pass
    start = time.perf_counter()
    rng._process((rng.Status.HMM, rng.Status.PATH, rng.Status.MATRIX))
    elapsed = time.perf_counter() - start
    summary = validate_timed_output(rng.timed_output)
    result = {
        "seconds": elapsed,
        "frames": args.frames,
        "pairs": args.pairs,
        "records": rng.temp1it,
        "retained": rng.hmmit,
        "nproc": args.nproc,
        "native_frames": summary.frames,
        "native_events": summary.reaction_events,
    }
    args.json.write_text(json.dumps(result, indent=2))
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
