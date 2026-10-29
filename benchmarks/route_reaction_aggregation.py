# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compare route/reaction aggregation against an independently installed base.

Each sample runs in a fresh process. Timings include real worker startup, result
reduction, and output writing; deterministic input preparation is excluded.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import multiprocessing
import os
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

CASES = {
    "route-small": ("route", "repeat", 8, 32, 1),
    "route-repeat": ("route", "repeat", 32, 32768, 1),
    "route-repeat-pool": ("route", "repeat", 32, 32768, 2),
    "route-unique": ("route", "unique", 4, 32768, 1),
    "route-sparse": ("route", "sparse", 64, 32768, 1),
    "route-mixed": ("route", "mixed", 16, 16384, 2),
    "reaction-repeat": ("reaction", "repeat", 524288, 2, 1),
    "reaction-unique": ("reaction", "unique", 8192, 2, 1),
    "reaction-high-degree": ("reaction", "high-degree", 8192, 2, 1),
}


def _digest(value):
    return hashlib.sha256(value).hexdigest()


def _usage():
    own = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    scale = 1 if sys.platform == "darwin" else 1024
    return {
        "cpu_seconds": own.ru_utime
        + own.ru_stime
        + children.ru_utime
        + children.ru_stime,
        "parent_peak_mib": own.ru_maxrss * scale / 2**20,
        "child_peak_mib": children.ru_maxrss * scale / 2**20,
    }


class _EvidenceSink:
    """Hash ordered participant evidence without retaining all event dictionaries."""

    def __init__(self):
        self.digest = hashlib.sha256()
        self.count = 0

    def write_events(self, events):
        for event in events:
            self.digest.update(json.dumps(event, sort_keys=True).encode())
            self.digest.update(b"\n")
            self.count += 1


def _run_sample(case, mode):
    import numpy as np

    from reacnetgenerator._matrix import _GenerateMatrix
    from reacnetgenerator._path import _CollectSMILESPaths
    from reacnetgenerator._reaction import ReactionsFinder
    from reacnetgenerator._step3state import _AtomFrameStore, _MoleculeNameTable

    stage, pattern, atoms, frames, nproc = CASES[case]
    molecule_count = atoms * frames + 2 if pattern in ("unique", "mixed") else 2
    if stage == "reaction":
        molecule_count = 10001 if pattern == "repeat" else atoms * 2
    names = _MoleculeNameTable(
        (np.arange(molecule_count) % 2).astype(np.uint8),
        np.array(["A", "B"], dtype=object),
    )
    with tempfile.TemporaryDirectory(prefix="rng-aggregation-") as temporary:
        root = Path(temporary)
        with _AtomFrameStore((atoms, frames), molecule_count, directory=root) as store:
            if stage == "route":
                for atom in range(atoms):
                    if pattern == "unique" or (
                        pattern == "mixed" and atom >= atoms // 2
                    ):
                        timeline = np.arange(frames) + atom * frames + 3
                    elif pattern == "sparse":
                        timeline = (np.arange(frames) // 8192) % 2 + 1
                    else:
                        timeline = np.arange(frames) % 2 + 1
                    store.atomeach[atom] = timeline
                store.prepare_active_transitions()
                collector = object.__new__(_CollectSMILESPaths)
                collector.N, collector.step, collector.nproc = atoms, frames, nproc
                collector.runHMM = mode
                collector.atomname, collector.atomtype = (
                    np.array(["H"]),
                    np.zeros(atoms, dtype=int),
                )
                collector.selectatoms, collector.mname = ["H"], names
                collector.atomroutefilename = str(root / "output.route")
                matrix = object.__new__(_GenerateMatrix)
                matrix.mname = names
                store.flush()
                initial = _usage()
                start = time.perf_counter()
                result = collector._printatomroute(store)
                counts = dict(matrix._getallroute(result))
                elapsed = time.perf_counter() - start
                final = _usage()
                semantic = {
                    "counts": sorted(
                        (str(a), str(b), int(n)) for (a, b), n in counts.items()
                    ),
                    "route": _digest(Path(collector.atomroutefilename).read_bytes()),
                    "active": _digest(store.active_transitions.tobytes()),
                }
            else:
                if pattern == "repeat":
                    store.atomeach[:, 0], store.atomeach[:, 1] = 10000, 10001
                else:
                    store.atomeach[:, 0] = np.arange(atoms) * 2 + 1
                    store.atomeach[:, 1] = (
                        np.arange(atoms) * 2 + 2 if pattern == "unique" else 2
                    )
                store.prepare_active_transitions()
                store.active_transitions[:] = 1
                store.flush()
                finder = object.__new__(ReactionsFinder)
                finder.step, finder.nproc, finder.mname = frames, nproc, names
                finder.printreactionevent = mode
                finder.reactionabcdfilename = str(root / "output.reactionabcd")
                finder.reactioneventfilename = str(root / "output.events.csv")
                sink = _EvidenceSink() if mode else None
                initial = _usage()
                start = time.perf_counter()
                finder.findreactions(None, None, matrix_store=store, timed_writer=sink)
                elapsed = time.perf_counter() - start
                final = _usage()
                semantic = {
                    "counts": _digest(Path(finder.reactionabcdfilename).read_bytes()),
                    "events": _digest(Path(finder.reactioneventfilename).read_bytes())
                    if mode
                    else None,
                    "evidence": sink.digest.hexdigest() if sink else None,
                    "event_count": sink.count if sink else None,
                }
    return {
        "wall_seconds": elapsed,
        "cpu_seconds": final["cpu_seconds"] - initial["cpu_seconds"],
        "parent_peak_mib": final["parent_peak_mib"],
        "child_peak_mib": final["child_peak_mib"],
        "semantic": semantic,
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "start_method": multiprocessing.get_start_method(),
            "source": sys.modules[_CollectSMILESPaths.__module__].__file__,
            "dependencies": {
                name: importlib.metadata.version(name)
                for name in (
                    "reacnetgenerator",
                    "numpy",
                    "pandas",
                    "h5py",
                    "lz4",
                    "scipy",
                    "rdkit",
                )
            },
        },
    }


def main():
    """Run alternating samples and refuse to summarize unequal outputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-python", type=Path)
    parser.add_argument("--candidate-python", type=Path, default=Path(sys.executable))
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=list(CASES))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--worker", choices=CASES, help=argparse.SUPPRESS)
    parser.add_argument(
        "--mode", type=int, choices=(0, 1), default=0, help=argparse.SUPPRESS
    )
    args = parser.parse_args()
    if args.worker:
        args.output.write_text(
            json.dumps(_run_sample(args.worker, bool(args.mode)), indent=2)
        )
        return
    if args.baseline_python is None or args.repeats < 1:
        parser.error("--baseline-python and a positive --repeats are required")
    environment = dict(
        os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1"
    )
    environment.pop("PYTHONPATH", None)
    script = str(Path(__file__).resolve())
    rows = []
    with tempfile.TemporaryDirectory(prefix="rng-aggregation-samples-") as temporary:
        for case in args.cases:
            for mode in (0, 1):
                samples = {"baseline": [], "candidate": []}
                semantic = None
                for repeat in range(args.repeats):
                    order = (
                        ("baseline", "candidate")
                        if repeat % 2 == 0
                        else ("candidate", "baseline")
                    )
                    for revision in order:
                        executable = getattr(args, f"{revision}_python").absolute()
                        output = Path(temporary) / "sample.json"
                        command = [
                            str(executable),
                            script,
                            "--worker",
                            case,
                            "--mode",
                            str(mode),
                            "--output",
                            str(output),
                        ]
                        # A neutral cwd and separate installed environments avoid
                        # accidentally importing candidate code for the baseline.
                        subprocess.run(
                            command,
                            cwd=temporary,
                            env=environment,
                            check=True,
                            stdout=subprocess.DEVNULL,
                        )
                        sample = json.loads(output.read_text())
                        if semantic is None:
                            semantic = sample["semantic"]
                        if sample["semantic"] != semantic:
                            raise RuntimeError(
                                f"Semantic mismatch: {case}, mode={mode}, {revision}"
                            )
                        samples[revision].append(sample)
                medians = {
                    revision: {
                        key: statistics.median(sample[key] for sample in values)
                        for key in (
                            "wall_seconds",
                            "cpu_seconds",
                            "parent_peak_mib",
                            "child_peak_mib",
                        )
                    }
                    for revision, values in samples.items()
                }
                environments = [values[0]["environment"] for values in samples.values()]
                if (
                    Path(environments[0]["source"]).resolve()
                    == Path(environments[1]["source"]).resolve()
                ):
                    raise RuntimeError(
                        "Baseline and candidate imported the same source"
                    )
                dependencies = [
                    {
                        key: value
                        for key, value in env["dependencies"].items()
                        if key != "reacnetgenerator"
                    }
                    for env in environments
                ]
                if dependencies[0] != dependencies[1] or any(
                    environments[0][key] != environments[1][key]
                    for key in ("python", "platform", "start_method")
                ):
                    raise RuntimeError("Baseline and candidate environments differ")
                row = {
                    "case": case,
                    "mode": mode,
                    "medians": medians,
                    "samples": samples,
                }
                rows.append(row)
                args.output.write_text(json.dumps(rows, indent=2))
                print(
                    f"{case} mode={mode}: {medians['baseline']['wall_seconds']:.4f}s -> {medians['candidate']['wall_seconds']:.4f}s",
                    flush=True,
                )


if __name__ == "__main__":
    main()
