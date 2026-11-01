# SPDX-License-Identifier: LGPL-3.0-or-later
"""Replay the 12 reported benchmark bodies with a local wall-clock harness.

This is a matched-host diagnostic, not a reproduction of CodSpeed's runner.
Fixture preparation and imports are outside the numerical kernel timers.
CLI/import tests intentionally time their own fresh subprocess startup.
"""

import argparse
import copy
import hashlib
import importlib.metadata
import importlib.util
import json
import os
import statistics
import sys
import time
from pathlib import Path

TREE = Path(__file__).resolve().parents[1]
PARAMS = json.loads((TREE / "tests/test.json").read_text())


def module(name, filename):
    """Load the same checked-in fixture and test body for either installation."""
    spec = importlib.util.spec_from_file_location(name, TREE / "tests" / filename)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main():
    """Run the reported test bodies with one pinned-environment wall-clock harness."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", default="sample")
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    os.environ["PATH"] = (
        str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", "")
    )
    from reacnetgenerator import _detect, _hmmfilter

    original = module("source_reacnetgen_tests", "test_reacnetgen.py")
    step3 = module("source_step3_tests", "test_step3_benchmark.py")
    cli = module("source_cli_tests", "test_cli.py")
    results = {}
    current = None

    def benchmark(function, *, fixed_iterations=None):
        """Warm one callable and retain seven repeated timing blocks."""
        function()
        start = time.perf_counter()
        function()
        pilot = time.perf_counter() - start
        iterations = fixed_iterations or max(
            1, min(10000, int(0.08 / max(pilot, 1e-7)))
        )
        samples = []
        for _ in range(7):
            start = time.perf_counter()
            for _ in range(iterations):
                function()
            samples.append((time.perf_counter() - start) / iterations)
        results[current] = {
            "seconds": statistics.median(samples),
            "samples": samples,
            "iterations": iterations,
        }
        print(current, results[current]["seconds"], flush=True)
        return function

    for index in (0, 1, 2):
        current = f"detect_param{index}"
        param = copy.deepcopy(PARAMS[index])
        param["rngparams"]["inputfilename"] = str(
            TREE / "tests" / param["rngparams"]["inputfilename"]
        )
        original.TestReacNetGen().test_benchmark_detect(benchmark, param)
    for index in (0, 2, 3):
        current = f"hmm_param{index}"
        param = copy.deepcopy(PARAMS[index])
        param["rngparams"]["inputfilename"] = str(
            TREE / "tests" / param["rngparams"]["inputfilename"]
        )
        original.TestReacNetGen().test_benchmark_hmm(benchmark, param)
    current = "atom_route"
    step3.test_benchmark_atom_route(benchmark)
    current = "large_periodic_frame"
    detect_tests = module("source_detect_tests", "test_detect.py")
    fixture = detect_tests.TestPeriodicOpenBabelCandidates()
    detector = detect_tests.TestPeriodicOpenBabelCandidates.detect_instance.__wrapped__(
        fixture
    )
    fixture.test_benchmark_large_periodic_frame(detector, benchmark)
    current = "lammps_parser"
    module(
        "source_parser_tests", "test_lammps_parser.py"
    ).test_benchmark_lammps_dump_parser(benchmark)
    current = "transition_graph"
    step3.test_benchmark_transition_graph(benchmark)
    for current, function in [
        ("cli_help", cli.test_cli),
        ("module_import", cli.test_bench_module_import),
    ]:
        benchmark(function, fixed_iterations=1)
    output = {
        "variant": args.variant,
        "results": results,
        "python": sys.version,
        "affinity": sorted(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity")
        else None,
        "source": _detect.__file__,
        "hmm_source_sha256": hashlib.sha256(
            Path(_hmmfilter.__file__).read_bytes()
        ).hexdigest(),
        "dependencies": {
            name: importlib.metadata.version(name)
            for name in [
                "numpy",
                "scipy",
                "openbabel",
                "h5py",
                "rdkit",
                "hmmlearn",
                "networkx",
            ]
        },
    }
    args.json.write_text(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
