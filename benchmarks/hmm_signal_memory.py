# SPDX-License-Identifier: LGPL-3.0-or-later
"""Sample Linux process-tree PSS/USS separately from HMM speed measurements.

Pass the interpreter and hmm_signal_window.py command after --. The sampler
adds --json and --markers. Its polling delay is 20 ms plus procfs read time;
sampled peaks may miss brief allocations and are not hard memory bounds.
"""

import argparse
import json
import os
import signal
import subprocess
import tempfile
import time
from contextlib import suppress
from pathlib import Path


def tree_memory(pid):
    """Sum proportional and unique memory across the currently visible tree."""
    pending, seen = [pid], set()
    pss = uss = count = 0
    keys = {"Pss", "Private_Clean", "Private_Dirty", "Private_Hugetlb"}
    while pending:
        current = pending.pop()
        if current in seen:
            continue
        seen.add(current)
        folder = Path("/proc") / str(current)
        try:
            for path in (folder / "task").glob("*/children"):
                try:
                    pending.extend(map(int, path.read_text().split()))
                except (OSError, ValueError):
                    # A worker can exit between listing and reading its threads.
                    pass
            fields = {
                line.split(":", 1)[0]: int(line.split()[1]) * 1024
                for line in (folder / "smaps_rollup").read_text().splitlines()
                if line.split(":", 1)[0] in keys
            }
            pss += fields.get("Pss", 0)
            uss += sum(fields.get(key, 0) for key in keys if key != "Pss")
            count += 1
        except (OSError, ValueError):
            # Exit races can remove mappings; sample the remaining live tree.
            pass
    return pss, uss, count


def main():
    """Launch one isolated benchmark and observe only its marked filter phase."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=600)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command or args.timeout <= 0:
        parser.error("Supply a benchmark command and a positive timeout")
    if not Path("/proc/self/smaps_rollup").exists():
        parser.error("Linux procfs with readable smaps_rollup is required")
    with tempfile.TemporaryDirectory(prefix="rng-hmm-memory-") as temporary:
        root = Path(temporary)
        command = [
            *command,
            "--markers",
            str(root),
            "--json",
            str(root / "result.json"),
        ]
        process = subprocess.Popen(command, start_new_session=True)
        peak_pss = peak_uss = samples = max_processes = 0
        start = time.monotonic()
        try:
            while process.poll() is None:
                if time.monotonic() - start > args.timeout:
                    raise TimeoutError("HMM memory benchmark timed out")
                if (root / "ready").exists() and not (root / "done").exists():
                    pss, uss, processes = tree_memory(process.pid)
                    peak_pss, peak_uss = max(peak_pss, pss), max(peak_uss, uss)
                    max_processes = max(max_processes, processes)
                    samples += 1
                time.sleep(0.02)
            if process.returncode:
                raise RuntimeError(f"Benchmark exited with {process.returncode}")
            if not peak_pss or not samples:
                raise RuntimeError("No valid memory samples; use a longer workload")
            result = {
                "benchmark": json.loads((root / "result.json").read_text()),
                "peak_tree_pss_mib": peak_pss / 2**20,
                "peak_tree_uss_mib": peak_uss / 2**20,
                "samples": samples,
                "max_processes": max_processes,
                "polling_delay_seconds": 0.02,
            }
        finally:
            if process.poll() is None:
                with suppress(ProcessLookupError):
                    os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    with suppress(ProcessLookupError):
                        os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
    args.output.write_text(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
