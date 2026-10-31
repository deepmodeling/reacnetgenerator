# SPDX-License-Identifier: LGPL-3.0-or-later
"""Export captured SMILES timing samples to the published timing CSV schema.

Inputs are one JSON object per file from smiles_batches.py, augmented by the
sampling wrapper with label, variant, round and load_average. New captures
should also retain the requested --records value as record_limit. The mapping
below preserves that value for the original campaign, whose JSON omitted it.
Only the Python standard library is required; no benchmarks are rerun.
"""

import argparse
import csv
import json
from pathlib import Path

# Preserve the campaign's row order and requested --records values. Zero for
# the real cases means use the original stream, not convert zero records.
RECORD_LIMITS = {
    "tiny": 16,
    "fork-below": 256,
    "fork-threshold": 512,
    "spawn-below": 8192,
    "real": 0,
    "real-medium": 11230,
    "real-large": 22460,
    "real-four": 32768,
    "repeated": 16384,
    "low-repeat": 16384,
    "mixed": 16384,
    "large-output": 16384,
    "path-real": 0,
    "path-large": 22460,
    "path-hmm": 0,
    "path-miso1": 22460,
    "path-miso2": 22460,
    "serial": 0,
    "tiny-four": 16,
}
FIELDS = (
    "host",
    "label",
    "variant",
    "round",
    "case",
    "stage",
    "record_limit",
    "records",
    "nproc",
    "miso",
    "hmm",
    "method",
    "wall",
    "parent_cpu",
    "child_cpu",
    "total_cpu",
    "parent_peak_mib",
    "child_peak_mib",
    "load_1m",
    "load_5m",
    "load_15m",
    "input_sha256",
    "path_source_sha256",
    "worker_source_sha256",
    "runner_source_sha256",
    "text_sha256",
    "manifest_sha256",
    "results_sha256",
    "python",
    "dependencies",
)
DIRECT_FIELDS = (
    "label",
    "variant",
    "round",
    "case",
    "stage",
    "records",
    "nproc",
    "miso",
    "hmm",
    "method",
    "wall",
    "parent_cpu",
    "child_cpu",
    "parent_peak_mib",
    "child_peak_mib",
    "input_sha256",
)


def timing_row(sample, host):
    """Flatten one captured timing sample without deriving missing measurements."""
    if sample["audit"] is not None:
        raise ValueError("Audit samples must remain separate from timing samples")
    if "record_limit" in sample:
        requested = sample["record_limit"]
    elif sample["label"] in RECORD_LIMITS:
        requested = RECORD_LIMITS[sample["label"]]
    else:
        raise ValueError("Unknown campaign label: supply record_limit explicitly")
    if type(requested) is not int or requested < 0:
        raise ValueError("record_limit must be a non-negative integer")
    load_1m, load_5m, load_15m = sample["load_average"]
    environment = sample["environment"]
    sources = environment["sources"]
    semantic = sample["semantic"]
    return {
        **{key: sample[key] for key in DIRECT_FIELDS},
        "host": host,
        "record_limit": requested,
        "total_cpu": sample["parent_cpu"] + sample["child_cpu"],
        "load_1m": load_1m,
        "load_5m": load_5m,
        "load_15m": load_15m,
        "path_source_sha256": sources["_path.py"],
        "worker_source_sha256": sources.get("_smilesworker.py", ""),
        "runner_source_sha256": sources["utils.py"],
        # The historical column holds a filename-to-hash JSON object, not a
        # second digest. Keep its encoding identical to the published CSV.
        "text_sha256": json.dumps(semantic["text"], sort_keys=True),
        "manifest_sha256": semantic["manifest"],
        "results_sha256": semantic["results"],
        "python": environment["python"],
        "dependencies": json.dumps(environment["dependencies"], sort_keys=True),
    }


def _sample_folder(value):
    """Parse a host label and a directory of captured JSON samples."""
    host, separator, folder = value.partition("=")
    if not separator or not host or not folder:
        raise argparse.ArgumentTypeError("Expected HOST=DIRECTORY")
    return host, Path(folder)


def main():
    """Write deterministic timing rows, rejecting empty or duplicate inputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--samples",
        type=_sample_folder,
        action="append",
        required=True,
        metavar="HOST=DIRECTORY",
        help="repeat for each host's timing directory",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    seen = set()
    for host, folder in args.samples:
        files = sorted(folder.glob("*.json"))
        if not files:
            parser.error(f"No JSON samples in {folder}")
        for path in files:
            try:
                row = timing_row(json.loads(path.read_text()), host)
                key = (host, row["label"], row["round"], row["variant"])
                if key in seen:
                    raise ValueError(f"Duplicate timing sample: {key}")
                seen.add(key)
                rows.append(row)
            except (OSError, KeyError, TypeError, ValueError) as error:
                parser.error(f"{path}: {error}")
    order = {label: index for index, label in enumerate(RECORD_LIMITS)}
    rows.sort(
        key=lambda row: (
            row["host"],
            order.get(row["label"], len(order)),
            row["label"],
            row["round"],
            row["variant"],
        )
    )
    with args.output.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
