# Bounded SMILES conversion batches

This change reuses conversion metadata in each worker and amortizes ordered
SMILES task transport. It does not cache names. The single-process path from
PR #2517 remains local; the shared runner, HMM, route and reaction scheduling
are unchanged.

## Scheduling and resource limits

Fork batches start at 256 records per requested worker. Spawn and forkserver
retain the existing ordered per-record path below 8,192 records per worker.
These conservative admission thresholds were checked at several workload
sizes, with two workers and additional four-worker controls. They are not a
universal optimal crossover or a new default worker count.

Normal batches contain at most 64 records, 256 KiB of serialized record input,
and a 4 MiB decoded-size estimate. The estimate uses LZ4 content lengths with
an expansion allowance; it is a planning proxy, not an RSS guarantee. Unknown
or oversized records become singleton batches. An indivisible input may
exceed the input/estimate limits. One lookahead input can coexist with a batch.

Each worker converts one record at a time, serializing independent results
without retaining a batch of decoded graphs. Results exceeding a 1 MiB prefix
stream directly to a private temporary file; the parent decodes one result at
a time. Large NumPy buffers bypass the in-memory prefix. Ordinary result/task
caps exclude the shared runner's framing and bounded serialization copies.
Per-worker atom metadata and individual decoded molecules still scale with
the trajectory. This change does not claim a hard bound on total process RSS.

The existing ordered disk-spooled runner admits at most `2 * nproc` batches.
A bounded pre-scan obtains the exact task count without retaining a full batch
index. Workers stop before the parent removes the private result directory,
including on early close, conversion/serialization errors and worker exit.

## Reproduce

Install the baseline and candidate into separate environments with identical
dependencies on each host. Run the same benchmark script using both installed
interpreters. Import and input preparation are outside the timer. The script
prints one JSON sample per fresh process, including source hashes and versions.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

# Once, using the baseline; retain this fixed ordered detection stream.
"$RNG_BASE_PYTHON" benchmarks/smiles_batches.py --prepare \
    --input /path/to/rp3.lammpstrj --state /path/to/fixed-records

"$RNG_BASE_PYTHON" benchmarks/smiles_batches.py \
    --state /path/to/fixed-records --method spawn --stage names
"$RNG_HEAD_PYTHON" benchmarks/smiles_batches.py \
    --state /path/to/fixed-records --method spawn --stage names

# Threshold, synthetic, complete PATH and HMM-on PATH examples.
"$RNG_HEAD_PYTHON" benchmarks/smiles_batches.py \
    --case repeated --records 16384 --method spawn
"$RNG_HEAD_PYTHON" benchmarks/smiles_batches.py \
    --case mixed --records 16384 --method fork
"$RNG_HEAD_PYTHON" benchmarks/smiles_batches.py \
    --case large-output --records 16384 --stage results --method fork
"$RNG_HEAD_PYTHON" benchmarks/smiles_batches.py \
    --state /path/to/fixed-records --stage path --method spawn --records 22460
"$RNG_HEAD_PYTHON" benchmarks/smiles_batches.py \
    --state /path/to/fixed-records --stage path --method fork --hmm
```

Repeat both interpreters for every configuration in AB/BA/AB/BA order, with
no concurrent tests or builds on the same host. Use `--miso 1/2`, `--nproc 4`
and each CSV configuration for the other controls. Synthetic `repeated`,
`low-repeat`, `mixed` and `large-output` cases require no external data.
Run `--audit` separately for names/results: it repeats serialization to count
application input payloads and actual result envelopes, and records spilled
bytes and peak unconsumed task count. Audit timings are not performance samples.
Initializer metadata is transferred once per worker and is excluded from the
per-task byte totals. Runner framing and filesystem bookkeeping are excluded.

The names timer covers the production name-writing entry point. The large
return timer includes incremental hashing of each returned record, identically
for both versions. The PATH timer covers naming, route/reaction processing and
HDF5 output. HMM-on PATH uses deterministic single-process HMM preparation
outside the timer; it measures PATH after HMM, not HMM speed.

All fixed-stream text hashes and official semantic manifests must match within
a configuration. Every HDF5 file is checked by the official validator. Full
DETECT-through-MATRIX preflight runs showed differing molecule/species order
from the existing unordered DETECT stage; they are excluded from these strict
comparisons and speed claims. Complete PATH is the integration measurement.
The production end-to-end regression suite is tested separately.

CPU fields are parent and terminated-child CPU deltas. RSS fields are process
high-water marks including imports/preparation, sampled before validation.
Child RSS is the largest child high-water mark, not a simultaneous sum of
worker RSS; neither RSS column measures total process-tree physical memory.

### Generate the timing CSV

The original campaign used temporary helpers `rng-pr2-campaign.py` to capture
each sample and `rng-pr2-report.py` to assemble `smiles_batches.csv`.
The portable [CSV exporter](smiles_batches_export.py) preserves that report's
field mapping and row order. It uses only the Python standard library and
converts saved samples without rerunning measurements.

Each input file contains the JSON object printed by `smiles_batches.py`, plus
the wrapper's `label`, `variant` (`base` or `head`), zero-based `round` and
`load_average` fields. The wrapper captured `os.getloadavg()` after each
benchmark subprocess exited. For example, capture one original-stream sample
from the repository root as follows; use a dedicated directory for timing JSON:

```sh
python - "$RNG_BASE_PYTHON" /path/to/fixed-records /path/to/macos-json/real-0-base.json <<'PY'
import json
import os
import subprocess
import sys
from pathlib import Path

interpreter, state, destination = sys.argv[1:]
command = [
    interpreter, "benchmarks/smiles_batches.py", "--state", state,
    "--method", "spawn", "--stage", "names", "--records", "0",
]
run = subprocess.run(command, check=True, text=True, stdout=subprocess.PIPE)
sample = json.loads(run.stdout.strip().splitlines()[-1])
sample.update(
    label="real", variant="base", round=0, record_limit=0,
    load_average=os.getloadavg(),
)
target = Path(destination)
target.parent.mkdir(parents=True, exist_ok=True)
target.write_text(json.dumps(sample, sort_keys=True) + "\n", encoding="utf-8")
PY
```

Capture the other configurations and alternating base/head rounds in the same
way, keeping the wrapper metadata consistent with the command. Then export:

```sh
python benchmarks/smiles_batches_export.py \
    --samples macos=/path/to/macos-json \
    --samples linux=/path/to/linux-json \
    --output /path/to/smiles_batches.csv
```

`record_limit` is the requested `--records` argument, while `records` is the
actual `rng.hmmit` reported by the benchmark. In the original real-stream
samples these are 0 and 5,615, respectively. The original JSON did not retain
the argument, so the exporter includes the campaign's label-to-argument mapping;
new captures should supply `record_limit` explicitly. Unknown labels require it.

`host` comes from `--samples`; `label`, `variant` and `round` come from the
wrapper. `total_cpu` sums the parent and child CPU deltas, and `load_1m`,
`load_5m` and `load_15m` unpack `load_average`. Source hashes, Python and
dependency versions come from `environment`; output hashes come from
`semantic`. Despite its name, `text_sha256` stores the JSON mapping of output
filenames to hashes. Missing worker-module hashes and inapplicable manifest or
result hashes become empty cells. Other columns retain their sample values.

The exporter handles the main timing CSV only and rejects audit samples.
Keep IPC audits and follow-up controls in their separate datasets. Re-exporting
the original captured JSON produces a CSV identical to the published file;
new measurements will have different timings, resource use and host load.

## Measurements (2026-10-31)

Baseline: `6e35522f136de132c9513cd382539b3e4dfa875e`. Production implementation:
`9ab4964e3392b575cee0a9a5cdab1b07599e7b8e`; later benchmark/test/report changes
do not alter the measured production code. All 35 baseline Python/Cython/C++/JSON package
source files were verified against the official commit on both hosts; the
36 candidate files were also checked across hosts.

macOS: Apple M4, 24 GiB, macOS 15.6.1 arm64, Python 3.12.7, spawn.
NumPy 2.5.3, RDKit 2026.3.6, hmmlearn 0.3.3, LZ4 4.4.5, h5py 3.14.0.
Linux: Intel Core i5-14400F, 16 available logical CPUs, approximately 31 GiB,
Linux 7.0.0-38-generic x86_64, Python 3.12.9, fork.
NumPy 2.5.4, RDKit 2026.9.1, hmmlearn 0.3.3, LZ4 4.4.4, h5py 3.14.0.
Both variants on each host share exactly the same dependency installation.
Comparisons are within each host; absolute values across hosts are not a
controlled comparison of operating systems or start methods.

The desktop remained active and the Linux host was shared and heavily loaded.
No CPU affinity or exclusive node was used. Four alternating fresh-process
samples per variant/configuration are retained, including timing outliers.

The real trajectory has 5 frames, 12,326 C/H/O atoms and 5,615 detected records.
Input SHA-256:
`c7d7dce618f52b8e96451b4194c67ce63c976d22380e3b6d8c936bcb17503e85`.
Larger real-record replays cycle the fixed detection stream in its original
order; duplicates are scheduling stress inputs, not additional observations.
Synthetic cases use 256 C/O atoms: short repeated chains, longer low-repeat
chains, alternating 2/8/32/128-atom chains, or two large frame arrays among
ordinary records. They are transport/workload controls, not a chemical-accuracy
or representative cross-trajectory benchmark.

All 304 timing samples matched their configured text/record hashes and, for
PATH, the official semantic manifest. All 24 separate audit samples matched
exactly the corresponding timing-sample output. Raw [timings and resource
measurements](smiles_batches.csv) and [IPC audits](smiles_batches_ipc.csv)
retain every sample. `record_limit=0` means the unmodified real stream.

Values below are medians. Ranges retain the minimum and maximum of all four
samples; negative changes mean less elapsed time. “Batch” describes the head;
baseline parallel names always use the existing per-record runner.

### macOS spawn

| Case                   | Workers | Head       | Wall seconds, base → head | Change | Base range (s) | Head range (s) |
| ---------------------- | ------: | ---------- | ------------------------: | -----: | -------------: | -------------: |
| Tiny chains            |       2 | Per-record |             2.225 → 2.204 |  -0.9% |    2.220–2.376 |    2.193–2.226 |
| 256 chains             |       2 | Per-record |             2.287 → 2.309 |  +1.0% |    2.265–2.323 |    2.305–2.353 |
| 512 chains             |       2 | Per-record |             2.382 → 2.376 |  -0.3% |    2.370–2.472 |    2.343–2.411 |
| 8,192 chains           |       2 | Per-record |             3.369 → 3.468 |  +2.9% |    3.253–3.574 |    3.292–3.740 |
| Real names             |       2 | Per-record |             4.090 → 4.112 |  +0.5% |    3.436–5.108 |    3.515–4.735 |
| 11,230 real names      |       2 | Per-record |             7.122 → 7.302 |  +2.5% |    6.876–7.613 |    7.288–7.368 |
| 22,460 real names      |       2 | Batch      |             9.182 → 2.645 | -71.2% |    8.484–9.621 |    2.379–2.716 |
| 32,768 real names      |       4 | Batch      |            13.999 → 3.720 | -73.4% |  13.734–16.260 |    3.661–3.754 |
| Repeated chains        |       2 | Batch      |             6.879 → 2.084 | -69.7% |    6.761–7.116 |    2.045–2.127 |
| Low-repeat chains      |       2 | Batch      |             7.454 → 2.869 | -61.5% |    7.429–7.540 |    2.694–2.939 |
| Mixed sizes            |       2 | Batch      |            12.666 → 5.375 | -57.6% |  12.486–13.755 |    5.336–5.526 |
| Large returns          |       2 | Batch      |             7.229 → 2.800 | -61.3% |    7.138–7.360 |    2.699–2.899 |
| Real PATH              |       2 | Per-record |             6.077 → 6.105 |  +0.5% |    5.979–6.344 |    6.065–6.308 |
| 22,460 PATH            |       2 | Batch      |             9.856 → 4.119 | -58.2% |   9.795–10.161 |    4.038–4.218 |
| Real PATH after HMM    |       2 | Per-record |             5.473 → 5.351 |  -2.2% |    5.358–5.566 |    5.283–5.569 |
| 22,460 PATH, miso=1    |       2 | Batch      |            10.558 → 4.780 | -54.7% |  10.533–10.624 |    4.773–4.911 |
| 22,460 PATH, miso=2    |       2 | Batch      |            10.687 → 4.811 | -55.0% |  10.445–10.852 |    4.694–4.988 |
| Serial real names      |       1 | Local      |             0.784 → 0.785 |  +0.1% |    0.769–0.854 |    0.768–0.791 |
| Tiny chains, 4 workers |       4 | Per-record |             5.164 → 5.158 |  -0.1% |    5.050–5.315 |    5.038–5.341 |

| Resource case     | Total CPU seconds, base → head | Parent peak MiB, base → head | Largest-child peak MiB, base → head |
| ----------------- | -----------------------------: | ---------------------------: | ----------------------------------: |
| 22,460 real names |                 19.260 → 5.061 |                260.9 → 259.4 |                        263.4 → 91.6 |
| Mixed sizes       |                22.478 → 12.820 |                261.5 → 265.5 |                        260.4 → 92.7 |
| Large returns     |                 15.748 → 6.709 |                278.4 → 275.8 |                       280.5 → 103.7 |
| 22,460 PATH       |                 19.014 → 6.357 |                273.6 → 271.0 |                        263.5 → 91.4 |

Observed one-minute host load ranged from 3.92 to 8.93.

### Linux fork

| Case                   | Workers | Head       | Wall seconds, base → head | Change | Base range (s) | Head range (s) |
| ---------------------- | ------: | ---------- | ------------------------: | -----: | -------------: | -------------: |
| Tiny chains            |       2 | Per-record |             0.066 → 0.064 |  -4.0% |    0.041–0.077 |    0.058–0.212 |
| 256 chains             |       2 | Per-record |             0.248 → 0.223 |  -9.7% |    0.217–0.257 |    0.183–0.274 |
| 512 chains             |       2 | Batch      |             0.786 → 0.103 | -86.9% |    0.321–1.228 |    0.092–0.112 |
| 8,192 chains           |       2 | Batch      |             6.520 → 0.868 | -86.7% |    6.099–7.828 |    0.827–1.666 |
| Real names             |       2 | Batch      |             5.050 → 0.649 | -87.1% |    4.952–5.101 |    0.618–0.689 |
| 11,230 real names      |       2 | Batch      |             9.692 → 1.836 | -81.1% |    9.561–9.904 |    1.089–2.131 |
| 22,460 real names      |       2 | Batch      |            20.042 → 2.842 | -85.8% |  19.939–20.576 |    2.167–3.546 |
| 32,768 real names      |       4 | Batch      |            28.943 → 3.382 | -88.3% |  27.247–29.810 |    2.622–3.934 |
| Repeated chains        |       2 | Batch      |            12.622 → 2.142 | -83.0% |  11.831–13.609 |    1.691–2.434 |
| Low-repeat chains      |       2 | Batch      |            16.242 → 3.544 | -78.2% |  15.831–16.613 |    2.755–3.927 |
| Mixed sizes            |       2 | Batch      |            25.526 → 7.345 | -71.2% |  24.991–26.581 |    6.575–7.894 |
| Large returns          |       2 | Batch      |            17.618 → 4.304 | -75.6% |  16.715–18.600 |    3.328–4.784 |
| Real PATH              |       2 | Batch      |            10.984 → 5.185 | -52.8% |   9.012–12.311 |    5.059–5.258 |
| 22,460 PATH            |       2 | Batch      |            24.128 → 6.721 | -72.1% |  22.958–26.012 |    6.442–7.488 |
| Real PATH after HMM    |       2 | Batch      |             5.994 → 1.661 | -72.3% |    5.540–6.241 |    1.589–2.649 |
| 22,460 PATH, miso=1    |       2 | Batch      |            24.715 → 8.905 | -64.0% |  23.597–26.040 |   8.521–10.805 |
| 22,460 PATH, miso=2    |       2 | Batch      |            24.930 → 7.711 | -69.1% |  22.129–25.382 |    7.168–8.678 |
| Serial real names      |       1 | Local      |             0.559 → 0.580 |  +3.8% |    0.545–1.437 |    0.563–0.594 |
| Tiny chains, 4 workers |       4 | Per-record |             0.087 → 0.092 |  +5.6% |    0.076–0.104 |    0.091–0.102 |

| Resource case     | Total CPU seconds, base → head | Parent peak MiB, base → head | Largest-child peak MiB, base → head |
| ----------------- | -----------------------------: | ---------------------------: | ----------------------------------: |
| 22,460 real names |                 27.414 → 4.460 |                232.8 → 232.3 |                       158.7 → 158.6 |
| Mixed sizes       |                35.573 → 16.751 |                234.5 → 235.6 |                       152.6 → 153.2 |
| Large returns     |                 25.071 → 8.936 |                239.8 → 238.7 |                       163.3 → 161.8 |
| 22,460 PATH       |                 31.585 → 7.939 |                242.0 → 241.8 |                       160.6 → 160.4 |

Observed one-minute host load ranged from 14.49 to 18.11.

### Follow-up Linux controls

The original four-pair run showed +5.6% (+4.9 ms) for tiny four-worker tasks
and +3.8% (+21 ms) for the serial control. Eight additional AB/BA pairs were
run sequentially after the main campaign and audit. The original samples
remain in the tables above; all 32 additional samples are retained in the
[control CSV](smiles_batches_controls.csv).

| Case      | Wall median, base → head (s) | Change | Base range (s) | Head range (s) |
| --------- | ---------------------------: | -----: | -------------: | -------------: |
| serial    |                0.562 → 0.576 |  +2.3% |    0.499–1.012 |    0.553–1.482 |
| tiny-four |                0.096 → 0.086 | -10.9% |    0.077–0.474 |    0.076–0.345 |

The tiny four-worker slowdown did not reproduce. The serial control remained
+2.3% in median, with mixed paired differences and wide ranges; it does not
establish precise equivalence or a small speedup. Neither code path uses
batching. These controls qualify the large, consistently separated batching
gains rather than establishing universal no-regression guarantees.

### Transport audit

| Host/case                | Tasks, base → head | Input MiB, base → head | Largest result KiB, base → head | Head maximum pending | Head streamed MiB |
| ------------------------ | -----------------: | ---------------------: | ------------------------------: | -------------------: | ----------------: |
| macos, Real names        |        5615 → 5615 |        822.63 → 822.63 |                     9.18 → 9.18 |                    4 |              0.00 |
| macos, 22,460 real names |        22460 → 351 |        3291.46 → 12.61 |                    9.18 → 59.42 |                    4 |              0.00 |
| macos, 32,768 real names |        32768 → 512 |        4802.04 → 18.37 |                    9.18 → 59.42 |                    8 |              0.00 |
| macos, Mixed sizes       |        16384 → 512 |         153.65 → 23.48 |                   9.14 → 100.91 |                    4 |              0.00 |
| macos, Large returns     |        16384 → 276 |         145.34 → 19.26 |                 2049.19 → 71.46 |                    4 |              4.00 |
| linux, Real names        |          5615 → 88 |          819.96 → 3.15 |                    9.18 → 54.85 |                    4 |              0.00 |
| linux, 22,460 real names |        22460 → 351 |        3279.08 → 12.61 |                    9.18 → 59.42 |                    4 |              0.00 |
| linux, 32,768 real names |        32768 → 512 |        4783.97 → 18.37 |                    9.18 → 59.42 |                    8 |              0.00 |
| linux, Mixed sizes       |        16384 → 512 |         137.71 → 23.48 |                   9.14 → 100.91 |                    4 |              0.00 |
| linux, Large returns     |        16384 → 276 |         129.40 → 19.26 |                 2049.19 → 71.46 |                    4 |              4.00 |

The large-return case intentionally includes indivisible compressed input
records larger than 256 KiB. These remain singleton tasks; the largest returned
frame arrays instead travel in two disk streams (about 2 MiB each). Normal
batch payloads remained within their caps, and the measured unconsumed task
window never exceeded four/eight for two/four workers.

These are host- and workload-specific observations. Small-path fluctuations
and the shared-host ranges are retained rather than reported as reliable
speedups. No name cache, HMM acceleration, total-RSS guarantee, or full-pipeline
speedup is inferred from these measurements. Windows and forkserver performance
were not measured.

## Validation

- Full non-GUI suite: 583 passed, 9 skipped, 6 deselected, 6 expected failures
  and 3 non-strict XPASS before two additional timeline-filter regressions.
- Expanded batching suite: 24 passed on both macOS spawn and Linux fork,
  including the two added parent-side timeline-filter cases and a decoded-size
  admission case whose compressed input fits below the wire budget.
- Full repository hooks, changed-file hooks and `ty==0.0.15` passed.
- macOS wheel and both Linux editable builds passed, including native extension
  and frontend compilation.

The earlier uncached prototype's approximately 30.2% fork PATH reduction and
28.7% small-spawn name-stage regression motivated this work and its fallback.
They are development evidence, not results for this patch. This benchmark
isolates the present scheduling change; name caching remains a separate task.
