# Bounded HMM scheduling for long signals

The legacy runner groups 100 HMM records into each pool task. A chunk retains
many compressed definitions/results, and small inputs still create the full
requested pool. For at least 100,000 analyzed frames this change uses at most
`min(nproc, ceil(record_count / 100))` effective workers: the maximum number of
concurrent predictions the legacy chunks could have supplied. One effective
worker runs in the parent. Larger pools submit individual records with at most
`2 * effective_workers` submitted-but-unconsumed records.

The original initialized converter is installed once per worker. Its RNG
context, including frame metadata, is not serialized again with every record.
This preserves the model and filter flags while amortizing context transfer.
The requested `nproc` remains unchanged for later stages. Short signals and
no-HMM processing retain their scheduling. Completed results remain unordered;
each molecule definition stays together with its origin and filtered signal.

This bounds record admission, not total RSS. One prediction's working arrays
and the initialized converter still scale with trajectory size. The worker cap
is necessary: merely reducing chunk size can activate many more HMM working
arrays simultaneously when there are few records.

## Method and measurements

Baseline: official `f213256926e113b8500befcea1ea6c76d784cadb`. Measured candidate
`_hmmfilter.py` SHA-256: `69b3ffe9198aa0bbcfde4a3c9a70c5f841ded0bd4e6678d2b74d2c985fdaead5`.
Both use frozen source overlays in one Python 3.12.12 environment, NumPy 2.5.4,
hmmlearn 0.3.3 and LZ4 4.4.5, with the same native extension. Git/file hashes
identify the implementation rather than the environment's installed version label.
Runs are sequential on one Linux AMD EPYC 7K62 host, CPUs 0–31, native math
threads limited to one, with no concurrent tests or builds.

The fixture constructs real `ReacNetGenerator` / `_HMMFilter` instances, including
a frame-sized timestep dictionary. Synthetic observations mix dense intervals
and periodic sparse indices. They exercise transport and HMM filtering, not
additional chemical trajectories. Model preparation, compressed input creation,
explicit pre-reading and output fingerprints are outside the timer; `_calhmm`
includes pool startup, decoding, prediction, compression, writing and cleanup.
Speed samples use fresh processes and three alternating AB/BA/AB pairs.

Medians in seconds:

| Workload                                          |  Base | Candidate | Change |
| ------------------------------------------------- | ----: | --------: | -----: |
| 100k frames, 40 records, 4 requested              | 0.291 |     0.228 | -21.6% |
| 1M frames, 40 records, 4 requested                | 2.726 |     2.110 | -22.6% |
| 1M frames, 400 records, 4 requested               | 8.424 |     6.854 | -18.6% |
| 1M frames, 40 records, 32 requested               | 2.936 |     2.080 | -29.2% |
| 1M frames, 400 records, 32 requested              | 8.564 |     6.944 | -18.9% |
| 1M frames, 40 records, 4 requested, filtered-only | 2.700 |     2.039 | -24.5% |
| 10k frames, 40 records, 4 requested               | 0.075 |     0.084 | +12.4% |

The short control adds about 9 ms in these samples despite retaining the legacy
scheduling path. This change does not claim a universal speedup.

The below-threshold, explicit single-process, spawn and forkserver controls use
one pair each and are compatibility checks, not speed claims. Spawn/forkserver
include both 40-record local cases and 400-record actual-pool cases at 100k
frames. [All 54 speed samples](hmm_signal_window_timings.csv) retain negative
controls and source hashes. The two-state model's probabilities are fixed and
identical across variants; no parameter fitting is timed.

Separate Linux audits sample summed process-tree PSS and USS during filtering,
with a 20 ms polling delay plus procfs read time. They also count pool admissions.
Each configuration has one base/candidate pair; sampled peaks can miss brief
allocations and do not establish a universal memory reduction. These timings
are excluded from the speed table. PSS avoids summing shared pages as full RSS
in every child. Ordinary per-process RSS columns include preparation and are
not process-tree totals.

| Workload                             | Base tree PSS MiB | Candidate tree PSS MiB | Outstanding pool records, base / candidate |
| ------------------------------------ | ----------------: | ---------------------: | -----------------------------------------: |
| 1M frames, 40 records, 4 requested   |             872.3 |                  535.8 |                                     40 / 0 |
| 1M frames, 400 records, 4 requested  |            2698.5 |                 1787.0 |                                    400 / 8 |
| 1M frames, 40 records, 32 requested  |             924.3 |                  535.7 |                                     40 / 0 |
| 1M frames, 400 records, 32 requested |            2913.5 |                 1792.7 |                                    400 / 8 |

Zero queued records means local execution, not skipped analysis. The 400-record
candidate cases use four workers and at most eight outstanding records.
[Audit data](hmm_signal_window_memory.csv) include USS, process counts and sample
counts. The portable [memory sampler](hmm_signal_memory.py) can repeat PSS/USS
sampling independently of the optional queue instrumentation.

## Correctness and scope

692 non-GUI tests pass, with 6 GUI deselections, 6 xfailed and 3 non-strict XPASS.
Seventeen added cases cover pool/local boundaries, real fork/spawn output
alignment, all origin/filter-output combinations, forced out-of-order completion,
worker errors/exits, failed consumption, and exclusion of a large RNG context
from individual task payloads. HMM parameters and prediction code are unchanged.

All 54 speed samples and eight audits retain identical normalized definition
and signal fingerprints within each configuration. Definitions, origins and HMM
signals are read together before normalization by atom ID, so misaligned streams
cannot be hidden by independent sorting.

Twelve integration outputs pass native validation, and all six base/candidate
pairs preserve exact labels, graphs, ranges, populations, counts and transition
evidence after internal ID/table-order normalization. Three prepared long-signal
fixtures run HMM through MATRIX at 100k/1M frames with `miso=0/1/2`; they contain
240/102/240 records and exercise the actual bounded pool with full RNG metadata.
These are synthetic plumbing checks. Two real 5,001-frame, 2,600-atom windows
run DETECT through MATRIX with HMM on, plus a 101-frame no-HMM control. Their
signals are below the new threshold and do not establish a long-trajectory
production speedup. [Integration samples](hmm_signal_window_integration.csv)
contain single-pair correctness timings, not benchmark claims.

## Reproduce

Use identical dependencies and this checkout's scripts with separately installed
base/candidate interpreters. Replace every placeholder below. Run fresh output
names sequentially and alternate version order; `0-31` must be an available CPU set.
Use a writable Matplotlib cache shared by the runs to avoid temporary-cache
creation during worker imports.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0
export MPLCONFIGDIR=/path/to/writable/matplotlib-cache
mkdir -p "$MPLCONFIGDIR"
RNG_BASE_PYTHON=/path/to/base/bin/python
RNG_HEAD_PYTHON=/path/to/head/bin/python

taskset -c 0-31 "$RNG_BASE_PYTHON" benchmarks/hmm_signal_window.py \
    --frames 1000000 --records 400 --nproc 32 --method fork --json base.json
taskset -c 0-31 "$RNG_HEAD_PYTHON" benchmarks/hmm_signal_window.py \
    --frames 1000000 --records 400 --nproc 32 --method fork --json head.json

# Separate audit; do not use its wall time as a speed sample.
taskset -c 0-31 "$RNG_HEAD_PYTHON" benchmarks/hmm_signal_memory.py \
    --output head-memory.json -- "$RNG_HEAD_PYTHON" \
    benchmarks/hmm_signal_window.py --frames 1000000 --records 400 \
    --nproc 32 --method fork --audit

# Long prepared integration; use both interpreters and compare their outputs.
taskset -c 0-31 "$RNG_HEAD_PYTHON" benchmarks/hmm_signal_pipeline.py \
    --frames 1000000 --pairs 34 --nproc 4 --miso 0 \
    --output head-output --json head-pipeline.json
taskset -c 0-31 "$RNG_BASE_PYTHON" benchmarks/hmm_signal_pipeline.py \
    --frames 1000000 --pairs 34 --nproc 4 --miso 0 \
    --output base-output --json base-pipeline.json
"$RNG_HEAD_PYTHON" benchmarks/equal_z_compare.py base-output head-output
```

Use 10,000 / 99,999 / 100,000 frames for the short and boundary controls, 40
records for low-task cases, and `--signals filtered` to suppress origin and
zero-filter outputs. Requested worker counts are 1, 4 or 32. Start-method
controls use `--method spawn` or `--method forkserver` on Linux; Windows/macOS
performance has not been measured. For 100k-frame integration use 80 pairs and
`--miso 1` or `--miso 2`. Real-pipeline controls use `equal_z_candidates.py --mode pipeline --hmm --nproc 32 --miso 1 --atom-types C Cl H O`, with the input and fresh
output paths specified as in the [fixed-window protocol](isomer_candidates.md).
Pre-read the trajectory before each invocation; omit `--hmm` for the pilot
no-HMM control.
