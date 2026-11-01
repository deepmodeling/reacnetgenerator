# HMM scheduling for long signals

The legacy runner groups 100 HMM records into each pool task, retaining large
input/result batches. For at least 100,000 analyzed frames this change uses at
most `min(nproc, ceil(record_count / 100))` effective workers: the maximum number
of concurrent predictions the legacy chunks could have supplied. One effective
worker runs in the parent; the requested `nproc` remains available to later stages.

For larger fork/spawn pools, submit individual records with at most
`2 * effective_workers` submitted-but-unconsumed records. Initialize each worker
with only the actual model, frame count and two output flags. Forkserver retains
100-record chunks because measured per-record transport overhead outweighs the
benefit on the 100k-frame fixture; its chunk callable also carries only prediction
state. The full RNG context is excluded from both transport paths.

Serial and worker paths share one conversion function. Passing the actual model
preserves changed prediction parameters and algorithm selection. Short signals
and no-HMM processing retain their scheduling. Results remain unordered; each
molecule definition stays together with its origin and filtered signal.

The admission bound applies to fork/spawn pools, not forkserver chunks or total
RSS. A prediction's working arrays still scale with frame count. The worker cap
avoids activating more simultaneous prediction arrays than the old chunks.

## Performance measurements at f840e83

Baseline: official `f213256926e113b8500befcea1ea6c76d784cadb`. Measured candidate `f840e83`
`_hmmfilter.py` SHA-256: `b5f9a61ac28747f9fd95010e304fca337bf8bcfa59b829dce4c61c447ccec6be`.
Both use frozen source overlays in one Python 3.12.12 environment, NumPy 2.5.4,
hmmlearn 0.3.3 and LZ4 4.4.5, with the same native extension. Git/file hashes
identify code rather than the environment's installed version label. All final
runs are sequential on one Linux AMD EPYC 7K62 host, CPUs 0–31, native math
threads limited to one; tests and checks finish before timing begins.

The fixture constructs real `ReacNetGenerator` / `_HMMFilter` objects with a
frame-sized timestep dictionary. Synthetic signals mix dense intervals and
periodic sparse indices; they are transport/filtering tests, not additional
chemical trajectories. Preparation, explicit input pre-reading and fingerprints
are outside the timer. `_calhmm` includes pool startup, prediction, decoding,
compression, writing and cleanup. Model probabilities match across variants.

Each row below uses three alternating fresh-process AB/BA/AB pairs. Medians in
seconds:

| Workload                                          |  Base | Candidate | Change |
| ------------------------------------------------- | ----: | --------: | -----: |
| fork, 100k frames, 400 records, 4 requested       | 0.963 |     0.859 | -10.7% |
| spawn, 100k frames, 400 records, 4 requested      | 2.799 |     2.168 | -22.5% |
| forkserver, 100k frames, 400 records, 4 requested | 2.817 |     2.235 | -20.7% |
| fork, 1M frames, 400 records, 32 requested        | 9.508 |     7.250 | -23.7% |
| fork, 10k frames, 40 records, 4 requested         | 0.074 |     0.082 | +10.6% |

The short control changes by +7.8 ms (+10.6%) despite retaining
legacy scheduling. No universal speedup is claimed. Below-threshold, boundary,
single-process, low-record-count, filtered-only and 1M-frame forkserver controls
use one pair each; they check correctness and do not establish speed claims.
[All 42 final timing samples](hmm_signal_window_timings.csv) retain controls and
source hashes. Windows/macOS performance has not been measured.
Forkserver total times vary across independent reviews; these fixture samples
do not establish a stable speedup for that start method. The subsequent iterator
cleanup below was not used to regenerate these timing or memory samples.

## Review-driven correction

Previous PR head `6a5e61e` passed a bound converter to each worker initializer.
`SharedRNGData.rng` made that converter carry the full RNG context, including
frame metadata. Earlier single-pair actual-pool samples already slowed down;
classifying them only as compatibility checks missed a regression that required
investigation. This issue is separate from the old CodSpeed hardware warning.

We repeated the 100k-frame / 400-record / four-worker fixture in three cyclic
orders: base/previous/trial, trial/base/previous, previous/trial/base. The initial
trial transferred only prediction state but still used per-record forkserver
transport. Medians in seconds:

| Method     |  Base | Previous PR | Initial prediction-only trial |
| ---------- | ----: | ----------: | ----------------------------: |
| fork       | 0.962 |       0.842 |                         0.819 |
| spawn      | 2.767 |       7.599 |                         2.115 |
| forkserver | 2.866 |       9.289 |                         4.049 |

This trial fixed spawn startup but left forkserver slower than the baseline.
Separate phase diagnostics put its remaining overhead primarily in task
consumption, so the final implementation keeps forkserver chunking. No global
multiprocessing defaults or shared runner internals are changed.
[All 27 review-stage samples](hmm_worker_state_review.csv) retain this rejected
intermediate design (`candidate` in that file) separately from final samples.

Using the real RNG constructor, the actual fork/spawn initializer envelope
shrinks from 744,447 to 809 bytes at 100k frames and from
9,747,480 to 809 bytes at 1M frames. This measures only
serialized initialization state, not total memory or record payloads.

## Memory and queue audits

Separate Linux audits sample summed process-tree PSS and USS during filtering,
with a 20 ms polling delay plus procfs read time. They also count admissions.
Each configuration has one pair; sampled peaks can miss brief allocations and
are not hard memory bounds. Audit timings are excluded from speed tables. PSS
avoids counting shared pages in full in every child. Ordinary per-process RSS
columns include preparation and are not process-tree totals.

| Records (1M frames, fork, 32 requested) | Base tree PSS MiB | Candidate tree PSS MiB | Outstanding records, base / candidate |
| --------------------------------------- | ----------------: | ---------------------: | ------------------------------------: |
| 40                                      |             869.9 |                  535.7 |                                40 / 0 |
| 400                                     |            2768.7 |                 1800.4 |                               400 / 8 |

Zero queued records means local execution, not skipped analysis. The 400-record
fork candidate uses four workers and at most eight outstanding records.
[Audit data](hmm_signal_window_memory.csv) include USS, process and sample counts.
The standalone [memory sampler](hmm_signal_memory.py) can repeat these audits.

## Correctness and scope

704 non-GUI tests pass: 6 GUI deselections, 6 xfailed and 3 non-strict XPASS.
The 29 scheduling cases cover boundaries, the forkserver chunk fallback, actual
fork/spawn/forkserver output alignment, output flags, forced reordering, worker
errors/exits, failed consumption, and payload sizes. Model-state tests preserve
modified probabilities and algorithm selection while the parent RNG contains
an unrelated unpicklable callback. Filtering calculations and retention rules
are preserved.

All 42 final timing samples and four audits have identical aligned definition
and signal fingerprints within each configuration. Streams are joined before
normalization by atom ID, so separate sorting cannot hide misalignment.

At `f840e83`, all 14 integration outputs pass native validation; seven base/candidate
pairs preserve exact labels, graphs, ranges, populations, counts and transition
evidence after internal ID/table-order normalization. Long prepared HMM→MATRIX
cases at 100k/1M frames cover `miso=0/1/2` with 240/102/240 records; an additional
100k-frame / 240-record case runs the full pipeline with forkserver. These are
synthetic plumbing checks. Two real 5,001-frame, 2,600-atom windows run DETECT
through MATRIX with HMM, plus a 101-frame no-HMM control. The real windows are
below the threshold and do not establish a production long-trajectory speedup.
[Integration samples](hmm_signal_window_integration.csv) are correctness timings,
not repeated speed measurements.

## Iterator cleanup and run metadata

All HMM result iterators are closed explicitly, including the forkserver chunk
fallback, short-signal pools and the local progress iterator. If writing a
result fails, closing `run_mp` terminates and joins its pool before the error
escapes. Cleanup does not depend on releasing an exception traceback or garbage
collecting the generator.

Four regression cases retain the write-error traceback while checking that no
new workers remain: short/long signals with spawn/forkserver. Before the fix,
both short-signal cases and the long forkserver case leaked workers. The existing
bounded long-signal spawn path already closed its iterator. All four now pass.
Two additional 100k-frame / 240-record HMM→MATRIX comparisons against `f840e83`
cover spawn with `miso=2` and forkserver with `miso=1`; all four outputs pass
native validation and the two pairs preserve exact semantic evidence.

The pipeline benchmark now records `method` and `miso` alongside `nproc`, so
its JSON result identifies the actual scheduling and MISO configuration. This
follow-up changes cleanup and metadata only; the performance tables above retain
their measured source revision.

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
measurements use `--method spawn` or `--method forkserver` on Linux; Windows/macOS
performance has not been measured. Use `hmm_signal_pipeline.py --method forkserver` to check that transport through
MATRIX too. For 100k-frame integration use 80 pairs and
`--miso 1` or `--miso 2`. Real-pipeline controls use `equal_z_candidates.py --mode pipeline --hmm --nproc 32 --miso 1 --atom-types C Cl H O`, with the input and fresh
output paths specified as in the [fixed-window protocol](isomer_candidates.md).
Pre-read the trajectory before each invocation; omit `--hmm` for the pilot
no-HMM control.

## Same-host diagnostic for the automated benchmark warning

This section records a historical diagnostic at PR head `6a5e61e`. These
kernel-only samples do not exercise pool initialization and cannot rule out
the spawn/forkserver regression described above.

The initial CodSpeed report for this PR reported a 28.19% aggregate degradation,
while flagging all 12 significantly changed benchmarks as using different
hardware and warning about Walltime outside Macro Runners. This does not identify
which changes are caused by code. The functions timed by those tests do not call
the new HMM scheduler; detection, route/reaction kernels and CLI/import paths are
unchanged at that head, as was `_getoriginandhmm` itself.

We replayed those 12 checked-in test bodies on CPU 0 with the same dependency
and source environments, in five alternating fresh-process pairs. The local
harness warms each callable, selects iterations targeting 80 ms, and records
seven timing blocks per process; CLI/import tests use one subprocess per block.
Preparation is outside numerical timers. This is a local wall-clock diagnostic,
not a CodSpeed score or a recreation of its runner.

| Test body            | Median elapsed-time change |
| -------------------- | -------------------------: |
| detect_param0        |                    -0.636% |
| detect_param1        |                    +0.747% |
| detect_param2        |                    +0.477% |
| hmm_param0           |                    -0.048% |
| hmm_param2           |                    +0.009% |
| hmm_param3           |                    +0.133% |
| atom_route           |                    -0.507% |
| large_periodic_frame |                    -0.502% |
| lammps_parser        |                    +0.494% |
| transition_graph     |                    +0.450% |
| cli_help             |                    +0.024% |
| module_import        |                    -0.005% |

Changes range from −0.64% to +0.75%; the reported aggregate slowdown was not
reproduced on this host. This does not prove identical performance on every
runner. [All 120 test-body samples](hmm_kernel_audit.csv) retain the five pairs.
Run the portable harness with both installations and the same CPU affinity,
using the environment settings above, and alternate their order:

```sh
taskset -c 0 "$RNG_BASE_PYTHON" benchmarks/hmm_kernel_audit.py \
    --variant base --json base-kernels.json
taskset -c 0 "$RNG_HEAD_PYTHON" benchmarks/hmm_kernel_audit.py \
    --variant candidate --json head-kernels.json
```
