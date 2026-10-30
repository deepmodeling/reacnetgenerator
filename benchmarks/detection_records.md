# Streaming compact molecule detection records

The PR is based on upstream `8b0d8f40`, including species-population aggregation
(#2515). Detection stores compact frame indices and streams compressed records
through a bounded writer. Large collections of short records under fork use
bounded parallel compression; small, dense and spawn workloads retain the serial
path. Input parsing, bond perception, detection workers and frame dispatch are
unchanged.

The original 154 performance samples against `b2d30fb4` are preserved below.
The review follow-up separately compares that fixed baseline, the published
serial PR `20aa950d`, and the repaired implementation. It retains slower
exploratory results and separates them from final-revision evidence.

## Storage and ordering

The parent accumulates observed frames in unsigned `array` buffers, starting at
the native `H` width and widening to `I` or `Q` when necessary. Width selection
uses actual `itemsize`, handles jumps over intermediate widths, and rejects
unrepresentable values instead of wrapping. It preserves duplicates and worker
arrival order. It does not sort detections or replace the timestep dictionary or
per-source frame-count mapping.

After collection, the default path compresses one NumPy view in the parent and
immediately passes its bytes to the writer. Under fork with at least two requested
workers, at least 24,576 records, and every frame buffer at most 1 KiB, the save
stage uses at most four compression workers and four submitted-but-unconsumed
batches. Batch size derives from a 64 KiB raw-data budget and a 2,048-record cap;
serialized input/result envelopes are checked against 256 KiB/1 MiB limits before
entering their queues. The exact bounds and ordering protocol are detailed in the
review follow-up below.

The writer pairs definitions and frame payloads strictly in first-seen molecule
order with a 1 MiB byte budget. A single indivisible serial-path payload may exceed
that budget. Failed compression, misalignment or writing closes and removes the
incomplete owned temporary file before exposing a filename.

This is not a constant-memory detector. All unique molecule definitions, compact
occurrence buffers and timestep metadata remain resident until saving finishes.
Array growth/widening, pickling and compression can allocate storage proportional
to one molecule's observations. The budget bounds accumulated output buffering,
not a single record or total process RSS. Public outputs retain their existing
semantics; the private temporary frame arrays now use narrower unsigned dtypes.

## Reproduction and evidence boundaries

Run [detection_records.py](detection_records.py) with separate installations using
matched dependencies:

```sh
python benchmarks/detection_records.py \
    --baseline-python /path/to/baseline/.venv/bin/python \
    --candidate-python /path/to/candidate/.venv/bin/python \
    --repeats 5 --start-method spawn \
    --output /tmp/detection-spawn.json
```

The script alternates baseline/candidate order in fresh processes, fixes BLAS
threads to one, and verifies source paths and runtime versions. Preparation and
validation are outside the timed interval. CPU includes the parent and terminated
children. Parent and reported child RSS high-water marks are retained separately;
they include imports/input preparation and must not be summed into a concurrent
process-tree memory measurement.

Two kinds of workload answer different questions:

- **Replay:** inject already-detected definitions into the real accumulation/save
  path. Baseline compression still uses its actual process pool. These timings
  exclude frame parsing and detection and must not be described as complete
  DETECT timings.
- **Detect:** write valid bond trajectories of isolated H atoms or changing H pairs,
  then invoke the actual detector. Timings include input reading, parsing, connected components,
  detection worker startup/IPC, accumulation, compression and output writing.
  They do not include HMM, path analysis, HDF5 or rendering.

For real workers, arrival order can differ between runs in the baseline itself.
The semantic digest therefore includes complete atom/bond/order definitions,
frame multisets (preserving multiplicity), atom types, timesteps and source frame
counts. Deterministic replays additionally require the same order-sensitive
digest. Internal dtype widths and compressed bytes are deliberately not required
to match; downstream text and HDF5 equivalence are validated separately.

| Case               |  Frames | Definitions/atoms | Workers | Pattern                                                               |
| ------------------ | ------: | ----------------: | ------: | --------------------------------------------------------------------- |
| replay-small       |       8 |                 4 |       1 | Every definition in every frame                                       |
| replay-dense-1/2/4 |  65,536 |               128 |   1/2/4 | Every definition in every frame                                       |
| replay-wide        | 131,072 |                64 |       2 | Dense observations crossing the 16-bit boundary                       |
| replay-sparse      |  32,768 |               256 |       2 | Two cyclically selected definitions per frame                         |
| replay-unique      |   8,192 |            32,768 |       2 | Four new definitions per frame                                        |
| detect-small-1/2   |       8 |                 3 |     1/2 | Real bond parsing, isolated atoms                                     |
| detect-dense-1/2/4 |     512 |               128 |   1/2/4 | Real bond parsing, isolated atoms                                     |
| detect-changing-2  |     512 |               256 |       2 | Seeded changing perfect matchings; 28,237 distinct H-pair definitions |

## Original serial PR measurements (`20aa950d`)

Measurements used an Apple M4 with 24 GiB RAM, macOS 15.6.1 arm64, Python
3.12.7, NumPy 2.5.3 and LZ4 4.4.5, with matched pinned runtime dependencies.
No other validation workload from this task ran during the measurements.
[detection_records.csv](detection_records.csv) retains all **154 samples**,
including slower cases. Every sample passed its semantic comparison; all replays
also passed the order-sensitive comparison.

Spawn results are medians of five samples per revision:

| Case           | Wall seconds, base → candidate | CPU seconds, base → candidate | Parent peak MiB, base → candidate |
| -------------- | -----------------------------: | ----------------------------: | --------------------------------: |
| replay-small   |                  1.438 → 0.004 |                 1.422 → 0.004 |                     254.9 → 257.5 |
| replay-dense-1 |                  2.222 → 0.331 |                 2.200 → 0.326 |                     441.0 → 303.0 |
| replay-dense-2 |                  2.293 → 0.344 |                 3.938 → 0.342 |                     427.2 → 305.4 |
| replay-dense-4 |                  2.169 → 0.343 |                 4.018 → 0.334 |                     482.2 → 301.8 |
| replay-wide    |                  2.078 → 0.408 |                 2.173 → 0.393 |                     461.4 → 353.6 |
| replay-sparse  |                  2.093 → 0.018 |                 4.077 → 0.017 |                     256.0 → 257.1 |
| replay-unique  |                  2.069 → 0.201 |                 4.006 → 0.193 |                     277.2 → 269.3 |
| detect-small-1 |                  2.999 → 1.457 |                 2.921 → 1.451 |                     254.6 → 256.6 |
| detect-small-2 |                  4.561 → 2.667 |                 7.050 → 5.213 |                     254.4 → 257.1 |
| detect-dense-1 |                  5.842 → 3.146 |                 5.461 → 3.113 |                     254.9 → 256.8 |
| detect-dense-2 |                  6.159 → 2.811 |                11.343 → 5.537 |                     254.9 → 257.1 |
| detect-dense-4 |                  7.282 → 3.982 |               21.472 → 14.705 |                     253.7 → 257.4 |

Fork results below use three samples per revision, except `detect-changing-2`,
which uses five. These are **macOS fork measurements**, not Linux validation.

| Case              | Wall seconds, base → candidate | CPU seconds, base → candidate | Parent peak MiB, base → candidate |
| ----------------- | -----------------------------: | ----------------------------: | --------------------------------: |
| replay-dense-2    |                  0.683 → 0.356 |                 0.774 → 0.350 |                     432.0 → 298.1 |
| replay-unique     |                  0.174 → 0.274 |                 0.322 → 0.259 |                     278.6 → 268.1 |
| detect-small-2    |                  0.031 → 0.023 |                 0.025 → 0.019 |                     254.3 → 257.3 |
| detect-dense-2    |                  0.166 → 0.145 |                 0.376 → 0.364 |                     255.1 → 256.5 |
| detect-changing-2 |                  0.370 → 0.395 |                 0.754 → 0.673 |                     273.1 → 267.0 |

The large dense replays reduce logical frame payloads from 64 to 16 MiB; the
width-expanding replay uses 32 MiB. Their parent RSS falls about 23–37% under
spawn. Small/sparse replay RSS and the small real-trajectory RSS do not improve;
imports dominate these samples. Complete spawn DETECT wall times decrease about
42–54%, substantially reflecting the removed second pool. This is not an
end-to-end analysis speedup claim.

Fork exposes a real tradeoff: saving 32,768 one-visit records takes about
0.174 → 0.274 seconds, despite lower CPU and RSS. In the additional real bond
trajectory, a deterministic random perfect matching changes every frame (seed
20261030); full DETECT takes 0.370 → 0.395 seconds, about 7% longer, with about
11% less CPU and 2% less parent RSS. The one-parent compression path is retained
for its simpler resource lifetime and lower IPC/CPU costs. It is not claimed to
minimize wall time for every distribution; compression-thread alternatives and
additional scheduling thresholds are not established by these measurements.

## Correctness coverage

`tests/test_detection_records.py` checks progressive writing before later
compression; preservation of first-seen molecule order and unsorted/duplicate
frames; safe width expansion through unsigned 64-bit boundaries; explicit
overflow; strict alignment; cleanup on compression/write failure; byte-bounded
buffering including oversized individual payloads; empty results; closing worker
results before the source on parent-side failure; and real bond/dump/xyz/extxyz
workers. Existing suites cover repeated/empty sources, global stride, incomplete
source boundaries, variable cells, atom IDs/types, bond perception fallbacks and
oversized-component guards.

```sh
python -m pytest -q tests/test_detection_records.py tests/test_detect.py tests/test_timedoutput.py
python -m pytest -q tests -k 'not gui' --tb=short -ra
```

Before the rebase, the independent fourth-item branch and a temporary integration
with the third-item species-population change (`fa46d5d9`) each matched
`b2d30fb4` across eight end-to-end configurations: all 65 ordinary text artifacts were byte-for-byte
identical and all six HDF5 semantic manifests matched. These exercise HMM on/off,
isomer modes 0/1/2, one/two workers, frame/timestep filters, splitting, selected
atoms, origin output and the species-only entry point. The integration also
passed 129 detection-record, species-population, timed-output and aromatic-species
tests. The integration is a validation-only checkout, not a branch dependency.

Before the rebase, the installed Python 3.12.7 tox environment completed the full
non-GUI suite with **490 passed, 8 skipped, 6 deselected, 6 xfailed and 3 xpassed** (22 warnings;
322.45 seconds). The eight skips require fork-based initializer monkeypatching;
the six deselections exclude GUI cases. The XFAIL/XPASS cases use the repository's
existing expected-failure fixture. Ruff lint/format, Velin and the repository's
configured production-source ty check also passed.

## Validation after the upstream refresh

Initial PR revision `20aa950d` was rebased onto official `master` at
`8b0d8f402d3ff5a41f8421a1e48cd470220fd653` (#2515). The detection implementation,
regression tests, benchmark runner and 154 historical CSV samples are unchanged
from the pre-rebase commit `022c0f5e`.

Independent baseline/candidate installations use Python 3.12.7 and the same 75
dependency versions. Both packages were rebuilt with cached CMake 4.4.3; this
build-only constraint avoided an unavailable CMake 4.4.4 download. Eight complete
baseline/candidate runs again produced 65 byte-identical generated text files and
six identical HDF5 semantic manifests, covering the configurations listed above.
The baseline now includes the merged species-population implementation, so this
comparison checks the actual upstream integration.

At that initial PR revision, the macOS fork regression run passed **255 tests**, including detection,
record streaming, utility/worker initialization, indexed Step 3, species populations,
timed output and aromatic species.

At that initial PR revision, the full non-GUI suite passed **519 tests**, with **8 skipped,
6 deselected, 6 xfailed and 3 xpassed** (22 warnings). The eight fork-only skips
were exercised by the separate fork run; the six deselections exclude GUI tests.
The expected-failure fixture accounts for the XFAIL/XPASS results. All applicable
pre-commit hooks, the configured `ty==0.0.15` production-source check and
`git diff --check` passed. Validation ran locally on macOS; native Linux,
Windows and GUI testing were not performed locally. The historical performance
measurements above were retained without rerunning or relabeling them.

## Review follow-up: bounded batch compression

The reviewed serial revision removed parallel saving unconditionally. The repair
amortizes scheduling/IPC across short records while retaining a serial fallback
where process startup is unlikely to pay for itself. A plain bounded parallel
prototype helped the targeted cases but added CPU versus the serial PR. Profiling
then found exact duplicate frame histories inside individual batches: 75% of
records in `replay-unique` and about 30% in one `detect-changing-2` run.
Distinct molecule definitions may share a frame history; they remain distinct
records. The worker reuses immutable compressed bytes only for the exact
`(array typecode, raw frame bytes)` key. Order, duplicate visits and dtype are part
of that key. Each input still produces one output, and the cache dies with its
batch. This removes repeated serialization and lets the result envelope reuse
identical byte objects without changing the file's independent record format.

### Explicit memory and ordering limits

| Component                      | Limit                                                                                         |
| ------------------------------ | --------------------------------------------------------------------------------------------- |
| Raw frame data per task        | 64 KiB, calculated from the largest record                                                    |
| Records per task               | 2,048, also bounding metadata and empty records                                               |
| Serialized input per task      | 256 KiB, checked before submission                                                            |
| Serialized result per task     | 1 MiB, checked before returning to the result queue                                           |
| Submitted but unconsumed tasks | 4, independent of total trajectory length                                                     |
| Early ordered results on disk  | At most 3 batches; each contains at most a 1 MiB result envelope plus fixed protocol metadata |
| Batch-local history cache      | At most 2,048 keys, 64 KiB raw key data and 1 MiB unique compressed values                    |
| Output write buffer            | 1 MiB; the existing indivisible-record exception remains on the serial path                   |

The pool uses `run_mp(max_inflight=4, disk_ordered=True, unordered=False)` with
an exact batch count. Its submission window advances only after the oldest batch
has been consumed. A delayed first batch therefore cannot accumulate an unbounded
reorder cache or silently admit more work. Early files are removed on consumption
or failure. Explicitly closing the compression iterator cancels worker results
before leaving a failed write. No list of all compressed trajectory results is
restored.

These are payload limits, not a process-RSS ceiling. `run_mp` adds fixed envelope
metadata and transient serialization copies. Each decoded batch has at most 2,048
objects; dictionary/list/array metadata is additionally bounded by that count.
Workers can hold their current batch/cache while results wait for consumption.
All unique definitions, compact occurrence buffers and timestep metadata remain
input-dependent resident state, and fork processes inherit the parent image.
A collection containing any frame buffer above 1 KiB uses the lazy serial path.

### Focused three-way comparison

Use three separate installations with identical dependencies:

```sh
python benchmarks/detection_records.py \
    --baseline-python /path/to/b2d30fb4/.venv/bin/python \
    --previous-python /path/to/20aa950d/.venv/bin/python \
    --candidate-python /path/to/fix/.venv/bin/python \
    --repeats 6 --start-method fork \
    --cases replay-unique detect-changing-2 \
    --output /tmp/detection-review.json
```

All six revision-order permutations are used once per six-round block, balancing
each revision in each execution position. Every sample starts a fresh process
with single-threaded BLAS. Python 3.12.7, NumPy 2.5.3, LZ4 4.4.5 and all 75 installed
dependency versions match. The baseline and previous source hashes were checked
against Git revisions `b2d30fb4` and `20aa950d`; the final source SHA-256 is
`f61cf806bbc34c31a75e437ee993751ea4d919e69ecd586a3d1e11ee67fd28ae`.

The high-churn timing was noisy, so a second complete six-round block was added
for that case. Both blocks are retained; no samples were discarded. The table
uses six samples per revision for replay and all twelve for high churn. Columns
show **fixed baseline → serial PR → repair**:

| Case              | Samples/revision |                   Wall s |                    CPU s |       Parent peak MiB |        Child peak MiB |
| ----------------- | ---------------: | -----------------------: | -----------------------: | --------------------: | --------------------: |
| replay-unique     |                6 | 0.1668 → 0.1934 → 0.1263 | 0.2822 → 0.1761 → 0.1426 | 275.7 → 265.0 → 267.7 |     24.3 → 1.0 → 21.1 |
| detect-changing-2 |               12 | 0.4046 → 0.4115 → 0.3931 | 0.7991 → 0.7199 → 0.7747 | 272.0 → 264.2 → 267.2 | 153.2 → 152.9 → 153.1 |

Replay median wall time is 24.2% below the fixed baseline and 34.7% below the
serial PR. Every repaired replay sample (0.1181–0.1314 s) is faster than the
fastest fixed-baseline sample (0.1490 s). High-churn medians are 2.8% below the
fixed baseline and 4.5% below the serial PR; each six-round block also has a
lower repaired median than its matched baseline. Its distributions overlap
substantially (baseline 0.3568–0.7902 s; repair 0.3379–0.5129 s), so this is an
observed removal of the median regression, not evidence of a uniform or
statistically established high-churn speedup. Non-task desktop services were
active during these measurements; no other validation workload from this task
ran concurrently.

Tradeoffs remain visible: parent peak RSS is about 2.7/3.0 MiB above the serial PR
for replay/high churn, but below the fixed baseline. High-churn CPU is 7.6% above
the serial PR and 3.1% below the fixed baseline. Parent and child high-water marks
must not be summed into concurrent process-tree memory. These macOS fork results
do not establish Linux, full-analysis or production-trajectory speedups.

The final 24,576-record threshold control uses three alternating serial-PR/repair
pairs, with matching results and medians 0.2322 → 0.1075 s. Earlier non-memoizing
16,384- and 24,576-record probes did not show a useful wall-time gain and increased
CPU; these unsuccessful trials are retained. The threshold is conservative and
is not presented as an optimum for every frame-history distribution.

[detection_records_review.csv](detection_records_review.csv) contains all **60
final-revision samples** (54 targeted samples and six threshold-control samples),
plus **12 serial-path controls** from the preceding batch-only checkpoint. Those
controls cover fork dense replay and spawn unique replay/dense detection. Their
executed serial path is unchanged by the final batch-local cache refinement;
all output digests matched. They are correctness controls with only two samples
per revision, not new speedup claims. Their checkpoint source hashes and labels
remain distinct. The original dense/spawn evidence above is also retained.
[detection_records_exploration.csv](detection_records_exploration.csv) preserves
all **96 exploratory samples**, including rejected raw-packing and threshold
probes and the non-memoizing batch implementation. Prototype rows are design
experiments, not measurements of the final revision; empty source-hash fields
come from the original runner before hash recording was added.

### Repair validation

The focused fork suite passes 130 tests covering detection records, bounded batch
compression, worker utilities and timed output. New tests cover byte budgets for
large/empty records, exact-history reuse without dropping records, serial fallback,
64-bit widening through the actual collection loop, delayed first-batch delivery,
bounded submission and disk spooling, progressive writing, preserved prior files,
and cleanup on compression/write failure. Before implementation, 13 new behavior
tests failed; the cache test separately failed on the non-memoizing implementation.
The earlier 35-test spawn run passed with its fork-only scheduling case skipped.
Applicable pre-commit hooks, the configured `ty==0.0.15` check and `git diff --check`
pass. The earlier full-suite and eight-pipeline checks are historical evidence
for the initial PR revision, not a newly repeated full campaign for this repair.

The unrelated CodSpeed cross-runner report is not used to establish or negate
these same-machine measurements. Re-review must assess this repair and its
retained evidence before merge.
