# Bounded species population aggregation

This PR is based on upstream `b2d30fb4`, including molecule-presence reuse
(#2513) and repeated reaction-path aggregation (#2514). The population
implementation and its dedicated tests were unchanged when rebasing.

The performance measurements below were collected before the rebase, comparing
upstream `ecca5231abcbedfbd6e3a2806db30098ea3b8488` with candidate `9136cf32`.
That baseline already includes aromatic species normalization and #2513.
The comparison covers the complete `_GenerateMatrix._printspecies`
operation: reading stored observations, aggregating populations, and writing the
`.species` file. It does not measure the whole trajectory analysis pipeline.

## Behavior and resource limits

The output still counts **observed frame indices**, including repeated indices,
and preserves the order in which species first appear in each frame while
visiting molecules in file order. It does not substitute effective HMM presence.
The timestep text and output formatting are unchanged.

A lazy dictionary handles small workloads without allocating a Counter for every
empty frame. Its conservative estimate includes both populated frames and active
species cells, so a short trajectory with many species can also exceed its budget.
Larger workloads use compact, frame-major count/priority cells or block-partitioned
observation events, choosing the smaller final disk payload. A duplicate-heavy
single frame that would exceed the sparse sorting budget uses counts instead.

| Resource                     | Limit or selection rule                                          |
| ---------------------------- | ---------------------------------------------------------------- |
| Small dictionary             | 8 MiB estimated storage and at most 65,536 observations          |
| Count or sparse-sort block   | 16 MiB target; a dense block needs at least one species-wide row |
| Observation processing chunk | At most 65,536 event rows                                        |
| Observation spool            | At most 256 MiB logical payload; source is reread when exceeded  |
| Output formatting            | At most 1,200 entries per formatting batch                       |

Sorted unique molecule records contribute at most one count per frame; other
records contribute at most their length. The sum of these bounds determines a
safe count dtype, including when duplicate observations outnumber molecules.
The first contributing molecule determines each species' output priority.

The observation spool avoids another decode pass when it fits. Spool creation or
write failures discard only owned scratch files and reread the input. If the
final scratch layout cannot be created or populated, aggregation falls back to
bounded frame windows with sequential input rescans. Output errors propagate;
scratch files are still cleaned up.

These are aggregation workspace targets, **not a constant total RSS guarantee**.
One decoded molecule record, caller-owned metadata, the name table, compact
partition offsets, and OS-resident mapped pages remain outside the block estimate.
Partition metadata grows with the number of frame blocks. The final memory-mapped file can
grow with the data, and scratch space is not globally capped at the spool limit.
For exceptionally wide frames, even one dense count row can exceed 16 MiB.

## Reproducible comparison

Run [species_populations.py](species_populations.py) with separate baseline and
candidate installations using the same dependencies:

```sh
python benchmarks/species_populations.py \
    --baseline-python /path/to/baseline/.venv/bin/python \
    --candidate-python /path/to/candidate/.venv/bin/python \
    --repeats 5 \
    --output /tmp/species-populations.json
```

The driver alternates baseline/candidate order, starts a fresh process for every
sample, and fixes BLAS thread counts to one. It checks separate source paths,
matching Python/platform/NumPy versions, and SHA-256 equality of every output.
Input preparation and output hashing are outside the timed interval. CPU time is
process user plus system time; this operation does not create workers. Peak RSS
is the process high-water mark, including imports and fixture preparation, not
an allocation delta for this function.

The measurements below were collected on an Apple M4 with 24 GiB RAM,
macOS 15.6.1 arm64, Python 3.12.7 and
NumPy 2.5.3, using matched pinned runtime dependencies. No other validation job
from this task ran during the measurements. They are local samples, not a general
speedup guarantee. [species_populations.csv](species_populations.csv) retains
all 60 raw samples, output hashes, and candidate resource diagnostics.

All fixtures use seed `20261030` and timestep `100 + 10 * frame`:

| Case         |  Frames | Species / molecules | Observations                                                      |
| ------------ | ------: | ------------------: | ----------------------------------------------------------------- |
| small        |      32 |              4 / 12 | 8 random visits per molecule                                      |
| sparse       | 200,000 |            64 / 128 | 2 random visits per molecule                                      |
| sparse-large | 200,000 |           128 / 256 | 128 random visits per molecule                                    |
| dense        |  32,768 |            64 / 128 | Every molecule visits every frame                                 |
| mixed        |  32,768 |           128 / 256 | First 32 molecules visit every frame; others have 2 random visits |
| many-species |       8 |       4,096 / 8,192 | 2 random visits per molecule                                      |

Median measurements over five samples per revision:

| Case         | Wall seconds, base → candidate | CPU seconds, base → candidate | Peak RSS MiB, base → candidate | Candidate layout | Peak scratch payload MiB |
| ------------ | -----------------------------: | ----------------------------: | -----------------------------: | ---------------- | -----------------------: |
| small        |            0.000321 → 0.000337 |           0.000325 → 0.000341 |                  80.28 → 79.98 | small            |                    0.000 |
| sparse       |            0.151946 → 0.086657 |           0.151530 → 0.086660 |                  97.61 → 80.98 | small            |                    0.000 |
| sparse-large |            0.166218 → 0.146605 |           0.165592 → 0.141511 |                 101.81 → 87.67 | sparse           |                    0.375 |
| dense        |            0.844357 → 0.512452 |           0.837378 → 0.488596 |                 139.75 → 92.59 | dense            |                   20.000 |
| mixed        |            0.248023 → 0.265676 |           0.245702 → 0.231994 |                116.69 → 102.80 | sparse           |                    8.003 |
| many-species |            0.026805 → 0.028868 |           0.026797 → 0.028779 |                  80.88 → 81.02 | small            |                    0.000 |

The dense and long sparse cases improve both time and peak RSS in this run.
Mixed occupancy reduces RSS by about 12%, with wall time about 7% higher and CPU
time lower. Small and short many-species workloads have roughly unchanged RSS;
their median time increases are about 16 microseconds and 2.06 milliseconds.
These costs are retained in the report rather than described as universal gains.
All 60 samples produced byte-identical output for their corresponding case.

Scratch diagnostics count the maximum simultaneous logical observation and final
layout payload, including observations still in the bounded write buffer. They
exclude filesystem allocation overhead. `working_bytes` in the CSV describes the
count/sort block estimate or small-dictionary estimate; it is not measured RSS.

## Correctness validation

`tests/test_speciespopulations.py` compares the matrix-writing entry point with
an independent implementation of the original ordered Counters. It covers both
name-table representations; small, dense, and sparse paths; duplicate/unsorted
and negative frame indices; shared and Unicode names; empty timelines; sparse
partition refinement; short trajectories with many species; observation capacity,
creation and write failures; unavailable scratch storage; output failure cleanup;
and truncated observation records. Seeded randomized fixtures assert both layout
selection and exact output bytes.

```sh
python -m pytest -q tests/test_speciespopulations.py
python -m pytest -q tests -k 'not gui' --tb=short -ra
```

Before rebasing, the non-GUI suite passed **462 tests**, with 8 fork-only skips, 6 GUI
deselections, 6 expected failures and 3 non-strict expected-failure tests that
passed. Ruff 0.17.2, Velin 0.0.12, ty 0.0.15, Markdown formatting and whitespace
checks passed. Validation was local to macOS; Linux, Windows and GUI validation
were not performed locally.

The dedicated suite passed all **25 tests**. A separate integration tree combined
this change with route/reaction aggregation `5ab91b24`; its **117 targeted tests**
passed, including real worker processes and aromatic species normalization.

Eight complete pipeline configurations were compared with the upstream baseline,
both for this candidate alone and for the integration tree. Each comparison
matched **65 generated text/CSV/table files byte for byte**, including split
outputs, and **six public HDF5 semantic manifests**. The cases cover HMM on/off,
all three `miso` modes, one/two workers, origin retention, frame/timestep filters,
split/selectatoms, species-only execution, and ordinary text output. These checks
establish output equivalence for the fixtures, not whole-pipeline speedups.

## Validation after updating to upstream `b2d30fb4`

The branch was rebased onto the latest upstream on 2026-10-30, after #2514
merged. A range comparison confirmed that rebasing preserved the implementation
and tests without changes; this report records the new validation separately
from the original performance samples above.

The complete non-GUI suite passed **493 tests**, with 8 fork-only skips,
6 GUI deselections, 6 expected failures, and 3 non-strict expected-failure tests
that passed.

A targeted run using the `fork` multiprocessing start method passed **171 tests**
across utilities, indexed Step 3, route aggregation, reaction-pair aggregation,
species populations, and aromatic species normalization. This also exercises
tests that require fork and are skipped by the default macOS start method.

All eight complete pipeline configurations were rerun against `b2d30fb4`:
**65 generated text/CSV/table files matched byte for byte**, including split
outputs, and **six public HDF5 semantic manifests matched**. The six species
writing fixtures also produced identical output in a fresh-process correctness
pass against the updated baseline. That single pass ran alongside other checks
and is not a new performance measurement.

The comparison used separate source trees with matching Python 3.12.7 runtime
dependencies and reused the same compiled native extension, whose sources and
build configuration are unchanged between these revisions. This was a source
compatibility check, not a fresh wheel build. Ruff 0.17.2, Velin 0.0.12,
ty 0.0.15, Markdown formatting, and whitespace checks passed on the updated
branch. Validation remains local to macOS; Linux, Windows, and GUI validation
were not performed locally.
