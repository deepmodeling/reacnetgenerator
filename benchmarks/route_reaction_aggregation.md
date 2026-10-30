# Route and reaction aggregation

This change reduces repeated Step 3 results before retaining them in the parent
process. It is based directly on upstream `ecca5231abcbedfbd6e3a2806db30098ea3b8488`,
including aromatic-species normalization (#2512) and molecule-presence reuse
(#2513). The validation and measurements below were rerun after this update;
they replace the earlier comparison against `97efe9ff`.

- Workers count repeated molecule pairs in bounded blocks. Short or apparently
  unique routes keep the original pair array; compaction abandons its dictionary
  after 65,536 unique pairs. Narrow unsigned pairs use exact 64-bit integer keys
  for counting; wider or signed inputs retain both original fields.
- The parent accumulates species-pair counts directly. Without HMM, every
  occurrence on the first contributing atom counts; later atoms do not contribute
  that molecule pair. With HMM, each molecule pair counts once globally. Species
  names are resolved after these rules, and equal-species transitions are removed.
- Route output uses a 1 MiB byte budget in addition to the existing item limit.
  Ordering, formatting, split windows, atom selection, and nonzero-position route
  coordinates remain unchanged. One indivisible atom route can exceed the budget.
- Reaction graphs retain each molecule-pair edge once per transition. Every
  changed atom still contributes conflict evidence. Connected reaction instances
  and their participant IDs remain distinct even when their species names match.

These changes do not bound all memory. The route ownership map grows with unique
molecule pairs, a reaction graph grows with distinct edges, and individual atom
routes and their text are still materialized. Worker scheduling, count-only
transition batching, species population output, and molecule detection are
unchanged. Legacy array-based matrix and reaction entry points remain supported.

## Reproduce the comparison

Install the baseline and this revision into separate environments with matching
Python and dependency versions, then run from this checkout:

```sh
python benchmarks/route_reaction_aggregation.py \
    --baseline-python /path/to/baseline/.venv/bin/python \
    --candidate-python /path/to/candidate/.venv/bin/python \
    --repeats 5 --output /tmp/route-reaction-aggregation.json
```

The script alternates baseline/candidate order and starts a fresh interpreter for
every sample. It checks source locations and records environment versions. It
rejects unequal results or mismatched measured dependency versions. BLAS thread
counts are fixed to one. The script supports macOS and Linux (`resource` is used
for CPU time and peak RSS).

Timings include real indexed worker startup, result transfer, parent reduction,
and output writing. Route timings also include the matrix consumer's conversion
to species counts. Input construction, module imports in the parent, and output
hashing are excluded. CPU time sums parent and reaped children. Parent and child
peak RSS are process high-water marks including setup; child RSS is the maximum
reported child peak, not simultaneous total pool memory. Do not add these peaks.

All samples require exact agreement for route text hashes, species counts, and
raw active-transition flags, or reaction count text, event CSV, and ordered
participant evidence. Route mode 0/1 means no-HMM/HMM; reaction mode 0/1 means
counts-only/events with participant evidence.

## Correctness checks

- 468 installed-package non-GUI tests passed (8 fork-only skips, 6 GUI
  deselections, 6 pre-existing expected failures and 3 pre-existing unexpected
  passes). The suite includes the upstream aromatic-species and molecule-presence
  tests as well as all aggregation/indexed-worker/benchmark regressions.
- 114 utility/indexed-worker/aggregation tests passed under `fork`, including the
  8 tests skipped by the default `spawn` run.
- Eight complete pipeline configurations matched the latest baseline's text
  hashes and validated HDF5 semantic manifests, covering HMM/no-HMM, isomer modes,
  splits, selected atoms, species-only output, ordinary text, and timed evidence.
  Because #2513 is now in the baseline, these also test the combined behavior
  with molecule-presence reuse.
- Both revisions were independently built and installed with identical dependency
  versions. The candidate's installed Python modules matched the checkout exactly.
  Ruff lint/format, the pinned `ty==0.0.15` type checker, all applicable pre-commit
  hooks, and `git diff --check` passed.
- Regression cases cover first-atom multiplicity, repeated species names, bounded
  lookups, uint8/uint16/uint32 packing boundaries, uint64 fallback, compaction
  fallback, late conflict flags, cross-block connected components, and real
  `spawn`/`fork` workers.

The checks above ran on macOS. The local `fork` result is not Linux validation;
no dedicated Linux baseline/candidate pipeline comparison is claimed. The PR's
Linux CI result is reported separately on GitHub.

## Measurements

Remeasured on 2026-10-30 against `ecca5231` using an Apple M4 with 24 GiB RAM,
macOS 15.6.1 arm64, Python 3.12.7 and `spawn`. Both environments used NumPy 2.5.3,
pandas 3.0.6, h5py 3.14.0, lz4 4.4.5, SciPy 1.18.1 and RDKit 2026.3.6.
Each cell below is the median of five fresh-process samples. All 180 samples
passed their semantic comparison. [Raw samples](route_reaction_aggregation.csv)
include CPU time and semantic hashes. Other verification commands in this task
finished before measurement began.

Route inputs are 8 atoms × 32 frames (small), 32 × 32,768 (repeat, one or two
workers), 4 × 32,768 (unique), 64 × 32,768 (sparse, changing every 8,192 frames),
and 16 × 16,384 (half repeated, half unique; two workers). Other cases use one
worker. Reaction inputs have two frames and either 524,288 repeated atom edges,
8,192 independent unique edges, or 8,192 unique edges sharing one product.

These are synthetic stage benchmarks, not an end-to-end trajectory speedup claim.

| Case                 | Mode | Wall time (s), base → candidate | Parent peak RSS (MiB), base → candidate | Child peak RSS (MiB), base → candidate |
| -------------------- | ---- | ------------------------------- | --------------------------------------- | -------------------------------------- |
| route-small          | 0    | 0.0007 → 0.0006                 | 125.8 → 126.2                           | 1.0 → 1.0                              |
| route-small          | 1    | 0.0007 → 0.0006                 | 125.7 → 126.3                           | 1.0 → 1.0                              |
| route-repeat         | 0    | 0.6910 → 0.4767                 | 164.5 → 138.3                           | 98.3 → 96.9                            |
| route-repeat         | 1    | 0.8526 → 0.4923                 | 165.6 → 137.8                           | 98.7 → 97.8                            |
| route-repeat-pool    | 0    | 0.7277 → 0.4450                 | 164.8 → 138.9                           | 97.4 → 96.6                            |
| route-repeat-pool    | 1    | 0.9034 → 0.4725                 | 166.0 → 137.8                           | 97.4 → 97.1                            |
| route-unique         | 0    | 0.4563 → 0.4022                 | 175.8 → 149.0                           | 96.6 → 97.5                            |
| route-unique         | 1    | 0.4326 → 0.4851                 | 136.8 → 147.8                           | 97.4 → 97.5                            |
| route-sparse         | 0    | 0.3501 → 0.3440                 | 129.4 → 130.1                           | 90.7 → 90.7                            |
| route-sparse         | 1    | 0.3613 → 0.3720                 | 129.7 → 130.3                           | 90.5 → 91.0                            |
| route-mixed          | 0    | 0.9113 → 0.7942                 | 181.7 → 149.0                           | 94.1 → 93.4                            |
| route-mixed          | 1    | 0.8612 → 0.8191                 | 138.1 → 149.1                           | 93.9 → 93.7                            |
| reaction-repeat      | 0    | 0.5075 → 0.1660                 | 129.3 → 128.9                           | 107.9 → 53.2                           |
| reaction-repeat      | 1    | 0.4745 → 0.1698                 | 129.2 → 129.1                           | 112.1 → 53.9                           |
| reaction-unique      | 0    | 0.2301 → 0.2202                 | 126.6 → 127.1                           | 54.4 → 55.2                            |
| reaction-unique      | 1    | 0.2494 → 0.2503                 | 127.6 → 128.2                           | 58.6 → 59.2                            |
| reaction-high-degree | 0    | 0.1670 → 0.1705                 | 126.1 → 126.4                           | 49.2 → 50.8                            |
| reaction-high-degree | 1    | 0.1668 → 0.1705                 | 126.9 → 127.0                           | 49.1 → 50.8                            |

Repeated routes reduced wall time by 31–48% and parent peak RSS by
25.9–28.2 MiB. Repeated reaction edges reduced wall time by 64–67%;
the table also records the reduction in child peak RSS.

The unique/HMM route case changed median wall time by +12.1%
and parent peak RSS by +11.1 MiB. Mixed/HMM changed parent peak RSS
by +11.0 MiB. Global ownership state now lives
in a Python map rather than the baseline's event arrays followed by a NumPy
unique pass. Sparse routes and unique/high-degree reaction graphs do not establish
a general speedup; all medians, including regressions, are shown above. The small
serial case is below 1 ms. These tradeoffs and the five-sample noise limit any
general performance claim.
