# Bounded timeline range-read cache

Writing transition evidence checks that every participant exists at its reaction
frame. The original binary search reads one HDF5 scalar per comparison, even
when consecutive events revisit the same range-table pages. This change keeps
the same search and uses a writer-owned LRU of contiguous array blocks.

The three signed-int64 range columns share 16 blocks of at most 8,192 rows:
at most 1 MiB of cached numeric data, plus small array/dictionary headers.
Eviction happens before the next read. The lookup is created after the range
tables are flushed and discarded before file close, including on errors. It
does not cache a complete trajectory-sized index. The native validator keeps
using the original scalar lookup independently.

## Measurement protocol

Baseline: official `master` at `2d6c88b02242aa8ec221bcce5e5057d063b8979a`.
Both revisions use the same Python 3.12.12 environment and native extension:
Open Babel 3.2.1, h5py 3.14.0, NumPy 2.5.4, SciPy 1.18.1 and RDKit 2026.9.1.
The host is Linux on an AMD EPYC 7K62 with 48 cores; each pipeline process is
pinned to CPUs 0–31, requests 32 workers and limits native math libraries to
one thread. Timings run sequentially in fresh processes, alternating AB/BA/AB.

The two real windows contain 5,001 frames and 2,600 atoms each, sampled from
a 2-chlorophenol/O2 trajectory at 2,500 K: early 0–50 ps and middle
1,250.02–1,300.02 ps. Atom types are C, Cl, H, O; HMM is disabled, `miso=1`,
and stride is 1. Input bytes are warmed before each run. The timer covers
`ReacNetGenerator.run()` (DETECT through MATRIX), excluding imports, object
construction and output validation. A separate 101-frame pilot is excluded
from the medians. These are bounded windows on one Linux host, not a complete
production trajectory or a cross-platform speed claim.

Three-repeat medians (seconds):

| Window | Pipeline base | Pipeline cached | Change | Step 3 base | Step 3 cached |
| ------ | ------------: | --------------: | -----: | ----------: | ------------: |
| Early  |        48.690 |          42.081 | -13.6% |      24.246 |        17.724 |
| Middle |        57.870 |          51.535 | -10.9% |      30.291 |        23.581 |

Step 3 is 26.9% / 22.2% faster. The early pipeline ranges are 48.362–49.328 s
for base and 41.240–42.601 s for cached; middle ranges are 57.843–58.916 s
and 49.590–52.291 s. Parent peak RSS medians increase from 298.8 to 301.6 MiB
and from 328.5 to 329.5 MiB; largest-child peaks remain about 247 MiB.
The pilot measures 4.306 / 4.282 s and is not a speed claim.

Isolated presence-query medians, including lookup construction:

| Case     | Range rows | Queries | Scalar seconds | Cached seconds |
| -------- | ---------: | ------: | -------------: | -------------: |
| small    |          7 |       7 |       0.000390 |       0.000256 |
| repeated |    262,147 |   1,000 |       0.279463 |       0.021949 |
| random   |    262,147 |   1,000 |       0.288602 |       0.106199 |
| early    |  2,252,407 |   9,128 |       1.845186 |       1.019333 |
| middle   |  2,463,083 |   8,876 |       1.777224 |       0.800199 |

Separate audits replace 381,862 scalar range reads with 53,585 bounded slice
reads for early, and 374,097 with 40,921 for middle. Slices can copy more bytes
than scalar reads; this is a reduction in HDF5/Python calls, not a claim about
physical I/O bytes. The small control differs by about 0.13 ms. Both synthetic
large-table controls improve on this host, including low-locality random queries.
Standalone query timings exclude the worker-pool and writer interactions in
the full pipeline and should not be used to predict pipeline speedup.

Raw data: [pipeline samples](timeline_range_cache.csv) and
[query samples and separate audits](timeline_range_cache_queries.csv).

## Correctness and memory

The full non-GUI suite passes 643 tests, with 6 GUI deselections, 6 expected
failures and 3 non-strict XPASS results. New tests compare cached and scalar
queries against an independent membership oracle, covering empty tables,
missing IDs, gaps, inclusive endpoints and block boundaries. They also cover
shared capacity, LRU eviction, file isolation, repeated-query read counts,
read failures and rejection of an absent participant after warming the cache.
The latter preserves the previous complete destination on failure.

All 14 pipeline outputs pass the native HDF5 validator. All 7 base/candidate
pairs preserve frame mapping, exact names, atom membership, raw bond orders,
presence ranges, species populations, aggregate counts and transition evidence
after normalizing internal IDs and table ordering. All 30 query timing runs and
4 separate read audits match the original scalar query result for every input.

The array budget does not bound total process RSS. Pipeline measurements
record parent and largest-child peaks separately; they are not a concurrent
process-tree memory total. Query-replay RSS includes preparation and the
post-timing scalar validation. The cache's retained array bytes are reported
separately and remain within the stated budget.

## Reproduce

Use installed base/head environments with identical dependencies, or isolated
source snapshots with a common native extension. Run the same scripts in fresh
processes; alternate the order and repeat three times. Avoid concurrent builds
or tests on the measured CPUs.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Complete pipeline; run once for each environment and real window.
"$RNG_BASE_PYTHON" benchmarks/equal_z_candidates.py --mode pipeline \
    --input /path/to/window.lammpstrj --atom-types C Cl H O \
    --nproc 32 --miso 1 --output base-output --json base.json
"$RNG_HEAD_PYTHON" benchmarks/equal_z_candidates.py --mode pipeline \
    --input /path/to/window.lammpstrj --atom-types C Cl H O \
    --nproc 32 --miso 1 --output head-output --json head.json
"$RNG_HEAD_PYTHON" benchmarks/equal_z_compare.py base-output head-output

# Reuse one immutable timeline for isolated presence-query measurements.
"$RNG_HEAD_PYTHON" benchmarks/timeline_range_cache.py \
    --input base-output/timeline.h5 --json scalar.json
"$RNG_HEAD_PYTHON" benchmarks/timeline_range_cache.py --cached \
    --input base-output/timeline.h5 --json cached.json

# Synthetic controls need no external data. Repeat with --cached.
"$RNG_HEAD_PYTHON" benchmarks/timeline_range_cache.py \
    --case random --json random-scalar.json
```

The complete-pipeline commands use an explicit output directory; recorded
timings used the default input-prefixed output paths in isolated run directories.
Both variants use the same layout within each comparison. The `small`,
`repeated` and `random` query controls include frames outside their ranges;
missing-ID behavior is also covered by the independent unit oracle.
`random` exceeds the cache capacity. Add `--audit` to separate
runs to count HDF5 scalar/slice calls; audit timings are excluded from speed
claims, and read-call reductions do not imply physical disk-read reductions.
