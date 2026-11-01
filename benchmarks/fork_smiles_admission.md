# Fork SMILES batch admission at higher worker counts

With 4,509 molecule records and 32 fork workers, the old `256 * nproc`
threshold selects per-record transport; the same input with four workers uses
batches. Per-record tasks repeatedly serialize the collector, including atom
metadata. This change caps the fork admission threshold at
`256 * min(nproc, 4)` records for `nproc > 1`, retaining the previously measured
two/four-worker thresholds while admitting medium inputs with larger pools.

The requested pool size and `2 * nproc` in-flight limit are unchanged. Existing
batch input/result budgets, singleton exceptions, spill behavior and ordered
consumption remain in use. Spawn/forkserver retain 8,192 records per worker;
single-process execution stays local. This is a measured admission heuristic,
not a claim that 1,024 records is optimal for every molecule or machine.

## Method

Baseline: official master `1d3a23a88d5c7adf34e58a40db9b7722ac3889c6`.
The candidate changes only `_should_batch_smiles` in production code. Its
`_smilesworker.py` SHA-256 is `ee5ec642fc88964281199bf0bce847a794ddef9dbabb41206be3ee0c3149c286`.
Both variants use the same native extension and Python 3.12.12 environment:
NumPy 2.5.4, RDKit 2026.9.1, NetworkX 3.7, Open Babel 3.2.1 and h5py 3.14.0.
The campaign uses frozen source overlays in one dependency environment; Git/file
hashes identify the code rather than the installed package's build-version label.
Measurements run sequentially on one Linux AMD EPYC 7K62 host, CPUs 0–31,
with native math threads limited to one and no concurrent test/build campaign.

The two real windows each contain 5,001 frames and 2,600 atoms from the same
2-chlorophenol/O2 trajectory at 2,500 K: early 0–50 ps and middle
1,250.02–1,300.02 ps. They use 32 workers, `miso=1`, no HMM and stride 1.
Fixed-record naming replays use 4,509 / 9,918 definitions and include conversion,
grouping and name-file writing. Full pipeline timers cover `rng.run()` from
DETECT through MATRIX. Imports, fixture preparation and output validation are
outside the timers. Every process explicitly pre-reads its input; fresh-process
pairs alternate AB/BA/AB. The middle pipeline was extended to five pairs
(AB/BA/AB/BA/AB) after its first three pairs were 1.6% slower for the candidate;
all original samples are retained. RSS is a process high-water mark including setup;
largest-child RSS is not the sum of simultaneous worker memory.

## Results

Medians in seconds (three pairs, five for the middle pipeline):

| Measurement     |   Base | Candidate | Change |
| --------------- | -----: | --------: | -----: |
| Early naming    |  2.598 |     0.703 | -72.9% |
| Middle naming   |  7.322 |     7.317 |  -0.1% |
| Early pipeline  | 42.182 |    38.451 |  -8.8% |
| Middle pipeline | 49.190 |    49.564 |  +0.8% |

The middle window already uses batches in both versions and is an unchanged
path control. Ranges and resource measurements (base / candidate):

- Early: pipeline ranges 41.535–42.246 s / 38.386–38.916 s; Step 3 medians 17.332 / 13.586 s; parent peak RSS 300.3 / 299.9 MiB; largest-child peak RSS 230.3 / 230.2 MiB.
- Middle: pipeline ranges 48.679–49.967 s / 48.795–50.608 s; Step 3 medians 22.222 / 22.179 s; parent peak RSS 328.5 / 329.8 MiB; largest-child peak RSS 243.1 / 249.2 MiB.

Repeated/mixed controls use the existing `smiles_batches.py` synthetic inputs;
they are stress workloads, not additional chemical observations. Each label
states its record and worker count. `large1024_32` also returns two large frame
arrays and includes incremental output hashing in both timers.

| Control       |  Base | Candidate | Change |
| ------------- | ----: | --------: | -----: |
| repeat1024_32 | 0.627 |     0.321 | -48.8% |
| repeat4096_32 | 2.555 |     0.553 | -78.4% |
| mixed1024_32  | 1.907 |     1.468 | -23.0% |
| mixed4096_32  | 9.122 |     7.056 | -22.6% |
| large1024_32  | 0.696 |     0.373 | -46.4% |
| repeat1024_8  | 0.398 |     0.178 | -55.3% |

Six unchanged-path controls each have one pair: 1,023 records with 32 fork
workers, 1,024/4 and 512/2 fork configurations, one process, and four-worker
spawn/forkserver runs. These check output compatibility, not speed. The 512/2
single pair changes from 0.138 to 0.151 s (+9.4%); its scheduling is unchanged
and one pair cannot establish a performance difference. All raw results,
including these controls, are in [naming samples](fork_smiles_admission_names.csv).

Separate [IPC audits](fork_smiles_admission_audit.csv) show early task count
4,509 -> 170, and application task serialization 201,391,060 -> 10,375,834 bytes.
Initializer metadata is transferred once per worker and excluded from these
per-task totals. The middle input retains exactly 270 tasks and 16,073,243 input
bytes. The large-result candidate spills 4,196,750 bytes across two SMILES result
streams, preserving exact output hashes. These counters exclude the shared
ordered spool's own files and runner framing; audit timings are not speed samples.

## Correctness and limits

675 non-GUI tests pass, with 6 GUI deselections, 6 xfailed and 3 non-strict XPASS.
The 15 added cases check admission boundaries, unchanged spawn/forkserver guards,
requested pool size, and exact names/timeline rows for all three `miso` modes.
Existing tests cover oversized inputs/results, ordering, filters, failed SMILES,
worker exit, serialization errors and cleanup.

All 60 timed naming/result samples and 10 separate audits preserve their exact
ordered output fingerprints. All 24 full-pipeline outputs pass native validation;
all 12 base/candidate pairs preserve exact species labels, atom/bond definitions,
presence ranges, populations, reaction counts and transition evidence after only
normalizing internal IDs/table order. Integration includes HMM-on `miso=0` and
no-HMM `miso=2` pilots plus HMM-on early/middle-window `miso=1` pairs. Those single
pairs are correctness checks. [Pipeline samples](fork_smiles_admission_pipeline.csv)
retain all timings and actual admission decisions.

These results cover bounded windows, one Linux host and up to 32 workers.
They do not establish a universal speedup or a total RSS bound. Spawn/forkserver
correctness controls ran on Linux; macOS/Windows performance was not measured.

## Reproduce

Install the base and candidate separately with identical dependencies; use this
checkout's benchmark scripts with either installed interpreter. Replace each
placeholder with the actual path. Choose an available fixed CPU set on Linux;
`0-31` matches this campaign. Use fresh output names for every repetition and
run samples sequentially, alternating version order.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0
RNG_BASE_PYTHON=/path/to/base/bin/python
RNG_HEAD_PYTHON=/path/to/head/bin/python
RNG_INPUT=/path/to/window.lammpstrj
RNG_STATE=/path/to/new-fixed-state

# Prepare one common stream, outside all timers.
taskset -c 0-31 "$RNG_BASE_PYTHON" benchmarks/isomer_candidates.py --prepare \
    --input "$RNG_INPUT" --atom-types C Cl H O --nproc 32 \
    --state "$RNG_STATE" --json prepared.json

# Each invocation explicitly warms its record input outside the timer.
taskset -c 0-31 "$RNG_BASE_PYTHON" benchmarks/fork_smiles_admission.py \
    --state "$RNG_STATE" --nproc 32 --method fork --json base-names.json
taskset -c 0-31 "$RNG_HEAD_PYTHON" benchmarks/fork_smiles_admission.py \
    --state "$RNG_STATE" --nproc 32 --method fork --json head-names.json

# Synthetic controls: repeat with both interpreters and all reported sizes.
taskset -c 0-31 "$RNG_HEAD_PYTHON" benchmarks/fork_smiles_admission.py \
    --case mixed --records 1024 --nproc 32 --method fork --json mixed.json

# Full pipeline; warm the trajectory separately before each process.
cat "$RNG_INPUT" > /dev/null
taskset -c 0-31 "$RNG_BASE_PYTHON" benchmarks/equal_z_candidates.py --mode pipeline \
    --input "$RNG_INPUT" --atom-types C Cl H O --nproc 32 --miso 1 \
    --output base-output --json base-pipeline.json
cat "$RNG_INPUT" > /dev/null
taskset -c 0-31 "$RNG_HEAD_PYTHON" benchmarks/equal_z_candidates.py --mode pipeline \
    --input "$RNG_INPUT" --atom-types C Cl H O --nproc 32 --miso 1 \
    --output head-output --json head-pipeline.json
"$RNG_HEAD_PYTHON" benchmarks/equal_z_compare.py base-output head-output
```

Use `--case repeated` or `--case large-output` for the other synthetic workloads.
Use `--records 1023` for the below-threshold control, `--nproc 1/2/4/8` for the
worker controls, and `--method spawn/forkserver` for those startup controls.
Pass one concrete value per option. Repeat selected cases with `--audit` only
after timing runs have finished; do not combine their timings with the table.
For the HMM integration controls, add `--hmm` to the pipeline command and select
`--miso 0` for the pilot or `--miso 1` for the real windows; the no-HMM pilot uses
`--miso 2`.
