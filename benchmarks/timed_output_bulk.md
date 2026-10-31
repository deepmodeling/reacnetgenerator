# Bounded numeric appends in timed output

Large molecule definitions and fragmented presence histories currently make
one Python append call per numeric value. This patch fills the existing column
buffers from integer-array slices and bounded range chunks. Records shorter
than 32 values retain the scalar path. Schema 1.1, little-endian signed int64,
offsets, molecule/range order, maximal intervals, string columns and atomic
`.incomplete` publication remain unchanged.

The internal numeric API rejects non-integer arrays, non-vector shapes and
unsigned values beyond signed int64 before narrowing. It supports strided
and non-native-endian arrays. Conversion to Python scalars is limited to the
remaining row/byte capacity; it does not convert a whole large array to a list.
A failed later block may leave a prefix in the private incomplete file; it
never publishes that file over the previous destination.

Persistent columns retain the existing 8,192-row and 1 MiB encoded-byte caps.
Range scratch chunks use at most `min(_BLOCK_ROWS, _BLOCK_BYTES // 256)` rows
(4,096 with defaults), allowing space for Python tuples/integers and numeric
views within the existing byte allowance. Scratch references are released
before the next chunk; decoded atom/bond arrays are released before range
writing. Already-decoded individual records, Python object overhead and HDF5
caches remain outside encoded buffer accounting. This is not a total-RSS cap.

## Reproduce

Install official baseline `454d30816fda486d702a82df4ec3b59fa01d8d38` and the
candidate into separate environments with the same dependencies on each host.
Use the same benchmark script from this PR with both installed interpreters.
It reuses the existing `single_process.py` metric helpers and
`smiles_batches.py` fixed-input preparation helpers.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

# One baseline preparation for reproducible real-input stage comparisons.
"$RNG_BASE_PYTHON" benchmarks/smiles_batches.py --prepare \
    --input /path/to/rp3.lammpstrj --state /path/to/fixed-records

"$RNG_BASE_PYTHON" benchmarks/timed_output_bulk.py \
    --case real --state /path/to/fixed-records --stage writer
"$RNG_HEAD_PYTHON" benchmarks/timed_output_bulk.py \
    --case real --state /path/to/fixed-records --stage writer

"$RNG_HEAD_PYTHON" benchmarks/timed_output_bulk.py \
    --case real --state /path/to/fixed-records --stage full --hmm
"$RNG_HEAD_PYTHON" benchmarks/timed_output_bulk.py \
    --case fragmented --atoms 256 --frames 8192 --stage writer
"$RNG_HEAD_PYTHON" benchmarks/timed_output_bulk.py \
    --case fragmented --atoms 256 --frames 8192 --stage path
"$RNG_HEAD_PYTHON" benchmarks/timed_output_bulk.py \
    --case array --atoms 32769 --stage writer
```

Repeat **both** interpreters for every configuration in AB/BA/AB/BA order.
Use `--method fork` for Linux and `--method spawn` for macOS. All cases resolve
`nproc=1`; this PR changes no multiprocessing policy. `--case small` generates
a changing 8-atom, 32-frame H/H2 trajectory. `--input` can replace `--state` for
real input, at the cost of repeating DETECT preparation outside stage timers.
Full runs always start from the trajectory and include DETECT through MATRIX.

The writer timer covers metadata, frame/atom tables, molecule definitions,
presence ranges, final flush and atomic publication. Naming/filtering are
prepared outside this timer, and the writer-only file contains no reaction
events. PATH includes naming, routes, reactions, evidence and HDF5 writing;
DETECT/HMM preparation and the subsequent MATRIX stage are outside PATH.
HMM-on stage preparation is deterministic and single-process. Every run uses
the official semantic manifest, which performs full validation, and compares
all generated text hashes. No manifest rules or output readers are changed.

`--audit` instruments persistent column buffers and observes row/encoded-byte
peaks and flush counts. Run it separately: per-value instrumentation affects
runtime. This accounting excludes caller input arrays, bounded scratch objects
and HDF5 caches; actual process RSS is recorded independently. Parent and
terminated-child CPU deltas are separate. RSS is a high-water mark including
imports/preparation, sampled before validation. Child RSS is a largest-child
high-water mark, not concurrent process-tree memory; the RSS columns cannot
be added to obtain physical memory use.

## Fixed versions, inputs and conditions

Measured production implementation: `52719e76` (the array-lifetime correction
on top of `3872759f`). Subsequent validation/report-only commits do not change
the measured production file. All 36 Python/Cython/C++/JSON package source
files were checked against the official baseline and across the two hosts.
Exploratory measurements before the lifetime correction are excluded below.
The production `_timedoutput.py` SHA-256 hashes are:

- Baseline: `67da461c7caf24b32e3e0e93243a71a885006e93ef3dc53ce644e724741898bc`
- Candidate: `444df5ee5534bb706c7aef37773b664addbae2774267ee1db8f3ff28ee53afc1`

macOS: Apple M4, 24 GiB, macOS 15.6.1 arm64, Python 3.12.7;
NumPy 2.5.3, h5py 3.14.0, RDKit 2026.3.6, hmmlearn 0.3.3, LZ4 4.4.5.
Linux: Intel Core i5-14400F, 16 available logical CPUs, about 31 GiB,
Linux 7.0.0-38-generic x86_64, Python 3.12.9;
NumPy 2.5.4, h5py 3.14.0, RDKit 2026.9.1, hmmlearn 0.3.3, LZ4 4.4.4.
The variants share dependencies within each host. Desktop services and other
Linux users remained active; no CPU affinity or exclusive-node guarantee was
used. Do not treat absolute cross-host values as an OS/start-method comparison.

The real input has 5 frames and 12,326 C/H/O atoms, with 5,615 baseline detection
records. SHA-256:
`c7d7dce618f52b8e96451b4194c67ce63c976d22380e3b6d8c936bcb17503e85`.

Fragmented stress replays contain pairs of identical monomer definitions with
alternating presence, so each atom is represented once per frame without
chemical reaction events. They generate 512, 65,536 or 2,097,152 ranges.
These are synthetic ordered table inputs, not parser trajectories or additional
physical observations. The array case writes one 32,769-atom chain with 32,768
bonds and two occupied frames; its name is a placeholder, and no chemical
conversion performance is inferred. These cases isolate numeric-write costs.

## Measured results

Four fresh-process samples per variant and case, alternating AB/BA/AB/BA.
Wall times are median [minimum, maximum] in seconds; change is the ratio
of medians minus one (negative means faster). These are descriptive results,
not confidence intervals. See [all timing/resource samples](timed_output_bulk.csv)
and [separate buffer audits](timed_output_bulk_buffers.csv).

### macOS

| Case              |            Baseline wall s |           Candidate wall s | Change |
| ----------------- | -------------------------: | -------------------------: | -----: |
| small-writer      |    0.0075 [0.0048, 0.0109] |    0.0058 [0.0054, 0.0062] | -22.6% |
| small-full        |    2.2995 [2.2037, 2.9180] |    2.2972 [2.2391, 3.5491] |  -0.1% |
| small-full-hmm    |    2.2819 [2.2416, 2.3580] |    2.2130 [2.0721, 2.3433] |  -3.0% |
| real-writer       |    0.1802 [0.1608, 0.2508] |    0.1752 [0.1606, 0.3071] |  -2.8% |
| real-writer-hmm   |    0.1709 [0.1482, 0.1858] |    0.1412 [0.1355, 0.2510] | -17.4% |
| real-path         |    2.0895 [1.7887, 2.3478] |    1.9980 [1.9728, 2.3421] |  -4.4% |
| real-path-hmm     |    1.0405 [0.9087, 1.2006] |    1.0421 [0.8917, 1.1261] |  +0.2% |
| real-full         | 14.8723 [14.5209, 17.2552] | 14.5304 [14.1126, 14.8100] |  -2.3% |
| real-full-hmm     | 13.9101 [13.3683, 14.5144] | 14.0218 [13.5721, 14.4207] |  +0.8% |
| fragmented-tiny   |    0.0076 [0.0057, 0.0100] |    0.0073 [0.0063, 0.0089] |  -4.4% |
| fragmented-medium |    0.0870 [0.0845, 0.1158] |    0.0966 [0.0880, 0.1185] | +11.0% |
| fragmented-large  |    2.6656 [2.5415, 2.8445] |    1.1941 [1.1634, 1.2891] | -55.2% |
| fragmented-path   | 21.8987 [19.8714, 23.3935] | 19.3580 [18.7600, 19.7496] | -11.6% |
| array             |    0.0778 [0.0716, 0.0811] |    0.0138 [0.0120, 0.0161] | -82.3% |

### Linux

| Case              |            Baseline wall s |           Candidate wall s | Change |
| ----------------- | -------------------------: | -------------------------: | -----: |
| small-writer      |    0.0127 [0.0118, 0.0149] |    0.0125 [0.0121, 0.0129] |  -1.9% |
| small-full        |    0.3123 [0.2766, 0.3605] |    0.3073 [0.2676, 0.3427] |  -1.6% |
| small-full-hmm    |    0.2796 [0.2478, 0.3775] |    0.3459 [0.2658, 0.3996] | +23.7% |
| real-writer       |    0.1920 [0.1893, 0.2156] |    0.1916 [0.1898, 0.1922] |  -0.2% |
| real-writer-hmm   |    0.1407 [0.1397, 0.1435] |    0.1409 [0.1401, 0.1424] |  +0.1% |
| real-path         |    2.0944 [2.0785, 2.0956] |    2.0923 [2.0593, 2.1778] |  -0.1% |
| real-path-hmm     |    0.8203 [0.7625, 0.8308] |    0.8088 [0.7839, 0.8722] |  -1.4% |
| real-full         | 26.3186 [25.9785, 26.7589] | 26.2461 [26.2030, 26.5046] |  -0.3% |
| real-full-hmm     | 24.8346 [24.7120, 25.5078] | 25.0392 [24.8764, 25.4880] |  +0.8% |
| fragmented-tiny   |    0.0121 [0.0119, 0.0204] |    0.0120 [0.0119, 0.0121] |  -1.5% |
| fragmented-medium |    0.0840 [0.0833, 0.1021] |    0.0450 [0.0435, 0.0483] | -46.4% |
| fragmented-large  |    2.3356 [2.2986, 2.4040] |    0.9965 [0.9886, 1.1201] | -57.3% |
| fragmented-path   | 20.1233 [19.7327, 21.5767] | 18.4317 [18.2773, 18.9115] |  -8.4% |
| array             |    0.0647 [0.0638, 0.0703] |    0.0112 [0.0111, 0.0164] | -82.7% |

### Linux small-input follow-up

The initial HMM-on small full-pipeline control fluctuated enough to warrant
a separate eight-pair AB/BA follow-up, with the HMM-off full run as a matched
control. These runs followed the main timing campaign and buffer audits.
The initial results remain above and both campaigns are retained in the CSV.

| Case             |         Baseline wall s |        Candidate wall s | Change | Parent CPU s, base → head |
| ---------------- | ----------------------: | ----------------------: | -----: | ------------------------: |
| small-full       | 0.2552 [0.2411, 0.2787] | 0.2574 [0.2365, 0.2673] |  +0.8% |           0.1963 → 0.1898 |
| small-full-hmm   | 0.2387 [0.2274, 0.2669] | 0.2533 [0.2095, 0.2763] |  +6.1% |           0.1910 → 0.1919 |
| small-writer-hmm | 0.0069 [0.0063, 0.0071] | 0.0067 [0.0064, 0.0071] |  -3.1% |           0.0069 → 0.0067 |
| small-path-hmm   | 0.0578 [0.0465, 0.0678] | 0.0659 [0.0528, 0.0711] | +13.9% |           0.0172 → 0.0180 |

The writer/PATH diagnostic was run afterward with the same small HMM-on
input, eight pairs each, using `--case small --hmm` with `--stage writer`
or `--stage path`. These measurements isolate the modified writer from
DETECT/HMM preparation and MATRIX, without changing benchmark source.

### macOS medium-input follow-up

The initial 65,536-range writer sample was unfavorable. An eight-pair
AB/BA follow-up, after the main campaign and audits, is retained separately.

|         Baseline wall s |        Candidate wall s | Change | Parent CPU s, base → head |
| ----------------------: | ----------------------: | -----: | ------------------------: |
| 0.1012 [0.0905, 0.1128] | 0.1081 [0.0966, 0.1200] |  +6.8% |           0.0977 → 0.1030 |

### CPU and resident memory

Selected medians, baseline → candidate. Parent CPU excludes child CPU.
RSS includes setup/imports; raw samples retain child RSS, every case, ranges
and load averages. A load average is not a measure of exclusive CPU access.

| Host / case              | Parent CPU s |   Child CPU s | Parent peak MiB |  Child peak MiB |
| ------------------------ | -----------: | ------------: | --------------: | --------------: |
| macOS / small-full       |  0.23 → 0.21 |   2.02 → 2.07 | 256.03 → 255.35 | 252.32 → 251.59 |
| macOS / real-writer      |  0.18 → 0.17 |   0.00 → 0.00 | 257.16 → 256.06 |     1.00 → 1.00 |
| macOS / real-full        |  1.79 → 1.74 | 13.03 → 12.68 | 269.73 → 269.04 | 287.94 → 279.83 |
| macOS / real-full-hmm    |  1.20 → 1.15 | 12.52 → 12.81 | 266.95 → 267.70 | 282.91 → 285.34 |
| macOS / fragmented-large |  2.61 → 1.19 |   0.00 → 0.00 | 257.71 → 258.00 |     1.00 → 1.00 |
| macOS / fragmented-path  |  3.51 → 2.25 | 18.03 → 18.02 | 271.50 → 272.93 |   93.35 → 93.19 |
| macOS / array            |  0.08 → 0.01 |   0.00 → 0.00 | 255.95 → 255.48 |     1.00 → 1.00 |
| Linux / small-full       |  0.19 → 0.19 |   0.02 → 0.02 | 237.74 → 237.56 | 150.85 → 150.91 |
| Linux / real-writer      |  0.19 → 0.18 |   0.00 → 0.00 | 238.93 → 239.10 |   42.78 → 42.70 |
| Linux / real-full        |  2.33 → 2.32 | 23.96 → 23.91 | 245.55 → 245.43 | 175.79 → 175.89 |
| Linux / real-full-hmm    |  1.21 → 1.21 | 23.57 → 23.78 | 245.38 → 245.61 | 175.64 → 176.48 |
| Linux / fragmented-large |  2.24 → 0.93 |   0.00 → 0.00 | 238.75 → 239.34 |   42.69 → 42.80 |
| Linux / fragmented-path  |  3.76 → 2.48 | 17.61 → 17.30 | 247.15 → 247.80 | 169.38 → 169.96 |
| Linux / array            |  0.06 → 0.01 |   0.00 → 0.00 | 236.05 → 235.91 |   42.63 → 42.81 |

### Observed persistent buffers

Separate instrumented runs; encoded bytes exclude Python objects and scratch.
The observed row caps and unchanged flush counts preserve HDF5 block sizes.

| Host / case / variant           | Largest column rows | Largest column bytes | Total buffer peak bytes | Flushes | Largest numeric input rows |
| ------------------------------- | ------------------: | -------------------: | ----------------------: | ------: | -------------------------: |
| Linux / array / base            |                8192 |                65536 |                  196761 |      38 |                          0 |
| Linux / array / head            |                8192 |                65536 |                   65705 |      38 |                       8192 |
| Linux / fragmented-large / base |                8192 |                65536 |                  215119 |     783 |                          0 |
| Linux / fragmented-large / head |                8192 |                65536 |                  198732 |     783 |                       4096 |
| Linux / real-writer / base      |                8192 |                65536 |                  563047 |      30 |                          0 |
| Linux / real-writer / head      |                8192 |                65536 |                  563047 |      30 |                       8192 |
| macOS / array / base            |                8192 |                65536 |                  196805 |      38 |                          0 |
| macOS / array / head            |                8192 |                65536 |                   65749 |      38 |                       8192 |
| macOS / fragmented-large / base |                8192 |                65536 |                  215163 |     783 |                          0 |
| macOS / fragmented-large / head |                8192 |                65536 |                  198776 |     783 |                       4096 |
| macOS / real-writer / base      |                8192 |                65536 |                  563007 |      30 |                          0 |
| macOS / real-writer / head      |                8192 |                65536 |                  563007 |      30 |                       8192 |

### Interpretation and limits

The 2,097,152-range writer median falls 55.2% on macOS and 57.3% on Linux;
complete PATH falls 11.6% and 8.4%, respectively. The large chain writer falls
82.3% and 82.7%. These stress inputs establish the benefit for large numeric
writes; they do not establish additional chemical accuracy or a general
whole-pipeline speedup. The real trajectory's full-pipeline changes range from
-2.3% to +0.8%, with overlapping sample ranges.

Unfavorable results remain relevant. The macOS 65,536-range writer was +11.0%
in the initial campaign and +6.8% (6.9 ms) in its follow-up. The Linux small
HMM full run was initially +23.7%; the follow-up was +6.1% (14.6 ms), while
isolated HMM PATH was +13.9% (8.0 ms). The isolated HMM writer was -3.1%
(0.2 ms). These observations do not justify a universal no-regression claim,
nor do they identify the cause of the integration differences. They are not
removed or silently replaced with the better controls. No additional caching,
range-reader rewrite or threshold expansion was introduced to optimize them.

Audited persistent encoded-buffer peaks did not increase. Large fragmented
writer parent RSS medians rose by about 0.3 MiB on macOS and 0.6 MiB on Linux;
complete PATH rose by about 1.4 MiB and 0.7 MiB. These include preparation and
allocation variability, so they are not direct scratch-allocation measurements.

## Validation

The final macOS non-GUI run passed **613 tests**, with 9 skips, 6 deselections,
6 expected failures and 3 non-strict XPASS results. The Linux focused
`test_timedoutput.py`, `test_timedoutput_validation.py` and new bulk suite passed
**111 tests**. The new file contributes 28 cases, covering buffer boundaries,
large and strided inputs, endian conversion, invalid dtype/shape/value rejection,
signed overflow, HMM effective ranges, offsets/alignment, bounded generator
consumption and decoded-array release. Injected HDF5/validation failures leave
the old destination intact and a closed, validator-rejected incomplete sibling.

A final macOS wheel built from `52719e76`; Linux baseline and candidate native
editable installations built successfully. `ty==0.0.15` passed on the final
production implementation.

Full repository and changed-file pre-commit hooks passed, including formatting,
license headers, file-size checks and Python/Cython linting. `git diff --check`
passed. GUI execution and other Python/platform combinations were not rerun.

All **224 main timing samples**, **80 follow-up samples** and **12 separate
buffer audits** passed official validation and matched semantic manifests plus
text hashes between variants within each host/configuration. Dependency
versions match within each host; cross-host dependency versions differ.
CSV `campaign` separates initial timing, diagnostic follow-ups and audits.
No tests, builds or other benchmark campaigns ran concurrently on the same
host during these final measurements.

The earlier 52% writer-only prototype improvement motivated this PR; it is not
used as evidence for the final patch or as a whole-pipeline acceleration claim.
The measured integration results and small-input controls are retained even
where their changes are negligible or unfavorable.
