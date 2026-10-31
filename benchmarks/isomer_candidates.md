# Narrowing VF2 candidates during isomer grouping

For `miso=1`, every new canonical name previously scanned all established
structural groups with the same formula. A formula can contain many different
labelled connectivities, so most of those VF2 calls can only reject a candidate.

The change indexes possible matches by a hash of the multiset of labelled
one-hop neighborhoods in `Molecule.graph`. That graph already applies the
requested isomer policy. Isomorphic graphs must have the same key; equal keys
do not prove isomorphism. Every possible match still uses the original VF2
node/edge predicates. Collisions only increase comparisons. Keys never become
Species identities or replace the canonical name, frequency or tie-break rule.

Indexing begins only after a formula has two distinct structural groups.
Unique formulas and repeated variants of one topology stay unindexed. `miso=0`
does not use this grouping entry; `miso=2` retains the composition-only scan.
Representatives keep their original order within each bucket. Retained keys are
integers, while temporary key construction scales with one decoded molecule.
Bucket dictionaries still scale with the number of retained structural groups;
this is not a fixed-memory cache or a total RSS guarantee.

## Measurements

The baseline is official master `5c278fbf50f4f8857229131ea19cd9c851619655`.
Both versions use Python 3.12.12, NumPy 2.5.4, NetworkX 3.7, RDKit 2026.9.1,
Open Babel 3.2.1 and h5py 3.14.0, with the same native extension. Samples run
sequentially on one Linux AMD EPYC 7K62 host, pinned to CPUs 0–31 with native
math threads limited to one. Three fresh-process pairs alternate AB/BA/AB.

Each real window contains 5,001 frames and 2,600 atoms from a 2-chlorophenol/O2
trajectory at 2,500 K. The early window covers 0–50 ps and the middle window
1,250.02–1,300.02 ps. The pipeline uses 32 workers, no HMM, `miso=1` and stride 1.
The trajectory is explicitly read before each timed process. Its timer covers
`rng.run()` (DETECT through MATRIX), excluding imports, construction and final
output validation. Separate 101-frame pilots cover the default settings,
HMM with `miso=0`, and no HMM with `miso=2`; pilots are not speed claims.

Pipeline medians (seconds):

| Window |   Base | Candidate | Change | Step 3 base | Step 3 candidate |
| ------ | -----: | --------: | -----: | ----------: | ---------------: |
| Early  | 41.977 |    41.816 |  -0.4% |      17.200 |           17.805 |
| Middle | 50.266 |    48.551 |  -3.4% |      23.253 |           21.777 |

The early total is effectively unchanged within the observed run spread:
base 41.104–42.361 s and candidate 41.208–42.408 s. Its Step 3 median is
0.605 s slower. Middle totals are 49.313–50.270 s and 48.488–48.554 s;
its Step 3 median decreases by 1.476 s. Parent peak RSS medians increase
from 299.6 to 299.9 MiB and from 328.6 to 329.5 MiB.

The naming replay uses a fixed ordered detection stream for each real window,
with 4,509 / 9,918 molecule records. These records were exported from the same
detection implementation before this change. It times the production
`_printmoleculename` entry, including conversion, grouping and name-file writing.
The helper explicitly reads the record file before every measurement. Exact
ordered name-table and definition-text hashes must match across all repeats.

Naming medians (seconds):

| Case       | Records |     Base | Candidate | Change |
| ---------- | ------: | -------: | --------: | -----: |
| small      |       3 | 0.001604 |  0.001770 | +10.3% |
| singletons |      64 | 0.022295 |  0.022476 |  +0.8% |
| repeated   |      64 | 0.024154 |  0.024239 |  +0.3% |
| distinct   |      64 | 0.093060 |  0.073439 | -21.1% |
| collisions |      32 | 0.018793 |  0.019830 |  +5.5% |
| mode2      |      64 | 0.020764 |  0.020442 |  -1.5% |
| mode0      |      64 | 0.005207 |  0.005246 |  +0.7% |
| early      |    4509 | 2.679538 |  2.882804 |  +7.6% |
| middle     |    9918 | 9.043198 |  7.345370 | -18.8% |

The middle replay improves by 18.8%, while early naming is 7.6% slower
(about 0.203 s). The small control adds 0.165 ms and the collision control
adds 1.037 ms; these costs are retained in the report. Repeated-topology
admission adds about 0.084 ms (0.35%). The optimization benefits sufficiently
diverse formula groups, not every input.

The synthetic controls use explicit candidate labels instead of RDKit SMILES
to isolate grouping. They cover three small records, unique formulas, repeated
topology with varying bond orders, different labelled trees sharing a formula,
and a prism/K3,3 pair whose neighborhood keys collide. Mode-0 and mode-2 controls
exercise unchanged equivalence policies. These are structural workloads, not
additional chemical observations. Separate `--audit` samples count grouping
VF2 calls and are excluded from timing medians.

Grouping-call audits (not timed speed samples):

| Case       | Base VF2 calls | Candidate VF2 calls |
| ---------- | -------------: | ------------------: |
| early      |            577 |                 290 |
| middle     |          18224 |                7326 |
| distinct   |           1909 |                 212 |
| repeated   |             63 |                  63 |
| collisions |             46 |                  46 |

The middle window avoids 59.8% of grouping VF2 calls. Full VF2 remains the
final decision for every surviving candidate, including deliberately colliding
keys. Raw data: [pipeline samples](isomer_candidates_pipeline.csv),
[naming samples](isomer_candidates_names.csv), and
[call audits](isomer_candidates_audit.csv).

These measurements do not establish a universal speedup. Low-diversity inputs
can pay extra indexing overhead; cold storage and other platforms were not
benchmarked. Parent RSS and largest-child RSS are separate process peaks, not
the sum of concurrent memory use.

## Correctness

The full non-GUI suite passes 660 tests, with 6 GUI deselections, 6 expected
failures and 3 non-strict XPASS results. New regressions cover shuffled atom IDs
and insertion order, isolated nodes and self-loops, all three isomer policies,
forced hash collisions, a non-isomorphic regular-graph collision, lazy indexing,
frequency ties and exact name/timeline output. Existing aromatic and failed
SMILES fallback tests retain their original expectations.

All 20 pipeline outputs pass native validation. All 10 base/candidate pairs
match exact species labels, atom/bond definitions, presence ranges, populations,
counts and transition evidence after normalizing internal IDs and table order.
The additional HMM + `miso=1` pilot exercises the optimized path with effective
HMM presence. All 54 naming samples and 10 audits retain identical ordered
name-table and definition-text fingerprints within each configuration.

## Reproduce

Install the base and candidate into separate environments with identical
dependencies. Run the scripts from this checkout using those installed
interpreters. Replace every placeholder below with the actual executable,
trajectory and fresh state-directory paths. On Linux, use an available fixed
CPU set and record it; `0-31` matches the host above. Use fresh output names for
every repeat, and run samples sequentially without concurrent tests or builds.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0
RNG_BASE_PYTHON=/path/to/base/bin/python
RNG_HEAD_PYTHON=/path/to/head/bin/python
RNG_INPUT=/path/to/window.lammpstrj
RNG_STATE=/path/to/new-window-state

# Prepare one common, fixed record stream; preparation is outside the timer.
taskset -c 0-31 "$RNG_BASE_PYTHON" benchmarks/isomer_candidates.py --prepare \
    --input "$RNG_INPUT" --atom-types C Cl H O --nproc 32 \
    --state "$RNG_STATE" --json prepared.json

# Each invocation explicitly warms the record file before its timer starts.
taskset -c 0-31 "$RNG_BASE_PYTHON" benchmarks/isomer_candidates.py \
    --state "$RNG_STATE" --nproc 32 --miso 1 --json base-names.json
taskset -c 0-31 "$RNG_HEAD_PYTHON" benchmarks/isomer_candidates.py \
    --state "$RNG_STATE" --nproc 32 --miso 1 --json head-names.json

# Synthetic controls use one process. Repeat with each environment.
taskset -c 0-31 "$RNG_HEAD_PYTHON" benchmarks/isomer_candidates.py \
    --case collisions --nproc 1 --miso 1 --json collisions.json

# Complete pipeline; explicitly warm the trajectory before each measurement.
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

For the other synthetic cases use `small`, `singletons`, `repeated` or `distinct`.
Repeat the `distinct` case with `--miso 2` and `repeated` with `--miso 0` for
policy controls. Add `--audit` in separate runs to obtain grouping call counts.
The pipeline reproduction uses explicit output directories; the recorded
pipeline samples used default input-prefixed paths inside isolated run folders.
Both variants used the same layout within each comparison.
