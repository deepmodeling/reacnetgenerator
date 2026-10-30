# Single-process SMILES and HMM execution

This change removes worker-pool startup and record transport only when the
resolved `nproc` equals 1. SMILES conversion stays lazy and ordered, and HMM
filtering uses its already initialized model. Name caching, worker batching,
parallel scheduling and the shared `run_mp()` implementation are unchanged.

## Reproduce

Install the baseline and candidate into separate environments with identical
Python and dependency versions. Run this same script with each interpreter.
Without `--input`, it generates a small 8-atom, 32-frame H/H2 bond trajectory.
With an input, supply its file type and atom-type ordering explicitly.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

"$RNG_BASE_PYTHON" benchmarks/single_process.py --stage full
"$RNG_HEAD_PYTHON" benchmarks/single_process.py --stage full
"$RNG_HEAD_PYTHON" benchmarks/single_process.py --stage full --hmm

"$RNG_HEAD_PYTHON" benchmarks/single_process.py --stage names \
    --input /path/to/rp3.lammpstrj --input-type dump --atomname C H O
"$RNG_HEAD_PYTHON" benchmarks/single_process.py --stage hmm --hmm --origin \
    --input /path/to/rp3.lammpstrj --input-type dump --atomname C H O
```

Repeat **both** interpreters for every configuration, alternating AB/BA order
at least four times. Run performance samples sequentially, after tests and
builds finish. Compare `text_sha256`, `manifest_sha256`,
`record_signal_sha256`, input hashes and record counts before interpreting
any timing difference. `--miso 1` and `--miso 2` select merged names;
`--filtered` retains filtered-out HMM signals; `--nproc 2` exercises the
existing parallel branch. The default multiprocessing start method is spawn;
`--start-method` can select another platform-supported method.

Each invocation imports the **installed** package and reports its source paths,
source digests, dependency versions and platform. It emits one JSON sample.
Name and HMM stage timers exclude detection and earlier stages. The full timer
covers `ReacNetGenerator.run()` (DETECT through MATRIX), including HDF5 output;
network drawing and HTML report generation are excluded. Imports, input setup,
output hashing and official HDF5 validation are outside every timer.

CPU deltas include parent and terminated children, reported separately and
summed. Peak RSS includes imports and preparation but is sampled before
validation. Parent/child RSS values are separate high-water marks and cannot be
added to obtain concurrent process-tree physical memory. This removes transport
and process costs; it does not impose a new memory budget on individual graphs
or HMM signals.

The full comparison hashes every generated text artifact and the official
semantic manifest, retaining molecule IDs, order, ranges and event evidence.
The isolated HMM comparison hashes complete four-block molecule records paired
with their origin and effective HMM signals, allowing legal worker reordering.
No public validator or manifest rules are changed.

## Measurements (2026-10-30)

Baseline: `93f17b42b984e473a80c5da180468a981c1cb5e6`. Measured implementation:
`3784c7ad219004b5a62a6c727174412999c05151`; subsequent report-only changes do not
alter the measured code. All 35 tracked Python/Cython/C++/JSON package source
files in the baseline installation were checked against that official commit.

Apple M4, 24 GiB, macOS 15.6.1 arm64, Python 3.12.7; spawn, `nproc=1`,
`miso=0`, one BLAS thread. Dependencies: NumPy 2.5.3, RDKit 2026.3.6,
hmmlearn 0.3.3, LZ4 4.4.5 and h5py 3.14.0. Both environments use the same
dependency installation. Four fresh-process samples per variant/configuration,
strict AB/BA/AB/BA order, sequential execution. Desktop services remained
active; no CPU affinity or dedicated node was used.

The real input is `rp3.lammpstrj`: 5 frames, 12,326 atoms per frame, C/H/O,
5,615 detected molecule records. SHA-256:
`c7d7dce618f52b8e96451b4194c67ce63c976d22380e3b6d8c936bcb17503e85`.
Names and full runs enable the molecule timeline CSV; full runs also write
HDF5. The isolated HMM case enables origin output. The small case is the
script-generated H/H2 trajectory. These are one real trajectory and one
synthetic input, not a cross-system or chemical-accuracy benchmark.

All values below are medians; negative wall changes mean lower elapsed time.

| Case                           | Wall seconds, base → head | Change | Total CPU seconds, base → head | Parent peak MiB, base → head |
| ------------------------------ | ------------------------: | -----: | -----------------------------: | ---------------------------: |
| Real names, timeline enabled   |           3.8167 → 0.8751 | -77.1% |                4.5168 → 0.8707 |                262.5 → 270.5 |
| Real HMM, origin enabled       |           2.1458 → 0.4264 | -80.1% |                2.0451 → 0.4110 |                261.1 → 270.2 |
| Real DETECT → MATRIX, HMM off  |         17.6916 → 14.3520 | -18.9% |              18.3245 → 14.0793 |                271.3 → 270.3 |
| Real DETECT → MATRIX, HMM on   |         18.0629 → 13.5492 | -25.0% |              18.5615 → 13.4159 |                266.8 → 270.9 |
| Small DETECT → MATRIX, HMM off |           4.2567 → 2.4246 | -43.0% |                4.1987 → 2.3856 |                254.5 → 256.3 |
| Small DETECT → MATRIX, HMM on  |           6.6519 → 3.4103 | -48.7% |                6.4633 → 3.2255 |                251.6 → 254.1 |

| Case                           | Wall range, base (s) | Wall range, head (s) |
| ------------------------------ | -------------------: | -------------------: |
| Real names, timeline enabled   |        3.4413–6.5462 |        0.8125–1.5981 |
| Real HMM, origin enabled       |        2.0153–3.1161 |        0.3693–0.5280 |
| Real DETECT → MATRIX, HMM off  |      15.2985–28.6084 |      13.2993–17.8835 |
| Real DETECT → MATRIX, HMM on   |      16.7033–18.9156 |      13.3363–14.4964 |
| Small DETECT → MATRIX, HMM off |        4.1748–4.5255 |        2.3318–2.5358 |
| Small DETECT → MATRIX, HMM on  |        5.9482–9.3512 |        2.3885–4.6175 |

All 48 samples passed their configured equivalence checks. Every full run
passed the official validator and matched the official manifest and all text
hashes. HMM samples matched complete record/origin/effective-signal pairs.
The raw [CSV](single_process.csv) retains every sample, including timing
outliers, separate parent/child CPU and RSS, and output/source digests.

These results measure this patch alone, without name caching. They do not
reuse the earlier combined prototype speedups. Desktop timing variation is
visible in the ranges; these small sample counts do not establish universal
speedups. Linux, Windows and HPC performance were not measured.

## Validation

- Focused isomer, aromatic, molecule-presence and timed-output regression: 202 passed.
- Full non-GUI suite: 556 passed, 9 skipped, 6 deselected, 6 expected failures
  and 3 non-strict XPASS; 5 existing socket reducer tests were blocked by
  sandbox socket permissions. All 5 passed when rerun outside that sandbox
  (561 passing tests across the run and targeted retry).
- Full repository pre-commit hooks and changed-file hooks passed.
- `ty==0.0.15` check passed.
- Offline wheel build passed, including the native extension and webpack.
