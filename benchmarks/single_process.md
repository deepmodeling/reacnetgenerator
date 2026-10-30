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
