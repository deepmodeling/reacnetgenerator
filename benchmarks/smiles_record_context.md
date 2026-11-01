# Per-record SMILES conversion state

Small SMILES inputs use ordered per-record tasks. Previously each task serialized
the bound collector, including its RNG owner and frame-sized timestep dictionary.
The standard collector now initializes the existing minimal converter once per
worker and submits a module-level function with the four compressed record blocks.
Batch admission, worker count, the two-records-per-worker admission bound and
disk-ordered result consumption retain their existing settings. Custom collectors
and instance callables keep the bound-method path; single-process conversion stays
local. No name cache is introduced.

The initializer carries element names, per-atom names, a need-frames flag and a
frequency/name selector. Timeline filtering, timestep lookup, VF2 fallback and
global-frequency representative selection remain in the parent. This removes
trajectory metadata from tasks; it does not bound parent RNG storage, per-atom
metadata or an individual molecule's decoded frames.

## Matched measurements

Baseline is official `dd47d2fde3d5db6b8e671752dc3a45af2f6917e9`. Measured production
SHA-256 values:

| Module             | Baseline                                                           | Candidate                                                          |
| ------------------ | ------------------------------------------------------------------ | ------------------------------------------------------------------ |
| `_path.py`         | `465e6beed007e927849ff6b98f50b97dabc28f4e462cc40eabedf99fdbaf86e1` | `4bcbc923f968e0da54a0a013e3d0e7dfc1e6f36358c44bd900df87e8bc0b7360` |
| `_smilesworker.py` | `ee5ec642fc88964281199bf0bce847a794ddef9dbabb41206be3ee0c3149c286` | `423637178c177f839c5c2089a3f3fc8413773e8a032a92db81a637b629bfdb18` |

All runs are sequential fresh processes on one Linux AMD EPYC 7K62 host, CPUs
0–7, with native math threads limited to one. Both versions use the same Python
3.12.12 environment, NumPy 2.5.4, RDKit 2026.9.1, NetworkX 3.7, LZ4 4.4.5 and
native extension. Source hashes identify the measured code independently of
generated package version labels. Each row uses three alternating AB/BA/AB pairs.

Synthetic fixtures construct actual RNG/collector objects, retaining a timestep
dictionary proportional to frame count. They isolate transport overhead with
small molecule graphs; they are not additional chemical trajectories. Imports,
fixture creation, input pre-reading and output fingerprints are outside the timer.
Naming includes conversion, global grouping, output, worker startup and cleanup.
The large-output case instead times conversion and incremental result hashing.

| Workload                                         | Base seconds | Candidate seconds | Change |
| ------------------------------------------------ | -----------: | ----------------: | -----: |
| fork, 64-frame metadata, 400 records             |       0.1841 |            0.1421 | −22.8% |
| fork, 100k-frame metadata, 400 records           |       1.9796 |            0.1416 | −92.8% |
| fork, 1M-frame metadata, 400 records             |      17.5767 |            0.1552 | −99.1% |
| spawn, 100k-frame metadata, 400 records          |       3.9656 |            0.5608 | −85.9% |
| forkserver, 100k-frame metadata, 400 records     |       6.7653 |            0.9324 | −86.2% |
| single process, 100k-frame metadata, 400 records |       0.0542 |            0.0539 |  −0.6% |
| existing fork batch path, 1,024 records          |       0.1405 |            0.1511 |  +7.6% |
| fork, two large returned frame arrays            |       1.9102 |            0.1500 | −92.1% |
| real 5,001-frame window, first 400 records, fork |       0.2951 |            0.1498 | −49.3% |
| same real window, all 4,509 records, spawn       |       4.7304 |            1.5893 | −66.4% |

The existing batch control is 10.6 ms slower in median despite identical task
counts and bytes; no batch-path speedup is claimed. The short candidate ranges
from 0.1328 to 0.2714 seconds, including one slower sample. Forkserver startup also
varies: base 5.3810–6.9459 seconds, candidate 0.8130–1.5321 seconds. The real cases
replay fixed no-HMM detection records, not DETECT or an entire production run.
Windows/macOS and full production trajectories have not been measured.

[All 60 timing samples](smiles_record_context_timings.csv) retain controls, CPU,
separate parent/child RSS peaks, source hashes and output fingerprints. RSS includes
imports/preparation and is not a concurrent process-tree total or a memory limit.

## Actual task payloads

[Twelve separate audit runs](smiles_record_context_audits.csv) wrap the production
serialization function and count the bytes actually returned by its protocol-5
ForkingPickler. They do not serialize each input a second time. The additional
initializer diagnostic serializes only the converter initializer and its args;
it excludes shared runner synchronization and process-startup framing.

| Metadata frames, 400 records | Largest base task bytes | Largest candidate task bytes | Total base / candidate task bytes |
| ---------------------------- | ----------------------: | ---------------------------: | --------------------------------: |
| 64                           |                   7,449 |                        1,021 |               2,977,218 / 406,018 |
| 100,000                      |                 744,849 |                        1,025 |             297,936,196 / 406,596 |
| 1,000,000                    |               9,747,882 |                        1,025 |           3,899,149,396 / 406,596 |

The candidate converter initializer is 1,245 bytes for each case, installed once
per worker. The slight short-record size difference comes from frame values in
the records themselves. Spawn/forkserver reproduce the same 100k-frame task
bytes. The existing batch control retains 16 tasks and 968,274 input bytes on both
versions. These are application call payloads, excluding outer pool framing.

## Correctness

The full non-GUI regression suite passes: 749 passed, 6 GUI deselected, 6 xfailed
and 3 non-strict XPASS. Forty-five new cases exercise actual fork/spawn/forkserver
pools, all three `miso` modes, parent-side timeline filters, failed SMILES/VF2
fallback and representative counts, custom extensions, and initializer/task/exit/
serialization/consumer failures. Failure tests keep exception tracebacks alive
while checking worker cleanup. Unrelated unpicklable RNG state never reaches
standard workers.

All 72 timing/audit outputs have equal ordered fingerprints within each
configuration. Downstream integration covers prepared HMM→MATRIX streams with
100k frames under fork/miso0, spawn/miso1 and forkserver/miso2, a 64-frame fork
control, and a real 101-frame DETECT→MATRIX window. All ten HDF5 outputs pass the
native validator; five base/candidate pairs preserve exact species labels,
atom/bond identities, ranges, populations, counts and transition evidence after
normalizing internal IDs. [Integration timings](smiles_record_context_integration.csv)
are single correctness samples and do not establish whole-pipeline speedups.

Repository pre-commit hooks, configured `ty==0.0.15`, and the existing SMILES audit
benchmark smoke check pass. Local tests reuse the official worktree's already
built ignored HTML bundle and native extension; these files are not PR contents.

## Reproduce

Use this checkout's script with separately installed base/candidate environments
that have matching dependencies. Replace interpreter, input and cache paths.
Run each configuration sequentially in fresh processes, alternating AB/BA/AB.
Use an available CPU set; timings here used `taskset -c 0-7` on Linux.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0
export MPLCONFIGDIR=/path/to/writable/matplotlib-cache

"$RNG_BASE_PYTHON" benchmarks/smiles_record_context.py --method fork \
    --frames 100000 --records 400 --nproc 4 --miso 1 --json base.json
"$RNG_HEAD_PYTHON" benchmarks/smiles_record_context.py --method fork \
    --frames 100000 --records 400 --nproc 4 --miso 1 --json candidate.json
```

Use `--frames 64` / `1000000`, `--method spawn` / `forkserver`, `--nproc 1`,
`--records 1024`, or `--large-output --miso 0` for controls. Repeat in separate
output files with `--audit` for payload measurements; exclude audit elapsed times
from speed summaries. Hashes and resource usage are emitted alongside timings.

Prepare a real detection stream once with `isomer_candidates.py --prepare --input TRAJECTORY --atom-types C Cl H O --nproc 1 --state STATE`, using the actual type map.
Then pass `--state STATE --frames 5001 --records 400` for the ordered prefix, or
`--records 100000 --method spawn` to select all records in the measured window.
`--records` is capped at the prepared state's count. This fixture has 2,600 atoms
and 4,509 records; its raw trajectory and prepared streams are not included here.

The prepared synthetic integration can be replayed with both interpreters:

```sh
"$RNG_HEAD_PYTHON" benchmarks/hmm_signal_pipeline.py --frames 100000 \
    --pairs 20 --nproc 4 --method spawn --miso 1 --output head-output --json head-run.json
"$RNG_HEAD_PYTHON" benchmarks/equal_z_compare.py base-output head-output
```

Create `base-output` with the same first command and the base interpreter. Use
the method/miso/frame combinations in the integration CSV for the other checks.
The real window uses `equal_z_candidates.py --mode pipeline --input TRAJECTORY --atom-types C Cl H O --nproc 4 --miso 1 --output OUTPUT --json RESULT`, then the
same independent comparison helper. Neither benchmark starts MD or training.
