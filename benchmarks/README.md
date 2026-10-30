# Molecule presence

`molecule_presence.py` compares filtering, atom-frame matrix construction, and
optional HDF5 molecule output using whichever ReacNetGenerator version is installed
in the invoking Python environment. It works with the official baseline
`97efe9ff3c6a9985acce96ff501575089d23c7f0` and this change.

These are synthetic occupancy records for the changed stages, not chemical
trajectory benchmarks or measurements of the complete analysis pipeline.

## Workloads and measurements

| Case   | Atoms |  Frames | Molecule records | Presence                                   |
| ------ | ----: | ------: | ---------------: | ------------------------------------------ |
| small  |     8 |      16 |               16 | Two consecutive half-trajectories per atom |
| sparse |   128 | 131,072 |              128 | One observation every 1,000 frames         |
| dense  |   512 | 131,072 |                1 | All atoms present throughout               |
| mixed  |   128 | 131,072 |              128 | 32 dense records and 96 sparse records     |

Each fresh process emits JSON containing:

- Wall time for filtering, matrix construction, and optional timeline writing.
  `wall_s` is their sum. Fixture creation, semantic hashing, and HDF5 validation
  are excluded from these timings.
- CPU time summed over the parent and completed children for those stages.
- Peak parent RSS and maximum child RSS in MiB, using `getrusage`. These include
  interpreter imports, fixture preparation, mapped pages, and validation. They
  are separate process peaks, not aggregate concurrent memory consumption.
- Additional and total temporary-file bytes after filtering. Shared paths are
  counted once; additional bytes exclude the original molecule records. These
  values exclude matrix mappings, the final HDF5 file, and transient worker spools.
- Matrix, conflict, and optional presence-range digests. Molecule IDs are mapped
  back to their input definitions before hashing so worker scheduling does not
  produce false mismatches. The public HDF5 semantic-manifest API is unchanged.

The implementation bounds matrix selections to 1,048,576 cells and presence
scans to 65,536 rows. This bounds temporary matrix selections, not total analysis
memory: matrices, individual signals, and retained trajectory metadata still
scale with the input.

## Running a comparison

Use separate baseline and candidate environments with identical Python and
dependency versions. Run this same script with both interpreters from the
candidate checkout. For example:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

"$RNG_BASE_PYTHON" benchmarks/molecule_presence.py --case dense --nproc 1
"$RNG_HEAD_PYTHON" benchmarks/molecule_presence.py --case dense --nproc 1
```

Set `RNG_BASE_PYTHON` and `RNG_HEAD_PYTHON` to the respective environment
interpreters. Add `--timed` for HDF5, `--hmm` for effective HMM signals, or
`--nproc 2` for the parallel fallback. Repeat each configuration at least five
times, alternating which revision runs first. Compare all semantic digests,
not timings alone. RSS collection requires a platform with the Python
`resource` module, such as macOS or Linux.

## Correctness checks

- The non-GUI suite completed with 405 passed, 8 skipped, 6 deselected,
  6 expected failures, and 3 non-strict unexpected passes. The skips were
  fork-only tests; the deselected tests matched `gui`.
- The final worker-payload change was followed by 33 passing presence/pipeline
  tests with the default macOS spawn method, and 112 passing presence, indexed
  Step 3, and utility tests with fork. A real two-worker test forces dense
  records to finish in reverse order and verifies both matrix and HDF5 output.
- Eight complete runs against the fixed official baseline produced identical
  text hashes and, where enabled, public HDF5 semantic manifests. Configurations
  covered HMM on/off, origin retention, all three isomer modes, frame/timestep
  filters, split/selectatoms, species-only execution, and ordinary text output.
- Ruff 0.17.2, Velin at tag 0.0.12, and ty 0.0.15 passed. Native Linux, Windows,
  and GUI testing were not performed locally.

The shared presence reader also checks missing/extra signals, record counts,
truncated definitions, and signal length. Regression tests cover empty records,
duplicate indices, overlapping molecules, write-batch boundaries, failed dense
expansion, shared-file cleanup, scratch-directory fallback, and retained
metadata for HMM-filtered molecules.

## Recorded comparison: 2026-10-30

Apple M4, macOS 15.6.1, CPython 3.12.7; NumPy 2.5.3, h5py 3.14.0,
hmmlearn 0.3.3, LZ4 4.4.5, pandas 3.0.6, RDKit 2026.3.6, and Open Babel
3.2.1. The complete installed distribution sets, including development tools,
were compared and differed only in the ReacNetGenerator revision. BLAS thread
limits were set as above.

All numbers below are medians from five alternating baseline/candidate pairs.
Across 13 configurations and 130 fresh-process runs, every matrix/conflict
digest and every enabled HDF5 range digest matched. Arrows mean baseline to
candidate; times are seconds, RSS is MiB, and additional temporary size is bytes.

### No HMM, one worker, no HDF5

| Case   | Stage wall time | Stage CPU time | Matrix time     | Parent peak RSS | Additional temporary bytes |
| ------ | --------------- | -------------- | --------------- | --------------- | -------------------------- |
| small  | 1.402 → 0.014   | 1.393 → 0.005  | 0.0091 → 0.0093 | 249.1 → 253.8   | 14,624 → 0                 |
| sparse | 1.509 → 0.056   | 1.454 → 0.016  | 0.0627 → 0.0507 | 328.7 → 333.2   | 283,738 → 0                |
| dense  | 1.694 → 0.065   | 1.647 → 0.030  | 0.2758 → 0.0588 | 519.1 → 397.0   | 528,864 → 787              |
| mixed  | 1.604 → 0.093   | 1.526 → 0.045  | 0.0832 → 0.0715 | 386.6 → 360.9   | 12,854,519 → 25,216        |

The no-HMM stage gains include avoiding worker startup: the baseline creates a
pool even for one worker. Matrix-only timings separate that effect from the
bounded-write improvement. Small and sparse parent RSS did not decrease; this
change does not promise lower total memory for every workload.

### HDF5 enabled

| Case   | HMM | Workers | Stage wall time | Matrix time     | Parent peak RSS |
| ------ | --- | ------: | --------------- | --------------- | --------------- |
| small  | off |       1 | 1.423 → 0.018   | 0.0099 → 0.0092 | 252.9 → 258.3   |
| sparse | off |       1 | 1.643 → 0.186   | 0.0603 → 0.0540 | 334.4 → 336.7   |
| dense  | off |       1 | 1.816 → 0.181   | 0.2801 → 0.0603 | 519.3 → 396.3   |
| mixed  | off |       1 | 1.987 → 0.256   | 0.0869 → 0.0673 | 389.7 → 352.1   |
| small  | on  |       1 | 1.696 → 1.699   | 0.0117 → 0.0112 | 253.3 → 259.1   |
| sparse | on  |       1 | 2.287 → 2.280   | 0.0483 → 0.0511 | 323.9 → 334.7   |
| dense  | on  |       1 | 2.231 → 1.918   | 0.3342 → 0.0641 | 519.8 → 394.5   |
| mixed  | on  |       1 | 2.475 → 2.400   | 0.0706 → 0.0606 | 389.6 → 390.8   |
| mixed  | off |       2 | 2.613 → 2.099   | 0.0981 → 0.0633 | 381.5 → 354.4   |

HMM filtering itself is unchanged. Small/sparse HMM stage timings stayed near
the baseline, while the dense matrix benefited from bounded slice writes.
These results isolate the implemented stages and must not be extrapolated to
the speedup of a complete scientific analysis.
