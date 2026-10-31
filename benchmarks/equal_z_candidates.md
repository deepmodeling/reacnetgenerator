# Periodic candidates with component-local z-order checks

Baseline: official `master` at `8dced2340189db52bd056bb64e36990b2255a629`.
The candidate changes only `reacnetgenerator/_detect.py` in production code.

Rounded coordinates often give atoms in different molecules exactly the same z
coordinate. The baseline abandons periodic candidate acceleration for the entire
frame whenever any such tie exists. The candidate checks connectivity before
deciding whether the tie can affect bond insertion.

## Compatibility boundary

The accepted geometric candidate edges form an upper bound on the connectivity
that bond insertion and cleanup can retain. If tied atoms belong to different
candidate components, their unspecified global z-sort order cannot change the
relative candidate insertion order inside either component. The existing native
Open Babel fallback remains when a component contains tied atoms, including ties
connected only through intermediate atoms or across a periodic boundary.

The component check runs only for frames with z ties. It uses a sparse graph and
component labels, adding storage proportional to atoms plus accepted candidate
edges; it does not construct a dense atom-by-atom matrix. Candidate generation
retains its existing input-dependent cost. This is not a total-RSS bound.

Tests compare exact molecule records with the original Open Babel implementation.
They cover disconnected H/H2/O2 cases, isolated atoms, shuffled aromatic rings,
overvalent cleanup, zero-valence elements, phosphorus insertion rules, indirect
ties, periodic ties, and dispatch through the real 2048-atom threshold. Existing
intra-component tie fallbacks remain covered.

## Reproduction

Use the same dependency versions and native extension for base and candidate,
with separately resolved source trees. The measurement JSON records the actual
`_detect.py` hash and dependency versions. This campaign used Linux x86-64,
an AMD EPYC 7K62 host with 48 visible CPUs (one thread per core), Python 3.12.12,
Open Babel 3.2.1, h5py 3.14.0, NumPy 2.5.4, SciPy 1.18.1, and
RDKit 2026.9.1. The full dependency set was identical between variants.

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python benchmarks/equal_z_candidates.py --case independent --iterations 3 --json independent.json
python benchmarks/equal_z_candidates.py --case unique --iterations 3 --json unique.json
python benchmarks/equal_z_candidates.py --case connected --iterations 3 --json connected.json
python benchmarks/equal_z_candidates.py --case rounded --iterations 3 --json rounded.json
python benchmarks/equal_z_candidates.py --case small --iterations 3 --json small.json
```

Each kernel invocation first verifies exact molecule records against original
Open Babel, then times three calls through the production dispatch entry point.
Five fresh processes per variant/case run in alternating AB/BA order. The cases
contain 2048 atoms, except the 64-atom small-input control. `independent` places
z ties across separate H2 components; `unique` perturbs those z values;
`connected` introduces a tie inside one pair and must keep the fallback;
`rounded` uses fixed-seed mixed elements with rounded z coordinates.

Pipeline mode accepts a LAMMPS dump without copying it:

```sh
python benchmarks/equal_z_candidates.py --mode pipeline --input INPUT.dump \
    --atom-types C Cl H O --nproc 32 --miso 1 --output OUTPUT --json run.json
python benchmarks/equal_z_compare.py BASE_OUTPUT CANDIDATE_OUTPUT --json comparison.json
```

The real-input campaign uses two distinct contiguous windows from one trajectory:
each has 5001 frames and 2600 atoms. HMM is disabled, `miso=1`, with 32 workers
bound to CPUs 0-31; the input is read before each run to warm the file cache.
Each variant/window runs three times in alternating order. They are computational
repeats, not independent physical trajectories. A separate 101-frame preflight is
excluded from performance medians. Raw trajectory coordinates are not included
in this repository; the synthetic cases above are fully reproducible.

Pipeline times include DETECT through MATRIX, including normal outputs and the
schema-1.1 HDF5 timeline. Imports and post-run validation are excluded. Every
timeline passes the official validator. The comparison helper retains exact
species labels, atom membership, bond orders, effective presence, aggregate
counts, and transition participants/bond changes while excluding internal ID
assignment and table ordering. It does not equate different SMILES or formulas.
The ordinary official manifest is also retained; parallel DETECT may legally
change its ID-dependent column ordering.

RSS values are separate parent and largest-child process high-water marks and
must not be added into a simultaneous process-tree peak. Kernel RSS includes its
reference-validation preflight; pipeline RSS is captured before post-run
validation. Values reported by this Linux benchmark are in KiB.

## Results

All 50 kernel samples matched exact Open Babel molecule records. All 14 pipeline runs passed native HDF5 validation; all seven base/candidate pairs passed the published ID-independent comparison, including exact species labels and bond orders. The 101-frame preflight is excluded from the median tables.

| Kernel case | Base median seconds [min, max] | Candidate median seconds [min, max] | Wall change |
| ----------- | -----------------------------: | ----------------------------------: | ----------: |
| independent |  0.226695 [0.226213, 0.227033] |       0.031318 [0.031186, 0.031669] |      -86.2% |
| unique      |  0.031737 [0.031035, 0.032216] |       0.031352 [0.031140, 0.031718] |       -1.2% |
| connected   |  0.227029 [0.226350, 0.227298] |       0.228910 [0.228371, 0.229773] |       +0.8% |
| rounded     |  0.233835 [0.233471, 0.234110] |       0.022417 [0.022314, 0.022715] |      -90.4% |
| small       |  0.000785 [0.000774, 0.000792] |       0.000788 [0.000780, 0.000796] |       +0.4% |

The retained-fallback `connected` control is 1.9 ms slower in median (+0.8%). The small-input difference is about 3 microseconds (+0.4%); unique-z timings overlap. These controls do not establish universal non-regression.

| Real window, 5001 frames / 32 workers | Base median seconds [min, max] | Candidate median seconds [min, max] | Wall change |
| ------------------------------------- | -----------------------------: | ----------------------------------: | ----------: |
| early                                 |        94.806 [93.674, 96.029] |             49.410 [48.689, 49.950] |      -47.9% |
| middle                                |     103.979 [103.875, 105.084] |             57.224 [56.940, 57.247] |      -45.0% |

| Window | Step 1 median seconds, base → candidate | Candidate-eligible frames, base → candidate |
| ------ | --------------------------------------: | ------------------------------------------: |
| early  |                         66.460 → 21.408 |                         50/5001 → 4290/5001 |
| middle |                         66.403 → 19.832 |                         40/5001 → 4401/5001 |

Eligibility counts were collected separately, after timing, by passing all original frames through the unchanged parser and each variant's candidate-selection method. The candidate retains 711 early-window and 600 middle-window fallbacks. Six early-window frames previously sensitive to an unconditional tie-break also remain on the native fallback; their exact bond orders and presence intervals match the official baseline.

### Validation and limits

- Final non-GUI test suite: **635 passed, 6 deselected, 6 xfailed, 3 non-strict XPASS** on Linux/Python 3.12.12.
- Repository pre-commit checks and pinned `ty==0.0.15` passed. The wheel built successfully; an isolated install passed dependency checks and the independent-component dispatch smoke test, with its installed production source matching the benchmarked hash.
- The public comparison helper also rejects deliberate valid-schema timestep/bond-order changes and a changed species-population text entry; those negative controls are outside all timings.
- Three repeats per real window showed identical compared evidence within each variant. Ordered native manifests can differ because parallel detection assigns internal IDs in completion order.
- The same 64 pinned dependencies and the same newly built native `dps` extension were used for both frozen source trees. Only `_detect.py` differed in production sources.
- No other task-owned compute workload overlapped timed samples. Validation and eligibility audits ran afterward. The host was shared, and these samples are not a full-trajectory or cross-platform performance guarantee.
- The real windows are from the same 2600-atom trajectory; this is computational performance evidence, not a claim about chemical ground truth or independent replicates.

Production file hashes:

- Baseline: `f61cf806bbc34c31a75e437ee993751ea4d919e69ecd586a3d1e11ee67fd28ae`
- Candidate: `f4b86168e1ed23a3851ed215dc34640b62b029a75083f21c9f6de7bf6dfddc75`

All sample measurements, including unfavorable controls and the preflight, are in [equal_z_candidates.csv](equal_z_candidates.csv).
