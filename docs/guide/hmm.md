# Enable HMM filter

One can remove `--nohmm` in the command to enable HMM filter. However, HMM filter may also filter some important structures. One should read the original ReacNetGenerator paper{cite:p}`Zeng_PhysChemChemPhys_2020_v22_p683` and carefully adjust HMM parameters.

HMM filter should not be used to generate the special result by time step.

## Reproducibility and overlapping molecules

In the default mode, HMM filters each raw structure independently. With
`--miso 1` or `--miso 2`, it filters observations together according to the
[isomer policy](species.md#filtering-merging-and-trajectory-evidence).
Two independently filtered molecules sharing atoms can be present in the same
frame, even when their original observations were mutually exclusive.

With HMM enabled, retained records are ordered lexicographically by their atom
indices, bond endpoint pairs, and bond orders before assigning molecule IDs.
These numeric structure fields are already sorted by detection. Ordering uses
a disk-backed index and result spool; it does not retain all signals in memory.
In the default mode, temporary disk usage grows with the retained results.
Pooling instead spools the input structures before filtering. Each worker
processes the variants of one atom set together. The in-flight task limit
bounds the number of atom sets. Memory per task depends on the size of that
set's variant records and the trajectory length.

At an overlapping atom/frame cell, the largest molecule ID wins, regardless of
record arrival or write order. This is a deterministic tie-breaker, not a chemical
confidence score. The conflict flag remains set, and reaction-event inference
continues to exclude conflicting transitions. HDF5 molecule ranges describe the
effective presence of each raw structure, including overlaps between
independently filtered molecules.

This rule makes HMM atom routes independent of worker scheduling and worker
count. It can change molecule IDs, species ordering, and route positions compared
with older releases, including single-worker runs. Route positions count assigned
(nonzero) frames, rather than raw trajectory frames or physical timesteps.
