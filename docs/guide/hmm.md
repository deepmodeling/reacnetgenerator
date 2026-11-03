# Enable HMM filter

One can remove `--nohmm` in the command to enable HMM filter. However, HMM filter may also filter some important structures. One should read the original ReacNetGenerator paper{cite:p}`Zeng_PhysChemChemPhys_2020_v22_p683` and carefully adjust HMM parameters.

HMM filter should not be used to generate the special result by time step.

## Reproducibility and overlapping molecules

HMM filters each molecule independently. Two filtered molecules sharing atoms
can therefore be present in the same frame, even when their original
observations were mutually exclusive.

With HMM enabled, retained records are ordered lexicographically by their atom
indices, bond endpoint pairs, and bond orders before assigning molecule IDs.
These numeric structure fields are already sorted by detection. Ordering uses
a disk-backed index and result spool; it does not retain all signals in memory.
Temporary disk usage grows with the retained compressed results and structures,
while working memory still includes one complete molecule result.

At an overlapping atom/frame cell, the largest molecule ID wins, regardless of
record arrival or write order. This is a deterministic tie-breaker, not a chemical
confidence score. The conflict flag remains set, and reaction-event inference
continues to exclude conflicting transitions. HDF5 molecule ranges still describe
each molecule's independent HMM presence, including overlaps; they do not encode
the winner of this tie-breaker.

This rule makes HMM atom routes independent of worker scheduling and worker
count. It can change molecule IDs, species ordering, and route positions compared
with older releases, including single-worker runs. Route positions count assigned
(nonzero) frames, rather than raw trajectory frames or physical timesteps.
