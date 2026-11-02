# Immutable native participant metadata

Native timeline writing repeatedly searches the same immutable molecule-ID
bounds, inverts the species-name dictionary, and reads atom/bond definitions for
event participants. The writer now indexes species names once, retains up to
4096 ID bounds and 4096 windows of at most 64 intervals each, and caches atom/bond arrays with both an 8 MiB payload limit and
an 8192-entry limit. Oversized definitions bypass the array cache. Temporary read
arrays, Python object headers and the species-name index are additional; this
is not a process-memory limit. All caches belong to one open file.

Every event still checks presence at its actual frames, disjoint participants,
atom conservation, species consistency and inferred bond changes. Validation
creates a fresh range lookup over stored bytes, respecting a smaller caller
block size. The original scalar presence search remains the unit-test oracle.
No schema, output ordering, worker policy or scientific identity rule changes.

## Measurements

Baseline is official `ffefc395caaf0318139ea21959dc12d27ca837d2`. The three changed
production-module hashes identify the measured candidate in the
[24 recorded runs](native_participant_metadata_timings.csv). All use Python
3.12.12, NumPy 2.5.4, SciPy 1.18.1, OpenBabel 3.2.1, h5py 3.14.0, RDKit
2026.9.1, ASE 3.29.0, hmmlearn 0.3.3, LZ4 4.4.5 and NetworkX 3.7. The shared
native extension SHA-256 is
`ad2ed1dd40691fc73714e45907cc2b599a8e048cdec5ded37549c8c5a9bafd8c`.

The real LAMMPS dump contains 2600 atoms. The initial 201-frame window hash is
`0e3ce09a53bd0e489b4dbc4df8116a2c03371890527d72c35af0e3197ec1e537`; raw production inputs are not published here.
The two 5001-frame controls use initial and middle trajectory windows. They are
windows of one simulation, not independent scientific replicates.

Each timing is a fresh process, with 32 requested workers, fork, CPUs 0–31 and
one native math thread. Three-way dense-window order is official/candidate/fork,
candidate/fork/official, fork/official/candidate. Two-way controls alternate
AB/BA/AB. Timers cover DETECT through MATRIX, including worker startup/cleanup;
imports, initialization, input pre-reading and validation are outside timing.
Separate acceptance tasks run on CPUs 32–47, so these are same-host observations,
not isolated-machine or full-production guarantees.

| Full-pipeline window                 | Official median (s) | Candidate median (s) |
| ------------------------------------ | ------------------: | -------------------: |
| 201 frames, miso=0, HMM off          |              55.520 |               20.889 |
| 5001 initial frames, miso=1, HMM off |              53.576 |               51.094 |
| 5001 middle frames, miso=1, HMM on   |              33.256 |               33.898 |

The dense window improves 2.66 times (62.4% less time). The sparse-HMM control
observes 1.9% more time. All individual samples and control variation are retained
in the CSV.

Separate 201-frame audits sampled all descendants through `/proc` every 0.2 s.
Complete snapshots observed peak summed PSS of 1242.2 MiB for official and
1523.7 MiB for candidate (119 and 46 complete samples; up to 33 processes).
Short-lived workers can exit during a scan: the 68/87 transient read errors are
recorded and partial snapshots are excluded from the peak. Sampled PSS is an
observation and can miss an instantaneous maximum. Audit timings are excluded
from the three timing medians.

The deployed legacy fork is an additional context control: its dense-window
median is 26.989 seconds and observed PSS is 1002.3 MiB. It uses another native
schema and established identity/event conventions (82,542 events versus 82,001),
so it is not the upstream output-equivalence target. Candidate is faster on
this small window but uses more memory than that fork. This does not establish
a speedup on a complete production trajectory.

## Correctness

The full local suite reports 769 passed, 5 skipped, 6 xfailed and 3 non-strict
XPASS. All nine matched official/candidate timing pairs pass native validation
and agree on 12 normalized fingerprints covering frame mappings, exact raw
atom/bond definitions, names, ranges, populations, reaction counts and instance
evidence. Each dense result has 82,001 events. Memory-audit outputs are separate
observations and are not included in that nine-pair validation count.

New tests cover cache bounds and eviction, oversized/empty metadata, file-local
lifetime, cleanup, retryable I/O errors, invalid species IDs, and warm-cache
presence gaps, atom overlaps and conflicting bonds. The standalone benchmark
also completes a two-frame real-input smoke run and native semantic manifest.
Ruff 0.17.2 and the configured ty 0.0.15 pass.

## Reproduction

Use the same scientific environment and native extension for both source trees.
For example, with each `PYTHONPATH` pointing at one checkout and a fresh output
directory for every run:

```sh
PYTHONPATH=/path/to/source taskset -c 0-31 /path/to/python \
    benchmarks/native_participant_metadata.py \
    --trajectory /path/to/input.dump --atom-names C Cl H O \
    --frames 201 --nproc 32 --miso 0 --output-dir /path/to/fresh-result
```

[The script](native_participant_metadata.py) streams only the requested initial
window, refuses an existing output directory, records source/input hashes and
dependencies, and validates the output separately using a semantic manifest.
Use `--frames 5001`, and `--miso 1 --hmm` for the HMM control. This script does
not collect concurrent process-tree PSS; that requires a separate monitor.
Linux timings above do not establish cross-platform performance.

## Final interval windows and large-window acceptance

Range lookups retain the existing 1 MiB shared numeric blocks and at most 4096
molecule-ID bounds, plus 4096 windows of up to 64 intervals. Window arrays add
at most 4 MiB. A window answers only covered frames; unknown edges fall back to
the stored search, preserving inclusive endpoints and gaps. Fresh validator
state reads stored bytes independently. Smaller caller block sizes also limit
window reads. Cache budgets exclude Python headers and temporary arrays.

The final 201-frame observed PSS values are 1242.2 MiB official, 1523.7 MiB
candidate and 1002.3 MiB legacy. This small-window candidate observation is
higher; it is not evidence of a process-memory reduction.

One additional 5001-frame, miso=0/HMM-off pair used the same 32 requested
processes and CPUs 32–47 for both candidate and deployed fork:

| Version     | Pipeline seconds | Observed peak tree PSS (MiB) |
| ----------- | ---------------: | ---------------------------: |
| Candidate   |          459.054 |                       2748.7 |
| Legacy fork |          708.967 |                       2732.7 |

This observes 35.3% less time and 0.6% more PSS. It is one pair with concurrent
acceptance work, not a repeated whole-production speed claim. Independent native
validation and indexing are outside timing. A full 250001-frame production
trajectory has not been rerun.

The final large output passes full native validation, and its 2,239,475 event
instances agree with unmodified official output on all 12 normalized scientific
fingerprints. Interpretation metadata and dataset shapes/dtypes also agree.
Internal species/molecule IDs and row order can change between runs; raw HDF5
array hashes are not an equivalence criterion. Actual trajectory atom indices,
bond orders and times remain exact.
