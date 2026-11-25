# HDF5 timeline (schema 1.2)

`--show-molecule-time` and `--reaction-event` also write one compact HDF5
timeline when the PATH stage runs. Without `--output-dir`, its default path is
`<input>.timeline.h5`; with `--output-dir DIR`, it is `DIR/timeline.h5`.

Use `--timed-output FILE` to choose another path or to request a timeline
without enabling either legacy CSV switch:

```sh
reacnetgenerator -i trajectory.bond --type bond -a H He --nohmm --timed-output timeline.h5
```

In Python, pass `timed_output="timeline.h5"` to `ReacNetGenerator`, or through
`reacnetgenerator.run`. The returned artifact map includes `timeline`. Relative
paths resolve against the working directory; this explicit path is not relocated
by `output_dir`. Its parent directory must exist. A species-only `run_items`
request cannot produce a timeline and raises an error if combined with a timed
switch or `timed_output`.

The timeline contains molecule definitions and effective presence intervals,
aggregate reaction events, instance-level participants and inferred bond changes,
and frame/source/configuration metadata. It does not contain coordinates. It does
not change the existing text outputs. The legacy CSV switches continue to write
their CSV files during this compatibility period.

## Why validation is part of the format

A saved analysis result is useful beyond the process that wrote it only when a
user or another program can establish what it contains and whether its tables
still agree. Opening an HDF5 file proves neither. The timeline validator makes
schemas 1.0 through 1.2 a verifiable software boundary:

- CI and scientific regression tests can reject a result with broken offsets,
  dangling IDs, nonmaximal ranges, inconsistent reaction totals, or instance
  evidence that disagrees with the aggregate counts even when the file remains
  readable.
- Collaborators and archives can record a compact semantic manifest and compare
  results without depending on HDF5 compression, chunking, or local file paths.
- Downstream consumers, including visualization and analysis applications, can
  validate producer output before interpreting it and discover the contract from
  an installed JSON descriptor.
- Future schema migrations have a concrete 1.0 baseline against which changed
  structure and meaning can be reviewed.

This is an artifact-integrity and reproducibility check. It does not establish
chemical truth, transition states, barriers, kinetics, or agreement with the raw
trajectory. Source paths, sizes, and modification times are provenance hints;
they are not content hashes.

## Meaning of time and identity

- A frame is a zero-based **analyzed** frame. `frames/source_id` and
  `frames/source_frame` locate it in an input occurrence and the zero-based frame
  within that file, before `stepinterval` sampling. Sampling spans the concatenated
  input sequence. Repeated filenames have separate source IDs.
- `frames/timestep` preserves the parser's timestep value. XYZ/extxyz currently
  use the analyzed frame number as their timestep convention; no physical time
  unit is implied. Repeated or decreasing timestep values are allowed.
- Molecule IDs start at **1**, matching the internal atom-frame matrix. Species,
  reaction-type, source, atom-type and atom indices start at **0**. Dictionary IDs
  are local to a file, not stable scientific identifiers.
- Atom indices refer to RNG's canonical atom ordering, not arbitrary original
  trajectory atom IDs. The existing LAMMPS parser expects contiguous one-based
  source IDs; this feature does not extend that parser contract.
- Molecule ranges are closed `[start_frame, end_frame]` intervals of the signal
  used to construct the atom-frame matrix. HMM defines this signal when enabled.
  With HMM off, ranges follow bond hysteresis when configured or ordinary
  frame-local detection otherwise. The existing reaction conflict rules handle
  overlapping molecule signals. Ranges do not assert unique atom ownership.
- Legacy molecule CSV describes **observed** frames. With HMM enabled its rows
  need not match these effective ranges. CSV frame/timestep filters do not filter
  HDF5 ranges or events; their values remain recorded in configuration.
- An event at `transition=t` describes `t -> t+1`. `count` aggregates occurrences
  of one reaction type at that transition. Reaction recognition and species
  cancellation follow the existing ReacNetGenerator algorithm.

## Versioned file contract

Root attributes:

| Attribute               | Value/meaning                                                                                      |
| ----------------------- | -------------------------------------------------------------------------------------------------- |
| `format`                | `reacnetgenerator-timeline`                                                                        |
| `schema_version`        | UTF-8 string `1.2`                                                                                 |
| `status`                | `incomplete` during construction, `complete` on successful close                                   |
| `capabilities`          | JSON array `["molecule_ranges", "reaction_events", "transition_evidence"]`                         |
| `rng_version`           | Installed RNG version                                                                              |
| `created_utc`           | ISO 8601 UTC timestamp                                                                             |
| `configuration`         | JSON from `parameter_provenance()`: normalized constructor parameters and explicit parameter names |
| `atom_index_convention` | `zero-based RNG canonical atom order`                                                              |
| `molecule_range_basis`  | `HMM signal`, `observed signal`, or (in 1.2) `hysteresis signal`                                   |

Configuration includes HMM settings, sampling, cell/PBC options, detection
backend/cutoffs, species identification and selection settings. A null configured
cell means that the parser derives it from input; actual frame cell matrices are
not copied into this format. Input paths, sizes and modification times provide
provenance, not a content-integrity guarantee. No input hash is computed by default.

Use `--bond-hysteresis FORM BREAK` or the Python argument
`bond_hysteresis={"form": ..., "break": ...}` to track each pair across analyzed
frames. State starts with all pairs absent and continues across ordered input files
and `--split` partitions. Frames skipped by `stepinterval` do not update state.
Run independent trajectories separately. The opt-in detection pass runs in frame
order in one process. Later analysis stages retain the configured worker count.

- For coordinates in angstroms, the thresholds multiply the sum of ASE covalent
  radii. A pair forms at or below FORM and breaks at or above BREAK.
  Periodic inputs use the nearest image. For example, `--bond-hysteresis 1.25 1.35`
  specifies two distance thresholds. It replaces the single-cutoff ASE settings.
  Open Babel assigns orders on the retained graph, or `--use-ase` assigns order one.
  To adjust only a particular pair's break distance, add
  `--bond-break-mult H-H MULT` (repeatable for different pairs), or set
  `bond_hysteresis={"form": form, "break": outer, "pair_break": {"H-H": hh_outer}}`.
  Each override uses the same radius-sum units and must exceed the common formation
  multiplier. All other pairs keep the global outer multiplier. The optional map
  is recorded in the configuration. Older strict 1.2 readers reject files using it.
- For ReaxFF bond files the thresholds apply to the reported fractional bond
  order before rounding. A pair forms at or above FORM and breaks at or below
  BREAK (for example, `--bond-hysteresis 0.60 0.30`). Missing pairs break.
  When both directions report the same pair the larger value is used. The
  reporting cutoff in the simulation determines which orders are available.

The default HMM setting still applies. Use `--nohmm` to isolate bond hysteresis.
With both filters enabled the final ranges use `HMM signal`. The configuration
records the bond thresholds in either case. Readers predating the new enum can
reject files that use `hysteresis signal`.

If the source MD uses a fixed physical timestep, record it with
`--timestep-value 0.0001 --time-unit ps`, and record the coordinate unit with
`--length-unit angstrom`. The corresponding Python arguments are
`timestep_value`, `time_unit`, and `length_unit`. These optional values live in
`configuration.parameters`. Missing or `null` fields are unknown. This does
not alter the frame table, the capability list, or the schema version.

`timestep_value` is the positive finite duration of one source MD integration
step, before dump or analysis sampling. Multiply differences in `frames/timestep`
by this value only when those entries are source integration steps on one
consistent time axis. Use `fs`, `ps`, or `ns` for time and `angstrom`, `bohr`, or
`nm` for length. Supply a value only if it applies across the input files.

ASE and Open Babel bond perception interpret input distances in angstroms.
Convert coordinates and cell vectors to angstroms before coordinate-based
analysis, and use `--length-unit angstrom` to record their unit. For bond-file
input, the length metadata can describe the originating trajectory's unit.

All datasets are one-dimensional. Numeric columns are little-endian signed 64-bit
integers; text columns are variable-length UTF-8 strings. Columns in each ordinary
table have equal lengths. Dataset layout/chunking is not part of semantic identity.

| Group                 | Columns                                                                                        | Row meaning                                                                     |
| --------------------- | ---------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| `sources`             | `path` (text), `size_bytes`, `mtime_ns`                                                        | One input occurrence, in supplied order                                         |
| `frames`              | `source_id`, `source_frame`, `timestep`                                                        | One analyzed frame, in order                                                    |
| `atoms`               | `type`; `type_name` (text)                                                                     | `type` has one entry per atom; `type_name` is the separate atom-type dictionary |
| `species`             | `name` (text)                                                                                  | Unique species name                                                             |
| `species`             | `formula`, `ring_system` (optional text columns)                                               | Searchable descriptors aligned with `species/name`                              |
| `molecules`           | `species_id`                                                                                   | One molecule definition, row `molecule_id - 1`                                  |
| `molecules`           | `atom_offsets`, `atom_index`                                                                   | Offset-delimited atom lists                                                     |
| `molecules`           | `bond_offsets`, `bond_atom_index_1`, `bond_atom_index_2`, `bond_order`                         | Offset-delimited bonds with global canonical atom endpoints                     |
| `molecule_ranges`     | `molecule_id`, `start_frame`, `end_frame`                                                      | One maximal effective interval, ordered by molecule then start frame            |
| `reaction_types`      | `reactant` (text), `product` (text), `total_count`                                             | One unique pair of existing formatted reaction sides                            |
| `reaction_events`     | `transition`, `reaction_type_id`, `count`                                                      | One positive aggregate count, ordered by transition                             |
| `transition_evidence` | `transition`, `reaction_type_id`                                                               | One inferred connected reaction instance, ordered by transition                 |
| `transition_evidence` | `participant_offsets`, `participant_molecule_id`, `participant_side`                           | Offset-delimited molecule instances; side `0` is reactant and `1` is product    |
| `transition_evidence` | `bond_change_offsets`, `bond_atom_index_1`, `bond_atom_index_2`, `before_order`, `after_order` | Offset-delimited canonical bond differences for each instance                   |

Both offset columns start at zero and have `number_of_molecules + 1` entries.
For molecule ID `m`, its payload is `[offset[m-1]:offset[m]]`. The final offset
is the payload length. Bond endpoint/order columns have identical lengths. Empty
payloads and empty tables are valid. Reaction-type totals equal the sum of stored
event counts for that type. Reaction-side strings retain the existing output
notation; consumers should not assume splitting on `+` parses every possible
SMILES. Evidence participant and bond-change offsets also start at zero and end
at their payload lengths. Participants reference the existing molecule
definitions; their atom sets are disjoint within each side and conserved across
the transition. Each reactant is present at transition frame `t`, and each
product is present at `t + 1`, according to `molecule_ranges`. Regrouping
evidence rows by transition and reaction type exactly reproduces
`reaction_events/count`.

Schema 1.2 adds two independently optional UTF-8 columns, indexed by
`species_id`. Check for each dataset in the HDF5 file before reading it.
The reader and validator also accept schema 1.1 with transition evidence and
aggregate-only schema 1.0. Older readers that whitelist versions need an update
to open 1.2 files.

- `species/formula` counts recorded atoms using their configured type names.
  With carbon present, C and then H come first, followed by other labels in
  alphabetical order. Carbon-free formulas use alphabetical order throughout.
  A count of one is implicit. Only explicitly recorded hydrogens contribute.
- `species/ring_system` describes the undirected bond endpoints of a connected
  molecule. The writer computes descriptors when it first encounters each
  species name. With composition-only grouping (`miso=2`) it omits the entire
  ring column because one species ID can cover different graphs.
- Trees are `acyclic`. One cycle gives `monocyclic-N`. For a chain of edge-fused
  rings in a single outerplanar ring block, `fused-N-...` follows the interior
  faces in adjacency order, choosing the numerically larger sequence of the two
  directions (e.g. `fused-6-5-6`). Disjoint simple rings joined by bonds or
  chains give `linked-N-...`. Two simple rings sharing one atom give
  `spiro-N-M`. Linked and spiro sizes are in descending order. Remaining
  topologies, including bridged cages and branched fused systems, give
  `polycyclic`.

These fields support candidate filtering. Atom indices, bond order, aromaticity
and electronic state do not enter the ring label. Positional isomers and
chemically distinct ring systems can share both descriptors.

## Compact Python reading

```python
from reacnetgenerator.timedoutput import (
    read_metadata,
    iter_frames,
    iter_species,
    iter_molecules,
    iter_molecule_ranges,
    iter_reaction_types,
    iter_reaction_events,
    iter_transition_evidence,
)

metadata = read_metadata("timeline.h5")
for interval in iter_molecule_ranges("timeline.h5", block_rows=4096):
    print(interval.molecule_id, interval.start_frame, interval.end_frame)
for event in iter_reaction_events("timeline.h5", block_rows=4096):
    print(event.transition, event.reaction_type_id, event.count)
for evidence in iter_transition_evidence("timeline.h5", block_rows=4096):
    print(evidence.transition, evidence.participants, evidence.bond_changes)
```

`iter_frames`, `iter_molecules`, `iter_molecule_ranges`, `iter_reaction_events`
and `iter_transition_evidence` return frozen dataclass records. `iter_species`
yields `(species_id, name)`; `iter_reaction_types` yields
`(type_id, reactant, product, total_count)`. Definitions include `atom_index` and
`(atom1, atom2, order)` bond tuples. A `BondChange.kind` is `formed`, `broken` or
`order_changed`; the before/after orders remain available.
Source paths and atom-type tables can be accessed directly with h5py under the
public contract above.

## Validation and semantic comparison

Validate the complete structural and cross-table contract before consuming an
artifact:

```python
from reacnetgenerator.timedoutput import validate_timed_output

summary = validate_timed_output("timeline.h5")
print(summary.frames, summary.reaction_events)
```

The validator checks required attributes and datasets, strict JSON configuration,
agreement between `runHMM` and the declared range basis, one-dimensional column
types and alignment, frame mappings against the configured global sampling
stride, offset bounds, dictionary and frame references, molecule bond membership,
range ordering/maximality, event ordering and uniqueness, reaction-type totals,
evidence participant identity and atom conservation, exact bond differences, and
participant presence on the appropriate side of the transition, and
evidence-to-aggregate counts. Numeric scans use `block_rows`; one molecule
definition, evidence instance or dictionary can still exceed that working-memory
budget.

Repeated event participants are checked against immutable molecule tables after
the writer flushes them. The writer retains at most 8 MiB of atom/bond array
payload across 8192 cached entries; an oversized molecule bypasses that cache.
Python object headers and a newly read molecule are additional. Species names
are indexed once by species ID. Range lookups retain at most 1 MiB of numeric
blocks, 4096 molecule-ID bounds and 4096 windows of at most 64 intervals each.
Window arrays add at most 4 MiB; Python metadata is additional. A query outside
a cached window searches the stored ranges again, including gaps at its edges.
All caches
belong to one open file and are discarded when it closes; presence is still
checked for each event frame, including gaps between a molecule's ranges.

The validator creates a fresh range lookup over stored data and honors smaller
`block_rows` requests. These caches do not bypass atom conservation, participant
overlap, species, bond-change or aggregate-count validation, and do not change
the on-disk schema or event ordering.

Create and compare deterministic manifests in Python:

```python
from reacnetgenerator.timedoutput import (
    compare_semantic_manifests,
    semantic_manifest,
)

reference = semantic_manifest("reference.h5")
candidate = semantic_manifest("candidate.h5")
differences = compare_semantic_manifests(reference, candidate)
```

The default manifest hashes canonical dataset values, non-location configuration,
and interpretation metadata. HDF5 compression and chunk layout do not affect it.
It omits source paths/sizes/timestamps, creation time, RNG build version, and
path-valued configuration so that relocating an otherwise identical analysis does
not create a difference. Set `include_provenance=True` when those fields must also
match.
The manifest SHA-256 is a deterministic comparison key, not a signature or a
raw-input integrity guarantee.

The same operations are available to shell scripts and CI:

```sh
reacnetgenerator-check-timed-output timeline.h5
reacnetgenerator-check-timed-output timeline.h5 \
    --write-manifest timeline.manifest.json
reacnetgenerator-check-timed-output candidate.h5 \
    --compare-manifest timeline.manifest.json
```

Successful commands print a compact JSON summary and return 0. Invalid artifacts
or manifest differences are written to standard error and return 1; command-line
usage errors return 2. Manifest writes use an atomic sibling replacement.

`read_schema_descriptor()` returns the installed
`schemas/timed-output-schema.json` contract for tools that need to inspect the
schema without scraping this guide.

Numeric rows are read in blocks, ranges are not expanded into frame rows, and
counts are not expanded into individual events. Dictionary strings are read one
at a time. One molecule definition or one string can exceed the block budget.
Call `.close()` on an iterator when stopping early, or use `contextlib.closing`.
The compact readers perform local checks needed to read their requested table;
call `validate_timed_output()` when the complete artifact contract matters.

## Publication and resource behavior

The writer uses an exclusive temporary sibling ending in `.incomplete`. It
flushes and closes all datasets before `os.replace` publishes the destination.
A PATH failure preserves the prior destination and retains the incomplete sibling
for diagnosis. A replacement failure also preserves the destination; its sibling
may already have `status=complete`. Exceptions propagate to the caller.
Atomic replacement is not a guarantee of durability across power loss.

Each column buffers at most 8192 rows or about 1 MiB of scalar payload, except
that a single large string may exceed the byte budget. Existing decoded signal
records, species/reaction dictionaries, graphs and upstream atom-frame state still
scale with content. Evidence is returned one transition at a time and its columns
use the same bounded writer buffers; one connected reaction instance can still
scale with its participating molecules. This feature does not promise constant
total memory. HDF5 is written only in the parent process; worker processes never
receive a file handle.

Enabling the timeline selects the existing ordered event-analysis path. Aggregate
text counts retain their meaning, but equal-count reaction ordering may differ
from the default unordered count-only path. With the timeline disabled, that
existing scheduling path is unchanged.
