# Optional HDF5 timeline (schema 1.1)

Use `--timed-output FILE` to write one compact timeline when the PATH stage runs:

```sh
reacnetgenerator -i trajectory.bond --type bond -a H He --nohmm --timed-output timeline.h5
```

In Python, pass `timed_output="timeline.h5"` to `ReacNetGenerator`, or through
`reacnetgenerator.run`. The returned artifact map includes `timeline`. Relative
paths resolve against the working directory; this explicit path is not relocated
by `output_dir`. Its parent directory must exist. The default remains disabled.
A species-only `run_items` request cannot produce a timeline and raises an error
if combined with `timed_output`.

The timeline contains molecule definitions and effective presence intervals,
aggregate reaction events, instance-level participants and inferred bond changes,
and frame/source/configuration metadata. It does not contain coordinates. It does
not change the existing text outputs or enable the legacy CSV switches. Those
switches can still be selected separately.

## Why validation is part of the format

A saved analysis result is useful beyond the process that wrote it only when a
user or another program can establish what it contains and whether its tables
still agree. Opening an HDF5 file proves neither. The timeline validator makes
schemas 1.0 and 1.1 a verifiable software boundary:

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
  used to construct the atom-frame matrix: after HMM when enabled, observed
  presence otherwise. Overlapping molecule signals can exist and are handled by
  the existing reaction conflict rules. Ranges do not assert unique atom ownership.
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
| `schema_version`        | UTF-8 string `1.1`                                                                                 |
| `status`                | `incomplete` during construction, `complete` on successful close                                   |
| `capabilities`          | JSON array `["molecule_ranges", "reaction_events", "transition_evidence"]`                         |
| `rng_version`           | Installed RNG version                                                                              |
| `created_utc`           | ISO 8601 UTC timestamp                                                                             |
| `configuration`         | JSON from `parameter_provenance()`: normalized constructor parameters and explicit parameter names |
| `atom_index_convention` | `zero-based RNG canonical atom order`                                                              |
| `molecule_range_basis`  | `HMM signal` or `observed signal`                                                                  |

Configuration includes HMM settings, sampling, cell/PBC options, detection
backend/cutoffs, species identification and selection settings. A null configured
cell means that the parser derives it from input; actual frame cell matrices are
not copied into this format. Input paths, sizes and modification times provide
provenance, not a content-integrity guarantee. No input hash is computed by default.

All datasets are one-dimensional. Numeric columns are little-endian signed 64-bit
integers; text columns are variable-length UTF-8 strings. Columns in each ordinary
table have equal lengths. Dataset layout/chunking is not part of semantic identity.

| Group                 | Columns                                                                                        | Row meaning                                                                     |
| --------------------- | ---------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| `sources`             | `path` (text), `size_bytes`, `mtime_ns`                                                        | One input occurrence, in supplied order                                         |
| `frames`              | `source_id`, `source_frame`, `timestep`                                                        | One analyzed frame, in order                                                    |
| `atoms`               | `type`; `type_name` (text)                                                                     | `type` has one entry per atom; `type_name` is the separate atom-type dictionary |
| `species`             | `name` (text)                                                                                  | Unique species name                                                             |
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

Schema 1.1 is a compatible extension of schema 1.0 and advertises the
`transition_evidence` capability. Readers and the validator continue to accept
aggregate-only 1.0 files. Breaking changes require a new major version; unknown
versions fail rather than being guessed.

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
