# Optional HDF5 timeline (schema 1.0)

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
aggregate reaction events, and frame/source/configuration metadata. It does not
contain coordinates or instance-level reaction participants/bond-change evidence.
It does not change the existing text outputs or enable the legacy CSV switches.
Those switches can still be selected separately.

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
| `schema_version`        | UTF-8 string `1.0`                                                                                 |
| `status`                | `incomplete` during construction, `complete` on successful close                                   |
| `capabilities`          | JSON array `["molecule_ranges", "reaction_events"]`                                                |
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

| Group             | Columns                                                                | Row meaning                                                                     |
| ----------------- | ---------------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| `sources`         | `path` (text), `size_bytes`, `mtime_ns`                                | One input occurrence, in supplied order                                         |
| `frames`          | `source_id`, `source_frame`, `timestep`                                | One analyzed frame, in order                                                    |
| `atoms`           | `type`; `type_name` (text)                                             | `type` has one entry per atom; `type_name` is the separate atom-type dictionary |
| `species`         | `name` (text)                                                          | Unique species name                                                             |
| `molecules`       | `species_id`                                                           | One molecule definition, row `molecule_id - 1`                                  |
| `molecules`       | `atom_offsets`, `atom_index`                                           | Offset-delimited atom lists                                                     |
| `molecules`       | `bond_offsets`, `bond_atom_index_1`, `bond_atom_index_2`, `bond_order` | Offset-delimited bonds with global canonical atom endpoints                     |
| `molecule_ranges` | `molecule_id`, `start_frame`, `end_frame`                              | One maximal effective interval, ordered by molecule then start frame            |
| `reaction_types`  | `reactant` (text), `product` (text), `total_count`                     | One unique pair of existing formatted reaction sides                            |
| `reaction_events` | `transition`, `reaction_type_id`, `count`                              | One positive aggregate count, ordered by transition                             |

Both offset columns start at zero and have `number_of_molecules + 1` entries.
For molecule ID `m`, its payload is `[offset[m-1]:offset[m]]`. The final offset
is the payload length. Bond endpoint/order columns have identical lengths. Empty
payloads and empty tables are valid. Reaction-type totals equal the sum of stored
event counts for that type. Reaction-side strings retain the existing output
notation; consumers should not assume splitting on `+` parses every possible SMILES.

A compatible future evidence extension will use a minor version and explicit
capability. Breaking changes require a new major version. The present reader
accepts exactly `1.0`; it fails on an unknown version rather than guessing.

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
)

metadata = read_metadata("timeline.h5")
for interval in iter_molecule_ranges("timeline.h5", block_rows=4096):
    print(interval.molecule_id, interval.start_frame, interval.end_frame)
for event in iter_reaction_events("timeline.h5", block_rows=4096):
    print(event.transition, event.reaction_type_id, event.count)
```

`iter_frames`, `iter_molecules`, `iter_molecule_ranges` and `iter_reaction_events`
return frozen dataclass records. `iter_species` yields `(species_id, name)`;
`iter_reaction_types` yields `(type_id, reactant, product, total_count)`.
Definitions include `atom_index` and `(atom1, atom2, order)` bond tuples.
Source paths and atom-type tables can be accessed directly with h5py under the
public contract above.

Numeric rows are read in blocks, ranges are not expanded into frame rows, and
counts are not expanded into individual events. Dictionary strings are read one
at a time. One molecule definition or one string can exceed the block budget.
Call `.close()` on an iterator when stopping early, or use `contextlib.closing`.
Header/type/alignment checks are not a full structural or semantic validator;
a dedicated validator and semantic manifest are a later batch.

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
scale with content. This feature does not promise constant total memory. HDF5 is
written only in the parent process; worker processes never receive a file handle.

Enabling the timeline selects the existing ordered event-analysis path. Aggregate
text counts retain their meaning, but equal-count reaction ordering may differ
from the default unordered count-only path. With the timeline disabled, that
existing scheduling path is unchanged.
