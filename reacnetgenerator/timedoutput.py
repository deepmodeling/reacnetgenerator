# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compact, read-only access to ReacNetGenerator timeline schemas 1.0 and 1.1.

Iterators read numeric columns in blocks and retain one variable-size definition
at a time. Closing an iterator releases its file handle. IDs are file-local;
consumers should join on IDs rather than infer scientific equality from them.
"""

import json
from contextlib import contextmanager
from dataclasses import dataclass
from operator import index

import h5py

from . import _timedoutputcontract as _contract

SCHEMA_VERSION = _contract.SCHEMA_VERSION
SUPPORTED_SCHEMA_VERSIONS = _contract.SUPPORTED_SCHEMA_VERSIONS
TimedOutputValidationError = _contract.TimedOutputValidationError
ValidationSummary = _contract.ValidationSummary


@dataclass(frozen=True)
class Frame:
    """An analyzed frame and its zero-based source-file frame."""

    frame: int
    source_id: int
    source_frame: int
    timestep: int


@dataclass(frozen=True)
class MoleculeRange:
    """One closed interval of a molecule's effective analysis signal."""

    molecule_id: int
    start_frame: int
    end_frame: int


@dataclass(frozen=True)
class ReactionEvent:
    """An aggregated reaction type at transition frame -> frame + 1."""

    transition: int
    reaction_type_id: int
    count: int


@dataclass(frozen=True)
class Molecule:
    """One definition; atom indices use RNG's zero-based canonical ordering."""

    molecule_id: int
    species_id: int
    atom_index: tuple[int, ...]
    bonds: tuple[tuple[int, int, int], ...]


@dataclass(frozen=True)
class TransitionParticipant:
    """One concrete molecule instance on one side of an inferred reaction."""

    side: str
    molecule_id: int
    species: str
    atom_index: tuple[int, ...]
    bonds: tuple[tuple[int, int, int], ...]


@dataclass(frozen=True)
class BondChange:
    """One inferred connectivity change for a canonical atom pair."""

    atom1: int
    atom2: int
    before_order: int
    after_order: int

    @property
    def kind(self):
        """Classify the change while retaining the before/after bond orders."""
        if self.before_order == 0:
            return "formed"
        if self.after_order == 0:
            return "broken"
        return "order_changed"


@dataclass(frozen=True)
class TransitionEvidence:
    """Auditable evidence for one inferred connected reaction instance."""

    transition: int
    reaction_type_id: int
    reactant: str
    product: str
    participants: tuple[TransitionParticipant, ...]
    bond_changes: tuple[BondChange, ...]


@contextmanager
def _open(filename):
    with h5py.File(filename, "r") as file:
        if file.attrs.get("format") != "reacnetgenerator-timeline":
            raise ValueError("Not a ReacNetGenerator timeline")
        if file.attrs.get("schema_version") not in SUPPORTED_SCHEMA_VERSIONS:
            raise ValueError("Unsupported timeline schema version")
        if file.attrs.get("status") != "complete":
            raise ValueError("Timeline is incomplete")
        yield file


def read_metadata(filename):
    """Read configuration and format attributes without loading result tables.

    This performs header checks only, not full semantic validation. Call
    :func:`validate_timed_output` when the complete artifact contract matters.
    """
    with _open(filename) as file:
        result = dict(file.attrs)
        for key in ("configuration", "capabilities"):
            result[key] = json.loads(result[key])
        return result


def validate_timed_output(filename, *, block_rows=8192):
    """Validate one complete timeline and return its verified table sizes."""
    from ._timedoutputvalidate import validate_timed_output as _validate

    return _validate(filename, block_rows=block_rows)


def semantic_manifest(filename, *, include_provenance=False, block_rows=8192):
    """Return a deterministic manifest for a validated timeline's meaning."""
    from ._timedoutputmanifest import semantic_manifest as _manifest

    return _manifest(
        filename,
        include_provenance=include_provenance,
        block_rows=block_rows,
    )


def compare_semantic_manifests(left, right):
    """Return JSON-pointer paths whose values differ between two manifests."""
    from ._timedoutputmanifest import compare_semantic_manifests as _compare

    return _compare(left, right)


def read_schema_descriptor():
    """Read the installed machine-readable descriptor for the current schema."""
    from ._timedoutputschema import read_schema_descriptor as _read

    return _read()


def _dataset(container, name):
    """Return one one-dimensional dataset or reject a malformed table header."""
    try:
        dataset = container[name]
    except KeyError as exc:
        raise ValueError(f"Missing column {name}") from exc
    if not isinstance(dataset, h5py.Dataset) or dataset.ndim != 1:
        raise ValueError(f"Invalid column shape for {name}")
    return dataset


def _numeric_dataset(container, name):
    """Return a timeline signed 64-bit integer dataset."""
    dataset = _dataset(container, name)
    if dataset.dtype.kind != "i" or dataset.dtype.itemsize != 8:
        raise ValueError(f"Invalid numeric column {name}")
    return dataset


def _rows(file, group, fields, block_rows):
    size = index(block_rows)
    if size <= 0:
        raise ValueError("block_rows must be a positive integer")
    datasets = [_numeric_dataset(file, f"{group}/{field}") for field in fields]
    length = len(datasets[0])
    if any(len(ds) != length for ds in datasets):
        raise ValueError(f"Misaligned columns in {group}")
    for start in range(0, length, size):
        blocks = [ds[start : start + size] for ds in datasets]
        for row in zip(*blocks, strict=True):
            yield tuple(int(x) for x in row)


def iter_frames(filename, *, block_rows=8192):
    """Yield frame mappings in analysis order, reading at most a block at once."""
    with _open(filename) as file:
        for frame, row in enumerate(
            _rows(file, "frames", ("source_id", "source_frame", "timestep"), block_rows)
        ):
            yield Frame(frame, *row)


def iter_molecule_ranges(filename, *, block_rows=8192):
    """Yield compact closed ranges in molecule order, without frame expansion."""
    with _open(filename) as file:
        for row in _rows(
            file,
            "molecule_ranges",
            ("molecule_id", "start_frame", "end_frame"),
            block_rows,
        ):
            yield MoleculeRange(*row)


def iter_reaction_events(filename, *, block_rows=8192):
    """Yield aggregate events in transition order, without repeating counts."""
    with _open(filename) as file:
        for row in _rows(
            file,
            "reaction_events",
            ("transition", "reaction_type_id", "count"),
            block_rows,
        ):
            yield ReactionEvent(*row)


def _transition_participant(
    molecule_id,
    side,
    *,
    species_ids,
    atom_offsets,
    atom_index,
    bond_offsets,
    bond_columns,
    names,
):
    """Read one referenced molecule definition for transition evidence."""
    row = molecule_id - 1
    if row < 0 or row >= len(species_ids):
        raise ValueError("Transition evidence references an invalid molecule_id")
    atom_start, atom_stop = (int(value) for value in atom_offsets[row : row + 2])
    bond_start, bond_stop = (int(value) for value in bond_offsets[row : row + 2])
    species_id = int(species_ids[row])
    if species_id < 0 or species_id >= len(names):
        raise ValueError("Transition evidence references an invalid species_id")
    return TransitionParticipant(
        side=side,
        molecule_id=molecule_id,
        species=names[species_id],
        atom_index=tuple(int(value) for value in atom_index[atom_start:atom_stop]),
        bonds=tuple(
            zip(
                *(
                    tuple(int(value) for value in column[bond_start:bond_stop])
                    for column in bond_columns
                ),
                strict=True,
            )
        ),
    )


def iter_transition_evidence(filename, *, block_rows=8192):
    """Yield instance participants and bond changes in transition order."""
    with _open(filename) as file:
        if file.attrs.get("schema_version") != SCHEMA_VERSION:
            raise ValueError("Timeline does not contain transition evidence")
        try:
            evidence = file["transition_evidence"]
        except KeyError as exc:
            raise ValueError("Timeline does not contain transition evidence") from exc
        if not isinstance(evidence, h5py.Group):
            raise ValueError("Invalid transition_evidence group")
        transitions = _numeric_dataset(evidence, "transition")
        reaction_ids = _numeric_dataset(evidence, "reaction_type_id")
        if len(transitions) != len(reaction_ids):
            raise ValueError("Misaligned columns in transition_evidence")
        participant_offsets = _numeric_dataset(evidence, "participant_offsets")
        participant_ids = _numeric_dataset(evidence, "participant_molecule_id")
        participant_sides = _numeric_dataset(evidence, "participant_side")
        change_offsets = _numeric_dataset(evidence, "bond_change_offsets")
        atom1 = _numeric_dataset(evidence, "bond_atom_index_1")
        atom2 = _numeric_dataset(evidence, "bond_atom_index_2")
        before = _numeric_dataset(evidence, "before_order")
        after = _numeric_dataset(evidence, "after_order")
        row_count = len(transitions)
        if len(participant_offsets) != row_count + 1:
            raise ValueError("Invalid transition evidence participant offsets")
        if len(change_offsets) != row_count + 1:
            raise ValueError("Invalid transition evidence bond-change offsets")
        if len(participant_ids) != len(participant_sides):
            raise ValueError("Misaligned transition evidence participants")
        if len({len(atom1), len(atom2), len(before), len(after)}) != 1:
            raise ValueError("Misaligned transition evidence bond changes")
        reaction_types = file["reaction_types"]
        reactants = _dataset(reaction_types, "reactant").asstr()
        products = _dataset(reaction_types, "product").asstr()
        molecules = file["molecules"]
        if not isinstance(molecules, h5py.Group):
            raise ValueError("Invalid molecules group")
        species_ids = _numeric_dataset(molecules, "species_id")
        atom_offsets = _numeric_dataset(molecules, "atom_offsets")
        atom_index = _numeric_dataset(molecules, "atom_index")
        bond_offsets = _numeric_dataset(molecules, "bond_offsets")
        bond_columns = tuple(
            _numeric_dataset(molecules, name)
            for name in ("bond_atom_index_1", "bond_atom_index_2", "bond_order")
        )
        names = _dataset(file, "species/name").asstr()
        for row, (transition, reaction_type_id) in enumerate(
            _rows(
                file,
                "transition_evidence",
                ("transition", "reaction_type_id"),
                block_rows,
            )
        ):
            if reaction_type_id < 0 or reaction_type_id >= len(reactants):
                raise ValueError(
                    "Transition evidence references an invalid reaction type"
                )
            participant_start, participant_stop = (
                int(value) for value in participant_offsets[row : row + 2]
            )
            change_start, change_stop = (
                int(value) for value in change_offsets[row : row + 2]
            )
            if not 0 <= participant_start <= participant_stop <= len(participant_ids):
                raise ValueError("Invalid transition evidence participant offsets")
            if not 0 <= change_start <= change_stop <= len(atom1):
                raise ValueError("Invalid transition evidence bond-change offsets")
            participants = []
            for molecule_id, side in zip(
                participant_ids[participant_start:participant_stop],
                participant_sides[participant_start:participant_stop],
                strict=True,
            ):
                side = int(side)
                if side not in (0, 1):
                    raise ValueError("Invalid transition evidence participant side")
                participants.append(
                    _transition_participant(
                        int(molecule_id),
                        "reactant" if side == 0 else "product",
                        species_ids=species_ids,
                        atom_offsets=atom_offsets,
                        atom_index=atom_index,
                        bond_offsets=bond_offsets,
                        bond_columns=bond_columns,
                        names=names,
                    )
                )
            bond_changes = tuple(
                BondChange(*(int(value) for value in values))
                for values in zip(
                    atom1[change_start:change_stop],
                    atom2[change_start:change_stop],
                    before[change_start:change_stop],
                    after[change_start:change_stop],
                    strict=True,
                )
            )
            yield TransitionEvidence(
                transition=transition,
                reaction_type_id=reaction_type_id,
                reactant=reactants[reaction_type_id],
                product=products[reaction_type_id],
                participants=tuple(participants),
                bond_changes=bond_changes,
            )


def iter_molecules(filename, *, block_rows=8192):
    """Yield definitions, retaining a block of IDs and one molecule's graph.

    A single molecule's atom and bond payload can exceed ``block_rows``. The
    component-size guard, not this reader's block size, bounds that graph.
    """
    with _open(filename) as file:
        molecules = file["molecules"]
        if not isinstance(molecules, h5py.Group):
            raise ValueError("Invalid molecules group")
        species_ids = _numeric_dataset(molecules, "species_id")
        atom_offsets = _numeric_dataset(molecules, "atom_offsets")
        atom_indices = _numeric_dataset(molecules, "atom_index")
        bond_offsets = _numeric_dataset(molecules, "bond_offsets")
        bond_columns = tuple(
            _numeric_dataset(molecules, name)
            for name in ("bond_atom_index_1", "bond_atom_index_2", "bond_order")
        )
        expected_offsets = len(species_ids) + 1
        if (
            len(atom_offsets) != expected_offsets
            or len(bond_offsets) != expected_offsets
        ):
            raise ValueError("Invalid molecule offset length")
        if len({len(dataset) for dataset in bond_columns}) != 1:
            raise ValueError("Misaligned molecule bond columns")
        for molecule_id, (species_id,) in enumerate(
            _rows(file, "molecules", ("species_id",), block_rows), 1
        ):
            atom_start, atom_stop = (
                int(x) for x in atom_offsets[molecule_id - 1 : molecule_id + 1]
            )
            bond_start, bond_stop = (
                int(x) for x in bond_offsets[molecule_id - 1 : molecule_id + 1]
            )
            atoms = tuple(int(x) for x in atom_indices[atom_start:atom_stop])
            bonds = tuple(
                zip(
                    *(
                        tuple(int(x) for x in dataset[bond_start:bond_stop])
                        for dataset in bond_columns
                    ),
                    strict=True,
                )
            )
            yield Molecule(molecule_id, species_id, atoms, bonds)


def iter_species(filename):
    """Yield ``(species_id, name)`` pairs, reading one UTF-8 name at a time."""
    with _open(filename) as file:
        names = _dataset(file, "species/name").asstr()
        for species_id in range(len(names)):
            yield species_id, names[species_id]


def iter_reaction_types(filename):
    """Yield ``(type_id, reactant, product, total_count)`` dictionary entries."""
    with _open(filename) as file:
        types = file["reaction_types"]
        if not isinstance(types, h5py.Group):
            raise ValueError("Invalid reaction_types group")
        reactants = _dataset(types, "reactant")
        products = _dataset(types, "product")
        totals = _numeric_dataset(types, "total_count")
        if len(reactants) != len(totals) or len(products) != len(totals):
            raise ValueError("Misaligned columns in reaction_types")
        reactant_text = reactants.asstr()
        product_text = products.asstr()
        for type_id in range(len(totals)):
            yield (
                type_id,
                reactant_text[type_id],
                product_text[type_id],
                int(totals[type_id]),
            )
