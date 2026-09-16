# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compact, read-only access to ReacNetGenerator timeline schema 1.0.

Iterators read numeric columns in blocks and retain one variable-size definition
at a time. Closing an iterator releases its file handle. IDs are file-local;
consumers should join on IDs rather than infer scientific equality from them.
"""

import json
from contextlib import contextmanager
from dataclasses import dataclass
from operator import index

import h5py

SCHEMA_VERSION = "1.0"


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


@contextmanager
def _open(filename):
    with h5py.File(filename, "r") as file:
        if file.attrs.get("format") != "reacnetgenerator-timeline":
            raise ValueError("Not a ReacNetGenerator timeline")
        if file.attrs.get("schema_version") != SCHEMA_VERSION:
            raise ValueError("Unsupported timeline schema version")
        if file.attrs.get("status") != "complete":
            raise ValueError("Timeline is incomplete")
        yield file


def read_metadata(filename):
    """Read configuration and format attributes without loading result tables.

    This performs header checks only, not full semantic validation. Dataset IDs
    and offsets are public schema fields; a separate validator is planned.
    """
    with _open(filename) as file:
        result = dict(file.attrs)
        for key in ("configuration", "capabilities"):
            result[key] = json.loads(result[key])
        return result


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
    """Return a schema 1.0 signed 64-bit integer dataset."""
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
        names = file["species/name"].asstr()
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
