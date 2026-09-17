# SPDX-License-Identifier: LGPL-3.0-or-later
"""Bounded structural and cross-table validation for timeline schema 1.0."""

import json
import sys
from collections import Counter
from datetime import datetime, timezone
from operator import index

import h5py

from .timedoutput import (
    SCHEMA_VERSION,
    TimedOutputValidationError,
    ValidationSummary,
)

_CAPABILITIES = ["molecule_ranges", "reaction_events"]
_ATOM_CONVENTION = "zero-based RNG canonical atom order"
_RANGE_BASES = {"HMM signal", "observed signal"}


def _reject_json_constant(value):
    raise ValueError(f"Non-standard JSON constant {value}")


def _error(message):
    raise TimedOutputValidationError(message)


def _dataset(file, path):
    try:
        dataset = file[path]
    except KeyError:
        _error(f"Missing dataset {path}")
    if not isinstance(dataset, h5py.Dataset) or dataset.ndim != 1:
        _error(f"{path} must be a one-dimensional dataset")
    return dataset


def _numeric(file, path):
    dataset = _dataset(file, path)
    byteorder = dataset.dtype.byteorder
    little_endian = byteorder == "<" or (byteorder == "=" and sys.byteorder == "little")
    if dataset.dtype.kind != "i" or dataset.dtype.itemsize != 8 or not little_endian:
        _error(f"{path} must contain little-endian signed 64-bit integers")
    return dataset


def _text(file, path):
    dataset = _dataset(file, path)
    string = h5py.check_string_dtype(dataset.dtype)
    if string is None or string.encoding != "utf-8" or string.length is not None:
        _error(f"{path} must contain variable-length UTF-8 strings")
    return dataset


def _aligned(table, datasets):
    lengths = {len(dataset) for dataset in datasets}
    if len(lengths) != 1:
        paths = ", ".join(dataset.name.lstrip("/") for dataset in datasets)
        _error(f"{table} has misaligned columns: {paths}")


def _integer_blocks(dataset, block_rows):
    for start in range(0, len(dataset), block_rows):
        for raw in dataset[start : start + block_rows]:
            yield int(raw)


def _text_blocks(dataset, block_rows):
    values = dataset.asstr()
    for start in range(0, len(dataset), block_rows):
        yield from values[start : start + block_rows]


def _bounded(dataset, *, lower, upper, block_rows):
    for position, value in enumerate(_integer_blocks(dataset, block_rows)):
        if value < lower or value >= upper:
            _error(
                f"{dataset.name.lstrip('/')}[{position}]={value} is outside "
                f"[{lower}, {upper})"
            )


def _nonnegative(dataset, block_rows):
    for position, value in enumerate(_integer_blocks(dataset, block_rows)):
        if value < 0:
            _error(f"{dataset.name.lstrip('/')}[{position}] must be nonnegative")


def _unique_text(dataset, block_rows):
    seen = set()
    for position, value in enumerate(_text_blocks(dataset, block_rows)):
        if value in seen:
            _error(f"{dataset.name.lstrip('/')}[{position}] duplicates {value!r}")
        seen.add(value)


def _offsets(dataset, *, payload, table, block_rows):
    """Check one offset column without retaining its complete contents."""
    previous = None
    for position, value in enumerate(_integer_blocks(dataset, block_rows)):
        if position == 0 and value != 0:
            _error(f"{table} must start at zero")
        if previous is not None and value < previous:
            _error(f"{table} must be nondecreasing")
        previous = value
    final = 0 if previous is None else previous
    if final != len(payload):
        _error(
            f"{table} endpoint {final} does not match {payload.name.lstrip('/')} "
            f"length {len(payload)}"
        )


def _root_metadata(file):
    required = {
        "format": "reacnetgenerator-timeline",
        "schema_version": SCHEMA_VERSION,
        "status": "complete",
        "atom_index_convention": _ATOM_CONVENTION,
    }
    for name, expected in required.items():
        if name not in file.attrs:
            _error(f"Missing root attribute {name}")
        if file.attrs[name] != expected:
            _error(f"Invalid root attribute {name}: expected {expected!r}")
    range_basis = file.attrs.get("molecule_range_basis")
    if range_basis not in _RANGE_BASES:
        _error("Invalid root attribute molecule_range_basis")
    for name in ("rng_version", "created_utc", "configuration", "capabilities"):
        if not isinstance(file.attrs.get(name), str) or not file.attrs[name]:
            _error(f"Invalid root attribute {name}")
    try:
        created = datetime.fromisoformat(file.attrs["created_utc"])
    except ValueError:
        _error("Invalid root attribute created_utc")
    if created.utcoffset() != timezone.utc.utcoffset(created):
        _error("Root attribute created_utc must identify UTC")
    try:
        capabilities = json.loads(
            file.attrs["capabilities"], parse_constant=_reject_json_constant
        )
    except (TypeError, ValueError):
        _error("Root attribute capabilities must be valid JSON")
    if capabilities != _CAPABILITIES:
        _error(f"Root attribute capabilities must equal {_CAPABILITIES!r}")
    try:
        configuration = json.loads(
            file.attrs["configuration"], parse_constant=_reject_json_constant
        )
    except (TypeError, ValueError):
        _error("Root attribute configuration must be valid JSON")
    if not isinstance(configuration, dict):
        _error("Root attribute configuration must be a JSON object")
    parameters = configuration.get("parameters")
    if not isinstance(parameters, dict):
        _error("configuration/parameters must be a JSON object")
    explicit = configuration.get("explicit_parameters")
    if not isinstance(explicit, list) or not all(isinstance(x, str) for x in explicit):
        _error("configuration/explicit_parameters must be an array of strings")
    run_hmm = parameters.get("runHMM")
    expected_basis = (
        "HMM signal"
        if run_hmm is True
        else "observed signal"
        if run_hmm is False
        else None
    )
    if expected_basis != range_basis:
        _error(
            "configuration/parameters/runHMM must be boolean and agree with "
            "molecule_range_basis"
        )
    stepinterval = parameters.get("stepinterval")
    if (
        not isinstance(stepinterval, int)
        or isinstance(stepinterval, bool)
        or stepinterval <= 0
    ):
        _error("configuration/parameters/stepinterval must be a positive integer")
    return configuration


def _validate_frames(file, *, source_count, stepinterval, block_rows):
    source_id = _numeric(file, "frames/source_id")
    source_frame = _numeric(file, "frames/source_frame")
    timestep = _numeric(file, "frames/timestep")
    _aligned("frames", (source_id, source_frame, timestep))
    _bounded(source_id, lower=0, upper=source_count, block_rows=block_rows)
    _nonnegative(source_frame, block_rows)
    previous_source = previous_frame = None
    for position, (source, frame) in enumerate(
        zip(
            _integer_blocks(source_id, block_rows),
            _integer_blocks(source_frame, block_rows),
            strict=True,
        )
    ):
        if previous_source is None:
            if frame != 0:
                _error(
                    "frames/source_frame[0] must start at zero for the declared "
                    "stepinterval"
                )
        elif source < previous_source:
            _error(f"frames/source_id[{position}] breaks source order")
        elif source == previous_source:
            if frame != previous_frame + stepinterval:
                _error(
                    f"frames/source_frame[{position}] is inconsistent with "
                    "configuration/parameters/stepinterval"
                )
        elif frame >= stepinterval:
            _error(
                f"frames/source_frame[{position}] is inconsistent with "
                "configuration/parameters/stepinterval at a source boundary"
            )
        previous_source, previous_frame = source, frame
    return source_id


def _validate_atoms(file, *, block_rows):
    atom_type = _numeric(file, "atoms/type")
    type_name = _text(file, "atoms/type_name")
    _unique_text(type_name, block_rows)
    _bounded(atom_type, lower=0, upper=len(type_name), block_rows=block_rows)
    return atom_type, type_name


def _validate_molecules(file, *, atom_count, species_count, block_rows):
    species_id = _numeric(file, "molecules/species_id")
    atom_offsets = _numeric(file, "molecules/atom_offsets")
    atom_index = _numeric(file, "molecules/atom_index")
    bond_offsets = _numeric(file, "molecules/bond_offsets")
    bond_atom_1 = _numeric(file, "molecules/bond_atom_index_1")
    bond_atom_2 = _numeric(file, "molecules/bond_atom_index_2")
    bond_order = _numeric(file, "molecules/bond_order")
    if len(atom_offsets) != len(species_id) + 1:
        _error(
            "molecules/atom_offsets must have one more row than molecules/species_id"
        )
    if len(bond_offsets) != len(species_id) + 1:
        _error(
            "molecules/bond_offsets must have one more row than molecules/species_id"
        )
    _aligned("molecules bond payload", (bond_atom_1, bond_atom_2, bond_order))
    _offsets(
        atom_offsets,
        payload=atom_index,
        table="molecules/atom_offsets",
        block_rows=block_rows,
    )
    _offsets(
        bond_offsets,
        payload=bond_atom_1,
        table="molecules/bond_offsets",
        block_rows=block_rows,
    )
    _bounded(species_id, lower=0, upper=species_count, block_rows=block_rows)
    _bounded(atom_index, lower=0, upper=atom_count, block_rows=block_rows)
    _bounded(bond_atom_1, lower=0, upper=atom_count, block_rows=block_rows)
    _bounded(bond_atom_2, lower=0, upper=atom_count, block_rows=block_rows)

    for molecule in range(len(species_id)):
        atom_start, atom_stop = (int(x) for x in atom_offsets[molecule : molecule + 2])
        bond_start, bond_stop = (int(x) for x in bond_offsets[molecule : molecule + 2])
        atoms = [int(x) for x in atom_index[atom_start:atom_stop]]
        if len(atoms) != len(set(atoms)):
            _error(
                f"molecules/atom_index contains duplicates for molecule {molecule + 1}"
            )
        members = set(atoms)
        for column in (bond_atom_1, bond_atom_2):
            for endpoint in (int(x) for x in column[bond_start:bond_stop]):
                if endpoint not in members:
                    _error(
                        f"{column.name.lstrip('/')} endpoint {endpoint} is not in "
                        f"molecule {molecule + 1}"
                    )
    return species_id


def _validate_ranges(file, *, molecule_count, frame_count, block_rows):
    molecule_id = _numeric(file, "molecule_ranges/molecule_id")
    start_frame = _numeric(file, "molecule_ranges/start_frame")
    end_frame = _numeric(file, "molecule_ranges/end_frame")
    _aligned("molecule_ranges", (molecule_id, start_frame, end_frame))
    _bounded(molecule_id, lower=1, upper=molecule_count + 1, block_rows=block_rows)
    previous = None
    rows = zip(
        _integer_blocks(molecule_id, block_rows),
        _integer_blocks(start_frame, block_rows),
        _integer_blocks(end_frame, block_rows),
        strict=True,
    )
    for position, row in enumerate(rows):
        molecule, start, end = row
        if start < 0 or end < start or end >= frame_count:
            _error(
                "molecule_ranges/start_frame or molecule_ranges/end_frame is "
                f"invalid at row {position}"
            )
        if previous is not None:
            old_molecule, old_start, old_end = previous
            if (molecule, start) < (old_molecule, old_start):
                _error(f"molecule_ranges row {position} breaks declared ordering")
            if molecule == old_molecule and start <= old_end + 1:
                _error(f"molecule_ranges row {position} is overlapping or not maximal")
        previous = row
    return molecule_id


def _validate_reactions(file, *, frame_count, block_rows):
    reactant = _text(file, "reaction_types/reactant")
    product = _text(file, "reaction_types/product")
    total_count = _numeric(file, "reaction_types/total_count")
    _aligned("reaction_types", (reactant, product, total_count))
    pairs = set()
    for position, pair in enumerate(
        zip(
            _text_blocks(reactant, block_rows),
            _text_blocks(product, block_rows),
            strict=True,
        )
    ):
        if pair in pairs:
            _error(f"reaction_types row {position} duplicates a reactant/product pair")
        pairs.add(pair)
    for position, total in enumerate(_integer_blocks(total_count, block_rows)):
        if total <= 0:
            _error(f"reaction_types/total_count[{position}] must be positive")

    transition = _numeric(file, "reaction_events/transition")
    reaction_type_id = _numeric(file, "reaction_events/reaction_type_id")
    count = _numeric(file, "reaction_events/count")
    _aligned("reaction_events", (transition, reaction_type_id, count))
    _bounded(
        reaction_type_id,
        lower=0,
        upper=len(total_count),
        block_rows=block_rows,
    )
    sums = Counter()
    previous_transition = None
    seen_at_transition = set()
    rows = zip(
        _integer_blocks(transition, block_rows),
        _integer_blocks(reaction_type_id, block_rows),
        _integer_blocks(count, block_rows),
        strict=True,
    )
    for position, (frame, type_id, value) in enumerate(rows):
        if frame < 0 or frame >= frame_count - 1:
            _error(
                f"reaction_events/transition[{position}]={frame} is outside the timeline"
            )
        if value <= 0:
            _error(f"reaction_events/count[{position}] must be positive")
        if previous_transition is not None and frame < previous_transition:
            _error(f"reaction_events/transition[{position}] breaks declared ordering")
        if frame != previous_transition:
            seen_at_transition.clear()
        if type_id in seen_at_transition:
            _error(
                f"reaction_events/reaction_type_id[{position}] duplicates type {type_id} "
                f"at transition {frame}"
            )
        seen_at_transition.add(type_id)
        sums[type_id] += value
        previous_transition = frame
    for type_id, expected in enumerate(_integer_blocks(total_count, block_rows)):
        if sums[type_id] != expected:
            _error(
                f"reaction_types/total_count[{type_id}]={expected} does not match "
                f"reaction_events sum {sums[type_id]}"
            )
    return transition, total_count


def validate_timed_output(filename, *, block_rows=8192):
    """Implement the public validator while keeping HDF5 details private."""
    try:
        size = index(block_rows)
    except TypeError as exc:
        raise TimedOutputValidationError(
            "block_rows must be a positive integer"
        ) from exc
    if size <= 0:
        raise TimedOutputValidationError("block_rows must be a positive integer")

    try:
        with h5py.File(filename, "r") as file:
            configuration = _root_metadata(file)
            source_path = _text(file, "sources/path")
            source_size = _numeric(file, "sources/size_bytes")
            source_mtime = _numeric(file, "sources/mtime_ns")
            _aligned("sources", (source_path, source_size, source_mtime))
            _nonnegative(source_size, size)

            frames = _validate_frames(
                file,
                source_count=len(source_path),
                stepinterval=configuration["parameters"]["stepinterval"],
                block_rows=size,
            )
            atoms, atom_types = _validate_atoms(file, block_rows=size)
            species = _text(file, "species/name")
            _unique_text(species, size)
            molecules = _validate_molecules(
                file,
                atom_count=len(atoms),
                species_count=len(species),
                block_rows=size,
            )
            ranges = _validate_ranges(
                file,
                molecule_count=len(molecules),
                frame_count=len(frames),
                block_rows=size,
            )
            events, reaction_types = _validate_reactions(
                file, frame_count=len(frames), block_rows=size
            )

            return ValidationSummary(
                schema_version=SCHEMA_VERSION,
                sources=len(source_path),
                frames=len(frames),
                atoms=len(atoms),
                atom_types=len(atom_types),
                species=len(species),
                molecules=len(molecules),
                molecule_ranges=len(ranges),
                reaction_types=len(reaction_types),
                reaction_events=len(events),
            )
    except TimedOutputValidationError:
        raise
    except (OSError, TypeError, ValueError) as exc:
        raise TimedOutputValidationError(str(exc)) from exc
