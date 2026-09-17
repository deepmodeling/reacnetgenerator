# SPDX-License-Identifier: LGPL-3.0-or-later
"""Single-owner, bounded-batch writer for the opt-in timeline format."""

import itertools
import json
import os
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

from ._timedoutputrange import molecule_present
from ._version import __version__
from .utils import bytestolist, get_timestep_value, read_compressed_block

SCHEMA_VERSION = "1.1"
_BLOCK_ROWS = 8192
_BLOCK_BYTES = 1024 * 1024


def _signal_ranges(values, block_rows=_BLOCK_ROWS):
    """Find closed presence intervals using bounded scan temporaries."""
    if block_rows <= 0:
        raise ValueError("Scan block size must be positive")
    previous = False
    begin = None
    for start in range(0, len(values), block_rows):
        block = np.asarray(values[start : start + block_rows], dtype=np.bool_)
        # HMM/origin records use (frames, 1); callers may also supply 1-D
        # signals. Flatten only a single value per frame, within this block.
        if block.size != len(block):
            raise ValueError("Expected one signal value per frame")
        block = block.reshape(-1)
        changes = np.empty(len(block), dtype=np.bool_)
        changes[0] = block[0] != previous
        changes[1:] = block[1:] != block[:-1]
        for local in np.flatnonzero(changes):
            frame = start + int(local)
            if block[local]:
                begin = frame
            else:
                assert begin is not None
                yield begin, frame - 1
                begin = None
        previous = bool(block[-1])
    if begin is not None:
        yield begin, len(values) - 1


class _Column:
    """Buffer one column, bounded by rows and encoded bytes."""

    def __init__(self, file, name, dtype):
        self.dataset = file.create_dataset(
            name, (0,), maxshape=(None,), dtype=dtype, chunks=(_BLOCK_ROWS,)
        )
        self.rows = []
        self.bytes = 0

    def append(self, value):
        """Append a scalar; a single large string may exceed the byte budget."""
        size = len(value.encode("utf-8")) if isinstance(value, str) else 8
        if self.rows and self.bytes + size > _BLOCK_BYTES:
            self.flush()
        self.rows.append(value)
        self.bytes += size
        if len(self.rows) >= _BLOCK_ROWS or self.bytes >= _BLOCK_BYTES:
            self.flush()

    def flush(self):
        """Write a batch without retaining historical rows."""
        if self.rows:
            start = len(self.dataset)
            self.dataset.resize((start + len(self.rows),))
            self.dataset[start:] = self.rows
            self.rows.clear()
            self.bytes = 0


class _TimedOutputWriter:
    """Publish only after PATH succeeds, preserving any previous complete file.

    State is parent-owned. Workers never receive an HDF5 handle. Failed builds
    retain an explicitly incomplete sibling for diagnosis; the destination is
    replaced only after every dataset and completion marker has been closed.
    """

    def __init__(self, rng):
        self.rng = rng
        self.filename = Path(rng.timed_output).expanduser()
        self.columns = {}
        self.species = {}
        self.reactions = {}
        self.totals = Counter()
        self.atom_offset = self.bond_offset = 0
        self.participant_offset = self.bond_change_offset = 0
        self.evidence_count = 0
        self._molecule_columns_flushed = False
        self.file = None
        self.partial = None

    def __enter__(self):
        """Open an exclusive sibling and initialize the public schema."""
        fd, self.partial = tempfile.mkstemp(
            prefix=f".{self.filename.name}.",
            suffix=".incomplete",
            dir=self.filename.parent,
        )
        os.close(fd)
        try:
            self.file = h5py.File(self.partial, "w", rdcc_nbytes=_BLOCK_BYTES)
            self.file.attrs.update(
                format="reacnetgenerator-timeline",
                schema_version=SCHEMA_VERSION,
                status="incomplete",
                capabilities=json.dumps(
                    [
                        "molecule_ranges",
                        "reaction_events",
                        "transition_evidence",
                    ]
                ),
                rng_version=__version__,
                created_utc=datetime.now(timezone.utc).isoformat(),
                configuration=json.dumps(
                    self.rng.parameter_provenance(), allow_nan=False
                ),
                atom_index_convention="zero-based RNG canonical atom order",
                molecule_range_basis="HMM signal"
                if self.rng.runHMM
                else "observed signal",
            )
            for group, fields in {
                "frames": ("source_id", "source_frame", "timestep"),
                "atoms": ("type",),
                "molecules": (
                    "species_id",
                    "atom_offsets",
                    "atom_index",
                    "bond_offsets",
                    "bond_atom_index_1",
                    "bond_atom_index_2",
                    "bond_order",
                ),
                "molecule_ranges": ("molecule_id", "start_frame", "end_frame"),
                "reaction_events": ("transition", "reaction_type_id", "count"),
                "transition_evidence": (
                    "transition",
                    "reaction_type_id",
                    "participant_offsets",
                    "participant_molecule_id",
                    "participant_side",
                    "bond_change_offsets",
                    "bond_atom_index_1",
                    "bond_atom_index_2",
                    "before_order",
                    "after_order",
                ),
                "reaction_types": ("total_count",),
                "sources": ("size_bytes", "mtime_ns"),
            }.items():
                for field in fields:
                    key = f"{group}/{field}"
                    self.columns[key] = _Column(self.file, key, "<i8")
            for key in (
                "species/name",
                "reaction_types/reactant",
                "reaction_types/product",
                "sources/path",
                "atoms/type_name",
            ):
                self.columns[key] = _Column(self.file, key, h5py.string_dtype("utf-8"))
            self._append("molecules", atom_offsets=0, bond_offsets=0)
            self._append(
                "transition_evidence",
                participant_offsets=0,
                bond_change_offsets=0,
            )
            self._write_frames()
            return self
        except BaseException:
            if self.file is not None:
                self.file.close()
            raise

    def _append(self, group, **values):
        for name, value in values.items():
            self.columns[f"{group}/{name}"].append(value)

    def _write_frames(self):
        counts = self.rng.source_frame_counts
        if counts is None or len(counts) != len(self.rng.inputfilename):
            raise ValueError("Timed output requires source frame mapping from DETECT")
        for path in self.rng.inputfilename:
            stat = os.stat(path)
            self._append(
                "sources",
                path=os.fspath(path),
                size_bytes=stat.st_size,
                mtime_ns=stat.st_mtime_ns,
            )
        for name in self.rng.atomname:
            self._append("atoms", type_name=str(name))
        for atomtype in self.rng.atomtype:
            self._append("atoms", type=int(atomtype))
        source = offset = 0
        for frame in range(self.rng.step):
            original = frame * self.rng.stepinterval
            while source < len(counts) and original >= offset + counts[source]:
                offset += counts[source]
                source += 1
            if source == len(counts) or frame not in self.rng.timestep:
                raise ValueError("Incomplete source frame or timestep mapping")
            self._append(
                "frames",
                source_id=source,
                source_frame=original - offset,
                timestep=int(get_timestep_value(self.rng.timestep[frame])),
            )
        if self.rng.step != len(range(0, sum(counts), self.rng.stepinterval)):
            raise ValueError("Source frame count does not match analyzed frames")

    def write_molecules(self, names):
        """Stream definitions and effective signal ranges in molecule-ID order.

        Use the very same signal records as the atom-frame matrix, including
        HMM-created/removed observations. Legacy CSV instead describes observed
        frames, so HMM-on CSV is deliberately not used as a range oracle.
        """
        count = 0
        with (
            open(self.rng.moleculetemp2filename, "rb") as definitions,
            open(
                self.rng.hmmfilename if self.rng.runHMM else self.rng.originfilename,
                "rb",
            ) as signals,
        ):
            records = itertools.zip_longest(*[read_compressed_block(definitions)] * 4)
            for count, (record, signal) in enumerate(
                zip(records, read_compressed_block(signals), strict=True), 1
            ):
                if count > self.rng.hmmit or any(x is None for x in record):
                    raise ValueError("Invalid molecule record alignment")
                name = str(names[count - 1])
                if name not in self.species:
                    self.species[name] = len(self.species)
                    self._append("species", name=name)
                self._append("molecules", species_id=self.species[name])
                for atom in bytestolist(record[0]):
                    self._append("molecules", atom_index=int(atom))
                    self.atom_offset += 1
                for pair, level in zip(
                    bytestolist(record[1]), bytestolist(record[2]), strict=True
                ):
                    self._append(
                        "molecules",
                        bond_atom_index_1=int(pair[0]),
                        bond_atom_index_2=int(pair[1]),
                        bond_order=int(level),
                    )
                    self.bond_offset += 1
                self._append(
                    "molecules",
                    atom_offsets=self.atom_offset,
                    bond_offsets=self.bond_offset,
                )
                values = bytestolist(signal)
                if len(values) != self.rng.step:
                    raise ValueError(
                        "Molecule signal length differs from analyzed frame count"
                    )
                for start, stop in _signal_ranges(values):
                    self._append(
                        "molecule_ranges",
                        molecule_id=count,
                        start_frame=start,
                        end_frame=stop,
                    )
        if count != self.rng.hmmit:
            raise ValueError("Molecule count differs from declared count")

    def write_events(self, events):
        """Store aggregate counts and one auditable row per inferred instance."""
        events = tuple(events)
        counts = Counter(
            (e["Timestep_Index"], e["Reactant"], e["Product"]) for e in events
        )
        for (transition, left, right), count in counts.items():
            pair = (left, right)
            if pair not in self.reactions:
                self.reactions[pair] = len(self.reactions)
                self._append("reaction_types", reactant=left, product=right)
            reaction_id = self.reactions[pair]
            self._append(
                "reaction_events",
                transition=transition,
                reaction_type_id=reaction_id,
                count=count,
            )
            self.totals[reaction_id] += count
        for event in events:
            pair = (event["Reactant"], event["Product"])
            self._write_transition_evidence(event, self.reactions[pair])

    def _flush_molecule_tables(self):
        """Make molecule definitions and ranges readable before event checks."""
        if self._molecule_columns_flushed:
            return
        for name, column in self.columns.items():
            if name.startswith(("molecules/", "molecule_ranges/")):
                column.flush()
        self._molecule_columns_flushed = True

    def _participant_atoms(self, molecule_ids):
        """Return the disjoint atom union for one reaction side."""
        self._flush_molecule_tables()
        species_ids = self.columns["molecules/species_id"].dataset
        offsets = self.columns["molecules/atom_offsets"].dataset
        atom_index = self.columns["molecules/atom_index"].dataset
        atoms = set()
        for molecule_id in molecule_ids:
            row = int(molecule_id) - 1
            if row < 0 or row >= len(species_ids):
                raise ValueError(
                    "Transition evidence references an invalid molecule_id"
                )
            start, stop = (int(value) for value in offsets[row : row + 2])
            current = {int(value) for value in atom_index[start:stop]}
            if atoms & current:
                raise ValueError(
                    "Transition evidence participants overlap in atom_index"
                )
            atoms.update(current)
        return atoms

    def _participant_bonds(self, molecule_ids):
        """Return the union of stored participant bonds for one reaction side."""
        self._flush_molecule_tables()
        offsets = self.columns["molecules/bond_offsets"].dataset
        atom1 = self.columns["molecules/bond_atom_index_1"].dataset
        atom2 = self.columns["molecules/bond_atom_index_2"].dataset
        order = self.columns["molecules/bond_order"].dataset
        bonds = {}
        for molecule_id in molecule_ids:
            row = int(molecule_id) - 1
            start, stop = (int(value) for value in offsets[row : row + 2])
            for left, right, level in zip(
                atom1[start:stop], atom2[start:stop], order[start:stop], strict=True
            ):
                pair = tuple(sorted((int(left), int(right))))
                level = int(level)
                if pair in bonds and bonds[pair] != level:
                    raise ValueError(
                        "Transition evidence participants contain conflicting bonds"
                    )
                bonds[pair] = level
        return bonds

    def _participant_species(self, molecule_ids):
        """Count stored species names for one reaction side."""
        self._flush_molecule_tables()
        species_ids = self.columns["molecules/species_id"].dataset
        names = {species_id: name for name, species_id in self.species.items()}
        counts = Counter()
        for molecule_id in molecule_ids:
            species_id = int(species_ids[molecule_id - 1])
            try:
                counts[names[species_id]] += 1
            except KeyError as exc:
                raise ValueError(
                    "Transition evidence references an invalid species_id"
                ) from exc
        return counts

    def _write_transition_evidence(self, event, reaction_type_id):
        """Append one reaction instance with participants and inferred bond changes."""
        reactants = tuple(int(value) for value in event["ReactantMoleculeIDs"])
        products = tuple(int(value) for value in event["ProductMoleculeIDs"])
        if (
            not reactants
            or not products
            or tuple(sorted(set(reactants))) != reactants
            or tuple(sorted(set(products))) != products
        ):
            raise ValueError(
                "Transition evidence participants must be nonempty, unique, and ordered"
            )

        self._flush_molecule_tables()
        range_ids = self.columns["molecule_ranges/molecule_id"].dataset
        range_starts = self.columns["molecule_ranges/start_frame"].dataset
        range_ends = self.columns["molecule_ranges/end_frame"].dataset
        transition = int(event["Timestep_Index"])
        for side, molecule_ids in enumerate((reactants, products)):
            frame = transition + side
            if any(
                not molecule_present(
                    molecule_id,
                    frame,
                    range_ids,
                    range_starts,
                    range_ends,
                )
                for molecule_id in molecule_ids
            ):
                raise ValueError(
                    "Transition evidence participant is absent from its reaction frame"
                )

        reactant_atoms = self._participant_atoms(reactants)
        product_atoms = self._participant_atoms(products)
        if reactant_atoms != product_atoms:
            raise ValueError("Transition evidence participants do not conserve atoms")
        reactant_species = self._participant_species(reactants)
        product_species = self._participant_species(products)
        pair = tuple(
            "+".join(sorted(side.elements()))
            for side in (
                reactant_species - product_species,
                product_species - reactant_species,
            )
        )
        if not all(pair) or pair != (event["Reactant"], event["Product"]):
            raise ValueError(
                "Transition evidence participant species disagree with reaction type"
            )

        self._append(
            "transition_evidence",
            transition=transition,
            reaction_type_id=int(reaction_type_id),
        )
        for side, molecule_ids in enumerate((reactants, products)):
            for molecule_id in molecule_ids:
                self._append(
                    "transition_evidence",
                    participant_molecule_id=molecule_id,
                    participant_side=side,
                )
                self.participant_offset += 1
        self._append("transition_evidence", participant_offsets=self.participant_offset)

        before = self._participant_bonds(reactants)
        after = self._participant_bonds(products)
        for atom_pair in sorted(before.keys() | after.keys()):
            before_order = before.get(atom_pair, 0)
            after_order = after.get(atom_pair, 0)
            if before_order == after_order:
                continue
            self._append(
                "transition_evidence",
                bond_atom_index_1=atom_pair[0],
                bond_atom_index_2=atom_pair[1],
                before_order=before_order,
                after_order=after_order,
            )
            self.bond_change_offset += 1
        self._append("transition_evidence", bond_change_offsets=self.bond_change_offset)
        self.evidence_count += 1

    def __exit__(self, exc_type, exc, traceback):
        """Close before atomic replacement; errors preserve the old destination."""
        assert self.file is not None and self.partial is not None
        try:
            if exc_type is None:
                if self.evidence_count != sum(self.totals.values()):
                    raise ValueError(
                        "Transition evidence count does not match reaction events"
                    )
                for reaction_id in range(len(self.reactions)):
                    self._append("reaction_types", total_count=self.totals[reaction_id])
                for column in self.columns.values():
                    column.flush()
                self.file.flush()
                self.file.attrs["status"] = "complete"
                self.file.flush()
        finally:
            self.file.close()
        if exc_type is None:
            os.replace(self.partial, self.filename)
