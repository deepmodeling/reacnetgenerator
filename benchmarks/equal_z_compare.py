# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compare matched benchmark outputs independently of parallel molecule IDs.

Each timeline first passes the official validator. Exact species names, atom
membership, bond orders, effective ranges and transition evidence are retained.
This helper targets bounded benchmark windows; it does not normalize chemistry.
"""

import argparse
import hashlib
import json
import struct
from collections import Counter
from pathlib import Path

import h5py

from reacnetgenerator.timedoutput import (
    iter_frames,
    iter_molecule_ranges,
    iter_molecules,
    iter_reaction_events,
    iter_reaction_types,
    iter_species,
    iter_transition_evidence,
    validate_timed_output,
)


def encoded(value):
    """Encode one canonical JSON record without platform-dependent whitespace."""
    return json.dumps(value, ensure_ascii=False, separators=(",", ":")).encode()


def digest(records):
    """Hash an ordered sequence of unambiguously delimited records."""
    result = hashlib.sha256()
    for record in records:
        value = encoded(record)
        result.update(len(value).to_bytes(8, "little"))
        result.update(value)
    return result.hexdigest()


def molecule_key(species, atoms, bonds):
    """Preserve an exact labelled graph while sorting its storage order."""
    return digest(
        [
            [
                species,
                sorted(atoms),
                sorted((min(a, b), max(a, b), order) for a, b, order in bonds),
            ]
        ]
    )


def participants(text):
    """Sort reaction participants without splitting bracketed formal charges."""
    values, current, depth = [], [], 0
    for character in text:
        if character == "[":
            depth += 1
        elif character == "]":
            depth -= 1
        if character == "+" and depth == 0:
            values.append("".join(current))
            current = []
        else:
            current.append(character)
    values.append("".join(current))
    return tuple(sorted(values))


def fingerprint(directory):
    """Validate and fingerprint one pipeline-mode output directory."""
    timeline = directory / "timeline.h5"
    summary = validate_timed_output(timeline)
    if summary.schema_version != "1.1":
        raise ValueError("This benchmark comparison requires schema 1.1")
    names = dict(iter_species(timeline))
    keys = {
        molecule.molecule_id: molecule_key(
            names[molecule.species_id], molecule.atom_index, molecule.bonds
        )
        for molecule in iter_molecules(timeline)
    }
    intervals = {molecule_id: hashlib.sha256() for molecule_id in keys}
    for interval in iter_molecule_ranges(timeline):
        intervals[interval.molecule_id].update(
            struct.pack("<qq", interval.start_frame, interval.end_frame)
        )
    types = {
        type_id: (participants(left), participants(right))
        for type_id, left, right, _ in iter_reaction_types(timeline)
    }
    events = Counter()
    for event in iter_reaction_events(timeline):
        events[(event.transition, *types[event.reaction_type_id])] += event.count
    evidence = []
    for event in iter_transition_evidence(timeline):
        evidence.append(
            (
                event.transition,
                participants(event.reactant),
                participants(event.product),
                sorted(
                    (part.side, keys[part.molecule_id]) for part in event.participants
                ),
                sorted(
                    (
                        min(change.atom1, change.atom2),
                        max(change.atom1, change.atom2),
                        change.before_order,
                        change.after_order,
                    )
                    for change in event.bond_changes
                ),
            )
        )
    populations = []
    with (directory / "species").open() as file:
        for line in file:
            fields = line.split()
            populations.append(
                (
                    int(fields[1].rstrip(":")),
                    sorted(
                        (name, int(count))
                        for name, count in zip(fields[2::2], fields[3::2], strict=True)
                    ),
                )
            )
    with h5py.File(timeline, "r") as file:
        atom_types = [
            file["atoms/type"][:].tolist(),
            file["atoms/type_name"].asstr()[:].tolist(),
        ]
    return {
        "frames": summary.frames,
        "atom_types": digest(atom_types),
        "frame_mapping": digest(
            (f.source_id, f.source_frame, f.timestep) for f in iter_frames(timeline)
        ),
        "definitions": digest(sorted(keys.values())),
        "presence": digest(
            sorted((keys[mid], h.hexdigest()) for mid, h in intervals.items())
        ),
        "aggregate_events": digest(sorted(events.items())),
        "transition_evidence": digest(sorted(evidence)),
        "species_populations": digest(populations),
    }


def main():
    """Compare two output directories, failing if any exact evidence differs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    base, candidate = fingerprint(args.base), fingerprint(args.candidate)
    differences = [name for name in base if base[name] != candidate[name]]
    result = {
        "equivalent": not differences,
        "differences": differences,
        "base": base,
        "candidate": candidate,
    }
    if args.json:
        args.json.write_text(json.dumps(result, indent=2))
    print(json.dumps(result), flush=True)
    raise SystemExit(bool(differences))


if __name__ == "__main__":
    main()
