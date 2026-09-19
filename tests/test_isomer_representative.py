# SPDX-License-Identifier: LGPL-3.0-or-later
"""Regression tests for representative selection when merging isomers."""

import csv
import itertools
from types import SimpleNamespace

import numpy as np
import pytest

from reacnetgenerator._path import _CollectSMILESPaths
from reacnetgenerator.utils import listtobytes


def _bond_key(bonds):
    return tuple(tuple(int(value) for value in bond) for bond in bonds)


def _write_molecule_records(path, records):
    with path.open("wb") as handle:
        for _, atoms, bonds, frames in records:
            handle.write(listtobytes(atoms))
            handle.write(listtobytes([(a, b) for a, b, _ in bonds]))
            handle.write(listtobytes([level for _, _, level in bonds]))
            handle.write(listtobytes(frames))


def _run_serial(_nproc, *, func, l, nlines, **kwargs):
    iterator = iter(l)
    while group := tuple(itertools.islice(iterator, nlines)):
        assert len(group) == nlines
        yield func(group)


def _name_converter(names):
    def convert_smiles(atoms, bonds):
        return names[_bond_key(bonds)]

    return convert_smiles


def _collect(
    tmp_path,
    monkeypatch,
    records,
    *,
    miso,
    print_timeline=False,
    convert_smiles=None,
):
    monkeypatch.setattr("reacnetgenerator._path.run_mp", _run_serial)
    molecule_records = tmp_path / "molecules.bin"
    _write_molecule_records(molecule_records, records)
    names = {_bond_key(bonds): name for name, _, bonds, _ in records}
    if convert_smiles is None:
        convert_smiles = _name_converter(names)

    atom_count = max(int(atom) for _, atoms, _, _ in records for atom in atoms) + 1
    frames = {int(frame) for _, _, _, values in records for frame in values}
    context = SimpleNamespace(
        hmmit=len(records),
        moleculefilename=str(tmp_path / "molecule-names.txt"),
        moleculetimelinefilename=str(tmp_path / "timeline.csv"),
        moleculetemp2filename=str(molecule_records),
        atomtype=np.zeros(atom_count, dtype=int),
        atomnames=np.full(atom_count, "C", dtype=object),
        atomname=["C"],
        miso=miso,
        convertSMILES=convert_smiles,
        printmoleculetime=print_timeline,
        _moleculeframefilter=None,
        _moleculetimestepfilter=None,
        timestep={frame: 100 * frame for frame in frames},
        nproc=1,
    )
    collector = object.__new__(_CollectSMILESPaths)
    collector.__dict__.update(vars(context))
    collector._printmoleculename()
    return collector


@pytest.mark.parametrize(
    ("miso", "low_bonds", "high_bonds", "atoms"),
    [
        (1, [[0, 1, 1]], [[0, 1, 2]], [0, 1]),
        (
            2,
            [[0, 1, 1], [1, 2, 1]],
            [[0, 1, 1], [1, 2, 1], [0, 2, 1]],
            [0, 1, 2],
        ),
    ],
)
@pytest.mark.parametrize("reverse", [False, True])
def test_miso_uses_order_independent_highest_frequency_representative(
    tmp_path, monkeypatch, miso, low_bonds, high_bonds, atoms, reverse
):
    """Modes 1 and 2 should choose the highest-frequency merged isomer."""
    records = [
        ("low", atoms, low_bonds, [0]),
        ("high", atoms, high_bonds, range(1, 11)),
    ]
    if reverse:
        records.reverse()

    collector = _collect(tmp_path, monkeypatch, records, miso=miso)

    assert list(collector.mname) == ["high", "high"]


def test_miso_aggregates_frequency_across_molecule_records(tmp_path, monkeypatch):
    """Separate molecule instances with one SMILES should share their count."""
    records = [
        ("B", [0, 1], [[0, 1, 2]], range(5)),
        ("A", [0, 1], [[0, 1, 1]], range(3)),
        ("A", [2, 3], [[2, 3, 1]], range(3, 6)),
    ]

    collector = _collect(tmp_path, monkeypatch, records, miso=1)

    assert list(collector.mname) == ["A", "A", "A"]


@pytest.mark.parametrize("reverse", [False, True])
def test_miso_breaks_frequency_ties_by_species_name(tmp_path, monkeypatch, reverse):
    """Equal-frequency representatives should not depend on record order."""
    records = [
        ("B", [0, 1], [[0, 1, 2]], range(3)),
        ("A", [0, 1], [[0, 1, 1]], range(3, 6)),
    ]
    if reverse:
        records.reverse()

    collector = _collect(tmp_path, monkeypatch, records, miso=1)

    assert list(collector.mname) == ["A", "A"]


def test_miso_zero_keeps_each_canonical_smiles(tmp_path, monkeypatch):
    """Frequency selection should not change the unmerged mode."""
    records = [
        ("low", [0, 1], [[0, 1, 1]], [0]),
        ("high", [0, 1], [[0, 1, 2]], range(1, 11)),
    ]

    collector = _collect(tmp_path, monkeypatch, records, miso=0, print_timeline=True)

    assert list(collector.mname) == ["low", "high"]
    with (tmp_path / "timeline.csv").open(newline="") as handle:
        assert [row["Species"] for row in csv.DictReader(handle)] == [
            "low",
            *("high" for _ in range(10)),
        ]


@pytest.mark.parametrize("miso", [0, 1])
def test_smiles_failure_reuses_vf2_fallback_name(tmp_path, monkeypatch, miso):
    """Both merged and unmerged paths should retain the VF2 fallback."""
    records = [
        ("unused", [0, 1], [[0, 1, 1]], [0]),
        ("unused", [2, 3], [[2, 3, 1]], [1, 2]),
    ]

    def fail_smiles(atoms, bonds):
        raise ValueError("forced SMILES failure")

    collector = _collect(
        tmp_path,
        monkeypatch,
        records,
        miso=miso,
        convert_smiles=fail_smiles,
    )

    assert list(collector.mname) == ["C2_unknownSMILES_0", "C2_unknownSMILES_0"]


@pytest.mark.parametrize(
    ("miso", "failed_bonds", "winner_bonds", "atoms"),
    [
        (
            1,
            ([[0, 1, 1]], [[0, 1, 2]]),
            [[0, 1, 3]],
            [0, 1],
        ),
        (
            2,
            (
                [[0, 1, 1], [1, 2, 1]],
                [[0, 1, 1], [1, 2, 1], [0, 2, 1]],
            ),
            [[0, 1, 2], [1, 2, 1]],
            [0, 1, 2],
        ),
    ],
)
@pytest.mark.parametrize("reverse", [False, True])
def test_failed_smiles_candidates_keep_individual_frequencies(
    tmp_path,
    monkeypatch,
    miso,
    failed_bonds,
    winner_bonds,
    atoms,
    reverse,
):
    """Distinct failed candidates must not pool counts before merging."""
    records = [
        ("unused-a", atoms, failed_bonds[0], range(0, 4)),
        ("unused-b", atoms, failed_bonds[1], range(4, 8)),
        ("winner", atoms, winner_bonds, range(8, 15)),
    ]
    if reverse:
        records.reverse()
    failed_keys = {_bond_key(bonds) for bonds in failed_bonds}

    def convert_smiles(atoms, bonds):
        if _bond_key(bonds) in failed_keys:
            raise ValueError("forced SMILES failure")
        return "winner"

    collector = _collect(
        tmp_path,
        monkeypatch,
        records,
        miso=miso,
        print_timeline=True,
        convert_smiles=convert_smiles,
    )

    assert list(collector.mname) == ["winner"] * 3
    assert [
        line.split(maxsplit=1)[0]
        for line in (tmp_path / "molecule-names.txt").read_text().splitlines()
    ] == ["winner"] * 3
    with (tmp_path / "timeline.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 15
    assert {row["Species"] for row in rows} == {"winner"}


def test_failed_smiles_reuses_full_structure_across_atom_ids(tmp_path, monkeypatch):
    """Repeated full structures should still aggregate after fallback."""
    records = [
        ("unused", [0, 1], [[0, 1, 1]], range(0, 3)),
        ("unused", [2, 3], [[2, 3, 1]], range(3, 6)),
        ("winner", [0, 1], [[0, 1, 2]], range(6, 11)),
    ]

    def convert_smiles(atoms, bonds):
        if bonds[0][2] == 1:
            raise ValueError("forced SMILES failure")
        return "winner"

    collector = _collect(
        tmp_path,
        monkeypatch,
        records,
        miso=1,
        convert_smiles=convert_smiles,
    )

    assert list(collector.mname) == ["C2_unknownSMILES_0"] * 3


def test_selected_representative_is_used_in_all_name_outputs(tmp_path, monkeypatch):
    """The name table, legacy file, and timeline should use one representative."""
    records = [
        ("low", [0, 1], [[0, 1, 1]], [0]),
        ("high", [0, 1], [[0, 1, 2]], [1, 2]),
    ]

    collector = _collect(tmp_path, monkeypatch, records, miso=1, print_timeline=True)

    assert list(collector.mname) == ["high", "high"]
    assert [
        line.split(maxsplit=1)[0]
        for line in (tmp_path / "molecule-names.txt").read_text().splitlines()
    ] == ["high", "high"]
    with (tmp_path / "timeline.csv").open(newline="") as handle:
        assert {row["Species"] for row in csv.DictReader(handle)} == {"high"}
