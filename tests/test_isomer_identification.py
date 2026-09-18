# SPDX-License-Identifier: LGPL-3.0-or-later
"""Tests for VF2 isomer matching and ``miso`` validation."""

import itertools
from types import SimpleNamespace

import numpy as np
import pytest

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._path import Molecule, _CollectMolPaths, _CollectSMILESPaths
from reacnetgenerator.commandline import main_parser
from reacnetgenerator.utils import listtobytes


def _molecule_context(miso, atomtypes, atomnames):
    return SimpleNamespace(
        atomtype=np.asarray(atomtypes),
        atomnames=np.asarray(atomnames),
        atomname=list(dict.fromkeys(atomnames)),
        miso=miso,
        n_unknown=0,
        convertSMILES=lambda atoms, bonds: f"bond-order-{bonds[0][2]}",
    )


def _molecule(context, atoms, bonds):
    return Molecule(context, np.asarray(atoms), np.asarray(bonds, dtype=int))


def test_miso_zero_distinguishes_bond_orders_and_atom_labels():
    """Normal VF2 matching should compare both node and edge attributes."""
    bond_context = _molecule_context(0, [0, 0], ["C", "C"])
    single = _molecule(bond_context, [0, 1], [[0, 1, 1]])
    double = _molecule(bond_context, [0, 1], [[0, 1, 2]])

    label_context = _molecule_context(0, [0, 1, 0, 0], ["C", "O", "C", "C"])
    carbon_oxygen = _molecule(label_context, [0, 1], [[0, 1, 1]])
    carbon_carbon = _molecule(label_context, [2, 3], [[2, 3, 1]])

    assert not single.isomorphic(double)
    assert not carbon_oxygen.isomorphic(carbon_carbon)


def test_miso_zero_matches_permuted_atom_indices():
    """VF2 matching should remain independent of trajectory atom indices."""
    context = _molecule_context(0, [0, 1, 1, 0], ["C", "O", "O", "C"])
    first = _molecule(context, [0, 1], [[0, 1, 2]])
    permuted = _molecule(context, [2, 3], [[2, 3, 2]])

    assert first.isomorphic(permuted)


def test_miso_one_ignores_bond_orders_but_preserves_labeled_connectivity():
    """Mode 1 should merge only structures differing in bond order."""
    bond_context = _molecule_context(1, [0, 0], ["C", "C"])
    single = _molecule(bond_context, [0, 1], [[0, 1, 1]])
    double = _molecule(bond_context, [0, 1], [[0, 1, 2]])

    label_context = _molecule_context(
        1,
        [0, 1, 0, 0, 0, 1],
        ["C", "O", "C", "C", "C", "O"],
    )
    oxygen_centered = _molecule(label_context, [0, 1, 2], [[0, 1, 1], [1, 2, 1]])
    carbon_centered = _molecule(label_context, [3, 4, 5], [[3, 4, 1], [4, 5, 1]])

    assert single.isomorphic(double)
    assert not oxygen_centered.isomorphic(carbon_centered)


def test_miso_two_compares_composition_only():
    """Mode 2 should ignore connectivity while retaining atom composition."""
    context = _molecule_context(
        2,
        [0, 0, 1, 0, 0, 1, 0, 0, 0],
        ["C", "C", "O", "C", "C", "O", "C", "C", "C"],
    )
    chain = _molecule(context, [0, 1, 2], [[0, 1, 1], [1, 2, 1]])
    triangle = _molecule(
        context,
        [3, 4, 5],
        [[3, 4, 1], [4, 5, 1], [3, 5, 1]],
    )
    different_composition = _molecule(
        context,
        [6, 7, 8],
        [[6, 7, 1], [7, 8, 1]],
    )

    assert chain.isomorphic(triangle)
    assert not chain.isomorphic(different_composition)


def _write_molecule_records(path, records):
    with path.open("wb") as handle:
        for atoms, bonds in records:
            handle.write(listtobytes(atoms))
            handle.write(listtobytes([(atom1, atom2) for atom1, atom2, _ in bonds]))
            handle.write(listtobytes([level for _, _, level in bonds]))
            handle.write(listtobytes([0]))


def _collector_context(tmp_path, *, convert_smiles, miso=0):
    molecule_file = tmp_path / "molecules.bin"
    _write_molecule_records(
        molecule_file,
        [
            ([0, 1], [(0, 1, 1)]),
            ([0, 1], [(0, 1, 2)]),
        ],
    )
    return SimpleNamespace(
        hmmit=2,
        moleculefilename=str(tmp_path / "molecule-names.txt"),
        moleculetimelinefilename=str(tmp_path / "timeline.csv"),
        moleculetemp2filename=str(molecule_file),
        atomtype=np.asarray([0, 0]),
        atomnames=np.asarray(["C", "C"]),
        atomname=["C"],
        miso=miso,
        convertSMILES=convert_smiles,
        printmoleculetime=False,
        _moleculeframefilter=None,
        _moleculetimestepfilter=None,
        nproc=1,
    )


def _run_serial(_nproc, *, func, l, nlines, **kwargs):
    iterator = iter(l)
    while group := tuple(itertools.islice(iterator, nlines)):
        assert len(group) == nlines
        yield func(group)


def test_direct_vf2_collector_keeps_different_bond_orders(tmp_path):
    """The non-SMILES collector should not merge distinct bond orders."""
    context = _collector_context(
        tmp_path,
        convert_smiles=lambda atoms, bonds: f"bond-order-{bonds[0][2]}",
    )
    collector = object.__new__(_CollectMolPaths)
    collector.__dict__.update(vars(context))

    collector._printmoleculename()

    assert list(collector.mname) == ["bond-order-1", "bond-order-2"]


def test_smiles_failure_fallback_keeps_different_bond_orders(tmp_path, monkeypatch):
    """SMILES failures should fall back to bond-order-aware VF2 matching."""

    def fail_smiles(atoms, bonds):
        raise ValueError("forced SMILES failure")

    monkeypatch.setattr("reacnetgenerator._path.run_mp", _run_serial)
    context = _collector_context(tmp_path, convert_smiles=fail_smiles)
    collector = object.__new__(_CollectSMILESPaths)
    collector.__dict__.update(vars(context))

    collector._printmoleculename()

    assert list(collector.mname) == ["C2_unknownSMILES_0", "C2_unknownSMILES_1"]


def test_smiles_collector_miso_one_merges_different_bond_orders(tmp_path, monkeypatch):
    """Successful SMILES collection should retain mode 1 merging semantics."""
    monkeypatch.setattr("reacnetgenerator._path.run_mp", _run_serial)
    context = _collector_context(
        tmp_path,
        convert_smiles=lambda atoms, bonds: f"bond-order-{bonds[0][2]}",
        miso=1,
    )
    collector = object.__new__(_CollectSMILESPaths)
    collector.__dict__.update(vars(context))

    collector._printmoleculename()

    assert list(collector.mname) == ["bond-order-1", "bond-order-1"]


@pytest.mark.parametrize("miso", [-1, 3])
def test_python_api_rejects_invalid_miso(miso, tmp_path):
    """The Python entry point should reject invalid modes immediately."""
    with pytest.raises(ValueError, match="miso must be one of 0, 1, or 2"):
        ReacNetGenerator(
            inputfilename=str(tmp_path / "bondless.bond"),
            inputfiletype="lammpsbondfile",
            atomname=["C"],
            miso=miso,
        )


def test_python_api_preserves_default_miso(tmp_path):
    """An omitted or explicit ``None`` mode should retain the default."""
    kwargs = {
        "inputfilename": str(tmp_path / "bondless.bond"),
        "inputfiletype": "lammpsbondfile",
        "atomname": ["C"],
    }

    assert ReacNetGenerator(**kwargs).miso == 0
    assert ReacNetGenerator(**kwargs, miso=None).miso == 0


@pytest.mark.parametrize("miso", ["-1", "3"])
def test_cli_rejects_invalid_miso(miso):
    """The CLI should reject unsupported modes during argument parsing."""
    parser = main_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["-i", "bondless.bond", "-a", "C", "--miso", miso])


@pytest.mark.parametrize("miso", ["0", "1", "2"])
def test_cli_accepts_supported_miso(miso):
    """The CLI should continue to accept each documented mode."""
    parser = main_parser()

    args = parser.parse_args(["-i", "bondless.bond", "-a", "C", "--miso", miso])

    assert args.miso == int(miso)


@pytest.mark.parametrize("miso", [-1, 3])
def test_bondless_molecule_rejects_invalid_miso(miso):
    """Defensive molecule validation should not depend on entering the bond loop."""
    context = _molecule_context(miso, [0], ["C"])

    with pytest.raises(ValueError, match="Unknown isomer identification method"):
        _molecule(context, [0], np.empty((0, 3), dtype=int))
