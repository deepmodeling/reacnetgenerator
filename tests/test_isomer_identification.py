# SPDX-License-Identifier: LGPL-3.0-or-later
"""Tests for VF2 isomer matching and ``miso`` validation."""

import io
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._path import Molecule, _CollectMolPaths, _CollectSMILESPaths
from reacnetgenerator.commandline import main_parser
from reacnetgenerator.utils import listtobytes, run_mp


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


def _unexpected_pool(*args, **kwargs):
    raise AssertionError("Single-process conversion must not start workers")


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

    monkeypatch.setattr("reacnetgenerator._path.run_mp", _unexpected_pool)
    context = _collector_context(tmp_path, convert_smiles=fail_smiles)
    collector = object.__new__(_CollectSMILESPaths)
    collector.__dict__.update(vars(context))

    collector._printmoleculename()

    assert list(collector.mname) == ["C2_unknownSMILES_0", "C2_unknownSMILES_1"]


def test_smiles_collector_miso_one_merges_different_bond_orders(tmp_path, monkeypatch):
    """Successful SMILES collection should retain mode 1 merging semantics."""
    monkeypatch.setattr("reacnetgenerator._path.run_mp", _unexpected_pool)
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


@pytest.mark.parametrize("nproc", [1, 2])
@pytest.mark.parametrize("miso", [0, 1, 2])
def test_smiles_entry_point_preserves_order(tmp_path, monkeypatch, nproc, miso):
    """Real conversion uses workers only for explicit parallel configurations."""
    context = _collector_context(tmp_path, convert_smiles=None, miso=miso)
    collector = object.__new__(_CollectSMILESPaths)
    collector.__dict__.update(vars(context))
    del collector.convertSMILES
    collector.nproc = nproc
    calls = []

    def workers(count, **kwargs):
        assert count == 2
        calls.append(count)
        return run_mp(count, **kwargs)

    monkeypatch.setattr("reacnetgenerator._path.run_mp", workers)
    collector._printmoleculename()
    assert calls == ([2] if nproc == 2 else [])
    assert list(collector.mname) == (
        ["[C][C]", "[C]=[C]"] if miso == 0 else ["[C]=[C]", "[C]=[C]"]
    )


def test_single_process_names_are_lazy_and_ordered(tmp_path, monkeypatch):
    """Requesting one result must not convert or read subsequent records."""
    context = _collector_context(tmp_path, convert_smiles=None)
    collector = object.__new__(_CollectSMILESPaths)
    collector.__dict__.update(vars(context))
    monkeypatch.setattr("reacnetgenerator._path.run_mp", _unexpected_pool)
    stream = io.BytesIO(Path(context.moleculetemp2filename).read_bytes())
    converted = []

    def convert(record):
        converted.append(record)
        return len(converted)

    results = collector._getSMILESresults(stream, convert)
    assert not converted and stream.tell() == 0
    iterator = iter(results)
    assert next(iterator) == 1
    assert len(converted) == 1 and stream.tell() < len(stream.getvalue())
    assert list(iterator) == [2]


@pytest.mark.parametrize("miso", [0, 1, 2])
def test_local_name_failure_closes_files(tmp_path, monkeypatch, miso):
    """Non-conversion failures propagate while caller-owned streams close."""
    handles = []

    def tracked_open(*args, **kwargs):
        handle = open(*args, **kwargs)
        handles.append(handle)
        return handle

    def fail(atoms, bonds):
        raise RuntimeError("name conversion interrupted")

    context = _collector_context(tmp_path, convert_smiles=fail, miso=miso)
    collector = object.__new__(_CollectSMILESPaths)
    collector.__dict__.update(vars(context))
    monkeypatch.setattr("reacnetgenerator._path.open", tracked_open, raising=False)
    monkeypatch.setattr("reacnetgenerator._path.run_mp", _unexpected_pool)
    with pytest.raises(RuntimeError, match="name conversion interrupted"):
        collector._printmoleculename()
    assert handles and all(handle.closed for handle in handles)
    assert Path(context.moleculetemp2filename).is_file()
