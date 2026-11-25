# SPDX-License-Identifier: LGPL-3.0-or-later
"""Exercise species descriptors through the bond-input and HDF5 interfaces."""

import h5py
import pytest
from rdkit import Chem

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._speciesdescriptors import hill_formula, ring_system
from reacnetgenerator.timedoutput import validate_timed_output
from tests.test_aromatic_species import _records, _write_bond_trajectory


@pytest.mark.parametrize("miso", [0, 1, 2])
def test_species_descriptor_columns(tmp_path, miso):
    """Keep composition and topology aligned under each identity policy."""
    # Propene and cyclopropane have the same explicit atom counts. The third
    # component exercises the central five-membered ring of dibenzofuran.
    mol = Chem.AddHs(Chem.MolFromSmiles("CC=C.C1CC1.c1ccc2oc3ccccc3c2c1"))
    Chem.Kekulize(mol, clearAromaticFlags=True)
    symbols, bonds = _records(mol)
    source = tmp_path / "rings.bond"
    types = _write_bond_trajectory(source, symbols, [bonds, bonds])
    path = tmp_path / "timeline.h5"
    rng = ReacNetGenerator(
        inputfilename=source,
        inputfiletype="bond",
        atomname=types,
        miso=miso,
        runHMM=False,
        nproc=1,
        output_dir=tmp_path,
        timed_output=path,
    )
    rng.run_items(["species", "reactions"])
    assert validate_timed_output(path).schema_version == "1.2"
    with h5py.File(path) as file:
        formulas = file["species/formula"].asstr()[:]
        if miso == 2:
            assert sorted(formulas) == ["C12H8O", "C3H6"]
            assert "species/ring_system" not in file
        else:
            rings = file["species/ring_system"].asstr()[:]
            assert sorted(zip(formulas, rings, strict=True)) == [
                ("C12H8O", "fused-6-5-6"),
                ("C3H6", "acyclic"),
                ("C3H6", "monocyclic-3"),
            ]


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCCC1C2CCCC2", "linked-5-5"),
        ("C1CCC2(CC1)CCCC2", "spiro-6-5"),
        ("C1CC2CCC1C2", "polycyclic"),
        ("C1CCC2CCCC2C1", "fused-6-5"),
    ],
)
def test_ring_system_is_independent_of_atom_numbering(smiles, expected):
    """Describe ring sharing consistently after reversing trajectory indices."""
    mol = Chem.MolFromSmiles(smiles)
    for order in (
        list(range(mol.GetNumAtoms())),
        list(reversed(range(mol.GetNumAtoms()))),
    ):
        symbols, bonds = _records(Chem.RenumberAtoms(mol, order))
        assert ring_system(len(symbols), bonds[:, :2]) == expected


def test_hill_order_with_and_without_carbon():
    """Put hydrogen after carbon only in carbon-containing formulas."""
    assert hill_formula(["Cl", "H", "Cl"]) == "Cl2H"
    assert hill_formula(["Cl", "H", "C", "H"]) == "CH2Cl"
