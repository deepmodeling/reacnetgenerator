# SPDX-License-Identifier: LGPL-3.0-or-later
"""Aromatic identity, conservative fallback and trajectory evidence regressions."""

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._detect import _Detect
from reacnetgenerator._hmmfilter import _HMMFilter
from reacnetgenerator._matrix import _GenerateMatrix
from reacnetgenerator._path import Molecule, _CollectPaths, _CollectSMILESPaths
from reacnetgenerator._species import species_molecule
from reacnetgenerator.timedoutput import (
    iter_molecule_ranges,
    iter_molecules,
    iter_reaction_events,
    iter_species,
    validate_timed_output,
)
from reacnetgenerator.tools import read_species


def _collector(symbols, miso=0):
    """Build only the collector state needed for SMILES and VF2 identity checks."""
    collector = object.__new__(_CollectSMILESPaths)
    collector.atomnames = np.asarray(symbols)
    collector.atomname = np.asarray(list(dict.fromkeys(symbols)))
    collector.atomtype = np.asarray(
        [list(collector.atomname).index(symbol) for symbol in symbols]
    )
    collector.miso = miso
    collector.n_unknown = 0
    return collector


def _halobenzene(shift=0, position=1):
    """Return explicit-H 1-chloro-2/3/4-fluorobenzene with either Kekule form."""
    symbols = ["C"] * 6 + ["F", "Cl"] + ["H"] * 4
    bonds = [(i, (i + 1) % 6, 2 if (i + shift) % 2 == 0 else 1) for i in range(6)]
    bonds.extend(((0, 6, 1), (position, 7, 1)))
    bonds.extend(
        (atom, hydrogen, 1)
        for hydrogen, atom in enumerate(
            (i for i in range(6) if i not in (0, position)), 8
        )
    )
    return symbols, np.asarray(bonds, dtype=int)


def _records(molecule):
    """Export atom symbols and indexed RDKit bond types for identity fixtures."""
    symbols = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    bonds = np.asarray(
        [
            (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), int(bond.GetBondType()))
            for bond in molecule.GetBonds()
        ],
        dtype=int,
    )
    return symbols, bonds


def _aryne_bond_change():
    """Return explicit-H C6H4 graphs differing in one triple/double bond order."""
    symbols = ["C"] * 6 + ["H"] * 4
    bonds = [(i, (i + 1) % 6, order) for i, order in enumerate((3, 1, 2, 1, 2, 1))]
    bonds.extend((i, i + 4, 1) for i in range(2, 6))
    aryne = np.asarray(bonds, dtype=int)
    diradical = aryne.copy()
    diradical[0, 2] = 2
    return symbols, aryne, diradical


def _identity(symbols, bonds):
    """Return SMILES and VF2 representations built from the same raw records."""
    collector = _collector(symbols)
    atoms = np.arange(len(symbols))
    return collector.convertSMILES(atoms, bonds), Molecule(collector, atoms, bonds)


def test_kekule_forms_and_permutations_have_one_identity():
    """Canonicalization must ignore Kekule choice as well as insertion order."""
    symbols, first = _halobenzene()
    _, second = _halobenzene(1)
    originals = [first.copy(), second.copy()]
    collector = _collector(symbols)
    atoms = np.arange(len(symbols))
    reference = Molecule(collector, atoms, first)
    names = set()
    random = np.random.default_rng(2510)
    for bonds in (first, second):
        for _ in range(16):
            atom_order = random.permutation(atoms)
            bond_order = random.permutation(bonds)
            names.add(collector.convertSMILES(atom_order, bond_order))
            assert reference.isomorphic(Molecule(collector, atom_order, bond_order))
    assert len(names) == 1
    name = names.pop()
    assert "[c]" in name
    np.testing.assert_array_equal(first, originals[0])
    np.testing.assert_array_equal(second, originals[1])


@pytest.mark.parametrize("smiles", ["c1ccc2ccccc2c1", "c1ccncc1", "c1cc[nH]c1"])
def test_aromatic_and_kekule_input_agree(smiles):
    """Fused rings and heteroaromatics use the same policy as benzene rings."""
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    symbols, aromatic = _records(molecule)
    Chem.Kekulize(molecule, clearAromaticFlags=True)
    _, kekule = _records(molecule)
    aromatic_name, aromatic_graph = _identity(symbols, aromatic)
    kekule_name, kekule_graph = _identity(symbols, kekule)
    assert aromatic_name == kekule_name
    assert aromatic_graph.isomorphic(kekule_graph)


def test_positional_isomers_remain_distinct():
    """Aromaticity normalization must not erase substituent positions."""
    identities = [_identity(*_halobenzene(position=position)) for position in (1, 2, 3)]
    assert len({name for name, _ in identities}) == 3
    for i, (_, graph) in enumerate(identities):
        for _, other in identities[i + 1 :]:
            assert not graph.isomorphic(other)


@pytest.mark.parametrize(
    ("left", "right"),
    [
        ("C=C", "[CH2][CH2]"),
        ("C1=CC=CC=C1", "[CH]1[CH]C=CC=C1"),
        ("CC=O", "C=CO"),
        ("[CH3]", "C"),
    ],
)
def test_bond_changes_tautomers_and_radicals_remain_distinct(left, right):
    """Same formula, genuine bond changes and hydrogen counts remain meaningful."""
    left_name, left_graph = _identity(*_records(Chem.AddHs(Chem.MolFromSmiles(left))))
    right_name, right_graph = _identity(
        *_records(Chem.AddHs(Chem.MolFromSmiles(right)))
    )
    assert left_name != right_name
    assert not left_graph.isomorphic(right_graph)


def test_heteroatom_radical_is_not_aromatized():
    """Perception must account for an N radical before assigning aromaticity."""
    symbols, bonds = _records(Chem.AddHs(Chem.MolFromSmiles("C1=C[N]C=C1")))
    molecule = species_molecule(symbols, np.arange(len(symbols)), bonds)
    assert not any(bond.GetIsAromatic() for bond in molecule.GetBonds())


def test_aromatic_flagged_triple_bond_is_preserved():
    """An aromatic flag must not erase RDKit's retained triple bond type."""
    symbols, aryne, diradical = _aryne_bond_change()
    original = aryne.copy()
    molecule = species_molecule(symbols, np.arange(len(symbols)), aryne)
    bond = molecule.GetBondBetweenAtoms(0, 1)
    assert bond.GetIsAromatic()
    assert bond.GetBondType() == Chem.BondType.TRIPLE
    aryne_name, aryne_graph = _identity(symbols, aryne)
    diradical_name, diradical_graph = _identity(symbols, diradical)
    assert aryne_name != diradical_name
    assert not aryne_graph.isomorphic(diradical_graph)
    np.testing.assert_array_equal(aryne, original)


@pytest.mark.parametrize("case", ["overvalent", "nitro", "metal", "kekulization"])
def test_unsupported_chemistry_retains_original_graph(case):
    """Failed/unsupported perception cannot leak partial changes or add charge."""
    symbols, original = _halobenzene()
    bonds = original.tolist()
    if case == "overvalent":
        symbols.append("H")
        bonds.append((0, 12, 1))
    elif case == "nitro":
        symbols[6] = "N"
        symbols.extend(("O", "O"))
        bonds.extend(((6, 12, 2), (6, 13, 2)))
    elif case == "metal":
        symbols[7] = "Cu"
    else:
        symbols = ["C"] * 5 + ["H"] * 5
        bonds = [(i, (i + 1) % 5, 12) for i in range(5)]
        bonds.extend((i, i + 5, 1) for i in range(5))
    atoms = np.arange(len(symbols))
    molecule = species_molecule(symbols, atoms, bonds)
    assert _records(molecule)[1].tolist() == [list(bond) for bond in bonds]
    assert all(atom.GetFormalCharge() == 0 for atom in molecule.GetAtoms())
    assert all(atom.GetNoImplicit() for atom in molecule.GetAtoms())
    name, graph = _identity(symbols, bonds)
    assert "_unknownSMILES" not in name
    assert graph.smiles == name
    assert [graph.graph[a][b]["level"] for a, b, _ in bonds] == [
        order for _, _, order in bonds
    ]


def test_smiles_failure_vf2_fallback_recognizes_aromatic_equivalence():
    """A serialization failure still groups equivalent aromatic graphs."""
    symbols, first = _halobenzene()
    _, second = _halobenzene(1)
    collector = _collector(symbols)
    atoms = np.arange(len(symbols))
    molecules = defaultdict(list)
    names = [
        collector._resolveSMILESname(None, atoms, bonds, molecules)
        for bonds in (first, second)
    ]
    assert names[0] == names[1]
    assert "_unknownSMILES_0" in names[0]


def _write_bond_trajectory(path, symbols, variants):
    """Write two 30-frame bond assignments in the ReaxFF bond-file format."""
    atomnames = list(dict.fromkeys(symbols))
    with path.open("w") as handle:
        for frame in range(60):
            bonds = variants[frame // 30]
            adjacency = [[] for _ in symbols]
            for atom1, atom2, order in bonds:
                adjacency[atom1].append((atom2, order))
                adjacency[atom2].append((atom1, order))
            handle.write(
                f"# Timestep {frame}\n#\n# Number of particles {len(symbols)}\n#\n"
                "# Max number of bonds per atom 3 with coarse bond order cutoff 0.300\n"
                "# Particle connection table and bond orders\n"
                "# id type nb id_1...id_nb mol bo_1...bo_nb abo nlp q\n"
            )
            for i, (symbol, neighbors) in enumerate(zip(symbols, adjacency)):
                row = [
                    i + 1,
                    atomnames.index(symbol) + 1,
                    len(neighbors),
                    *(atom + 1 for atom, _ in neighbors),
                    1,
                    *(order for _, order in neighbors),
                    sum(order for _, order in neighbors),
                    0,
                    0,
                ]
                handle.write(" ".join(map(str, row)) + "\n")
            handle.write("#\n")
    return atomnames


def _bond_set(bonds):
    """Canonicalize endpoint order when comparing retained raw bond records."""
    return frozenset((min(a, b), max(a, b), order) for a, b, order in bonds)


@pytest.mark.parametrize("smiles", [True, False])
@pytest.mark.parametrize("run_hmm", [True, False])
@pytest.mark.parametrize("miso", [0, 1, 2])
def test_kekule_switch_preserves_evidence_without_reaction(
    tmp_path, smiles, run_hmm, miso
):
    """The real pipeline uses one species ID and retains both raw bond records."""
    symbols, first = _halobenzene()
    _, second = _halobenzene(1)
    source = tmp_path / "kekule.bond"
    atomnames = _write_bond_trajectory(source, symbols, (first, second))
    timeline = tmp_path / "timeline.h5"
    rng = ReacNetGenerator(
        inputfilename=str(source),
        inputfiletype="bond",
        atomname=atomnames,
        nproc=2,
        runHMM=run_hmm,
        SMILES=smiles,
        miso=miso,
        output_dir=tmp_path,
        printreactionevent=True,
        printmoleculetime=True,
        needprintspecies=True,
        timed_output=timeline,
    )
    try:
        _Detect.gettype(rng).detect()
        _HMMFilter(rng).filter()
        _CollectPaths.getstype(rng).collect()
        _GenerateMatrix(rng).generate()
        assert rng.hmmit == 2
        assert len(set(rng.mname)) == 1
        assert Path(rng.reactionabcdfilename).read_text() == ""
        assert Path(rng.reactionfilename).read_text() == ""
        _, species_counts = read_species(rng.speciesfilename)
        assert set(species_counts) == set(rng.mname)
        np.testing.assert_array_equal(next(iter(species_counts.values())), np.ones(60))
        with open(rng.reactioneventfilename, newline="") as handle:
            assert list(csv.DictReader(handle)) == []
        with open(rng.moleculetimelinefilename, newline="") as handle:
            rows = list(csv.DictReader(handle))
        assert len(rows) == 60
        assert len({row["Species"] for row in rows}) == 1
        assert len({row["BondIDs"] for row in rows}) == 2
        assert len({row["AtomIDs"] for row in rows}) == 1

        validate_timed_output(timeline)
        assert len(list(iter_species(timeline))) == 1
        molecules = list(iter_molecules(timeline))
        assert len(molecules) == 2
        assert {molecule.species_id for molecule in molecules} == {0}
        assert {_bond_set(molecule.bonds) for molecule in molecules} == {
            _bond_set(first),
            _bond_set(second),
        }
        assert [
            (row.start_frame, row.end_frame) for row in iter_molecule_ranges(timeline)
        ] == [(0, 29), (30, 59)]
        assert list(iter_reaction_events(timeline)) == []
    finally:
        for key in (
            "moleculetempfilename",
            "moleculetemp2filename",
            "originfilename",
            "hmmfilename",
        ):
            filename = getattr(rng, key, None)
            if filename is not None:
                Path(filename).unlink(missing_ok=True)


@pytest.mark.parametrize("smiles", [True, False])
@pytest.mark.parametrize("run_hmm", [True, False])
def test_aromatic_triple_bond_change_retains_reaction(tmp_path, smiles, run_hmm):
    """Keep the bond-order transition in species IDs and all reaction outputs."""
    symbols, aryne, diradical = _aryne_bond_change()
    source = tmp_path / "aryne.bond"
    atomnames = _write_bond_trajectory(source, symbols, (aryne, diradical))
    timeline = tmp_path / "timeline.h5"
    rng = ReacNetGenerator(
        inputfilename=str(source),
        inputfiletype="bond",
        atomname=atomnames,
        nproc=1,
        runHMM=run_hmm,
        SMILES=smiles,
        miso=0,
        output_dir=tmp_path,
        printreactionevent=True,
        timed_output=timeline,
    )
    try:
        _Detect.gettype(rng).detect()
        _HMMFilter(rng).filter()
        _CollectPaths.getstype(rng).collect()
        _GenerateMatrix(rng).generate()
        assert rng.hmmit == 2
        assert len(set(rng.mname)) == 2
        with open(rng.reactioneventfilename, newline="") as handle:
            csv_events = list(csv.DictReader(handle))
        assert len(csv_events) == 1
        event = csv_events[0]
        assert int(event["Timestep_Index"]) == 29
        assert event["Reactant"] != event["Product"]
        expected_reaction = f"1 {event['Reactant']}->{event['Product']}"
        for filename in (rng.reactionabcdfilename, rng.reactionfilename):
            assert Path(filename).read_text().splitlines() == [expected_reaction]

        validate_timed_output(timeline)
        assert len(list(iter_species(timeline))) == 2
        molecules = list(iter_molecules(timeline))
        assert len(molecules) == 2
        assert len({molecule.species_id for molecule in molecules}) == 2
        assert {_bond_set(molecule.bonds) for molecule in molecules} == {
            _bond_set(aryne),
            _bond_set(diradical),
        }
        assert [
            (row.start_frame, row.end_frame) for row in iter_molecule_ranges(timeline)
        ] == [(0, 29), (30, 59)]
        hdf5_events = list(iter_reaction_events(timeline))
        assert len(hdf5_events) == 1
        assert hdf5_events[0].transition == 29
        assert hdf5_events[0].count == 1
    finally:
        for key in (
            "moleculetempfilename",
            "moleculetemp2filename",
            "originfilename",
            "hmmfilename",
        ):
            filename = getattr(rng, key, None)
            if filename is not None:
                Path(filename).unlink(missing_ok=True)
