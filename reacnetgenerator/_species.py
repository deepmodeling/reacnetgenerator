# SPDX-License-Identifier: LGPL-3.0-or-later
"""Build species graphs without changing the underlying trajectory records."""

import networkx as nx
from networkx.algorithms import isomorphism as iso
from rdkit import Chem, rdBase

# Limit perception to ordinary organic chemistry. In particular, a missing
# metal/ligand charge or coordination bond must not trigger charge/bond cleanup.
_AROMATIC_ELEMENTS = frozenset((1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 52, 53))
_AROMATIC_BOND_TYPES = frozenset(
    (
        Chem.BondType.SINGLE,
        Chem.BondType.DOUBLE,
        Chem.BondType.TRIPLE,
        Chem.BondType.AROMATIC,
    )
)
_AROMATIC_SANITIZE_OPS = (
    Chem.SanitizeFlags.SANITIZE_PROPERTIES
    | Chem.SanitizeFlags.SANITIZE_SYMMRINGS
    | Chem.SanitizeFlags.SANITIZE_KEKULIZE
    | Chem.SanitizeFlags.SANITIZE_FINDRADICALS
    | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
)


_ATOM_MATCH = iso.categorical_node_match("atom", None)
_BOND_MATCH = iso.categorical_edge_match("level", None)


def molecule_graph(atoms, atomtypes, bonds, miso):
    """Apply the selected atom-label and bond policy for VF2 comparison."""
    graph = nx.Graph()
    if miso not in (0, 1, 2):
        raise ValueError(f"Unknown isomer identification method: {miso}.")
    if miso != 2:
        for bond in bonds:
            graph.add_edge(bond[0], bond[1], level=bond[2] if miso == 0 else 1)
    for atom, atomtype in zip(atoms, atomtypes, strict=True):
        graph.add_node(atom, atom=atomtype)
    return graph


def isomorphic(left, right):
    """Compare graphs built under the same isomer policy."""
    return nx.is_isomorphic(left, right, node_match=_ATOM_MATCH, edge_match=_BOND_MATCH)


def species_molecule(atomnames, atoms, bonds):
    """Build a molecule with equivalent Kekule forms sharing aromatic bonds.

    Atom names follow the order of ``atoms``; bond endpoints use trajectory
    indices. All hydrogens must be explicit. Only aromatic flags and aromatic
    bond orders are transferred from a successfully validated copy. Charge
    cleanup, hydrogen adjustment and tautomer enumeration are never requested.
    Unsupported chemistry or failed perception retains the original graph.
    """
    mol = Chem.RWMol()
    indices = {}
    for name, number in zip(atomnames, atoms):
        atom = Chem.Atom(name)
        # Preserve radicals and missing hydrogens as in the trajectory.
        atom.SetNoImplicit(True)
        indices[number] = mol.AddAtom(atom)
    for atom1, atom2, level in bonds:
        mol.AddBond(indices[atom1], indices[atom2], Chem.BondType(int(level)))

    # Avoid sanitization/copying for acyclic species, including common small
    # reactive intermediates. Ring discovery only sets computed properties.
    Chem.FastFindRings(mol)
    if not mol.GetRingInfo().NumRings():
        return mol
    if any(atom.GetAtomicNum() not in _AROMATIC_ELEMENTS for atom in mol.GetAtoms()):
        return mol
    if any(bond.GetBondType() not in _AROMATIC_BOND_TYPES for bond in mol.GetBonds()):
        return mol

    candidate = Chem.Mol(mol)
    # A failed sanitization can partially mutate its argument. Discard that
    # copy, and keep the serializable raw graph instead of losing the record.
    # Such failures are expected for transient/charge-incomplete MD structures.
    with rdBase.BlockLogs():
        failed = Chem.SanitizeMol(
            candidate, sanitizeOps=_AROMATIC_SANITIZE_OPS, catchErrors=True
        )
    if failed:
        return mol
    for bond in candidate.GetBonds():
        if bond.GetIsAromatic():
            original = mol.GetBondWithIdx(bond.GetIdx())
            # RDKit can mark an aryne triple bond aromatic while retaining its
            # triple bond type. The aromatic flag alone does not imply AROMATIC.
            original.SetBondType(bond.GetBondType())
            original.SetIsAromatic(True)
            original.GetBeginAtom().SetIsAromatic(True)
            original.GetEndAtom().SetIsAromatic(True)
    return mol
