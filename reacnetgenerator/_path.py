# cython: language_level=3
# cython: linetrace=True
"""Collect paths.

To produce a reaction network, every molecule (species) should be treated as a
node in the network. Therefore, all detected species are indexed by canonical
SMILES to guarantee its uniqueness. Isomers are also identified according to
SMILES codes. The VF2 algorithm can be also used to identify isomers, which is
an option in ReacNetGen. After filtering out noise, the reaction path of atoms
and the number of intermolecular reactions can be calculated.

Reference:
[1] Landrum, G. RDKit: Open-Source Cheminformatics Software 2016.
[2] Cordella, L. P.; Foggia, P.; Sansone, C.; Vento, M. A (Sub)Graph
Isomorphism Algorith for Matching Large Graphs. IEEE Trans. Pattern Analysis
and Machine Intelligence 2004, 26, 1367-1372.
"""

import itertools
from abc import ABCMeta, abstractmethod
from collections import Counter, defaultdict

import networkx as nx
import networkx.algorithms.isomorphism as iso
import numpy as np
from rdkit import Chem
from tqdm import tqdm

from ._reaction import ReactionsFinder
from .utils import WriteBuffer, bytestolist, listtostirng, run_mp, SharedRNGData


class _CollectPaths(SharedRNGData, metaclass=ABCMeta):
    def __init__(self, rng):
        SharedRNGData.__init__(self, rng, ["runHMM", "N", "step", "atomname", "originfilename", "hmmfilename", "moleculefilename",
                                           "moleculetemp2filename", "atomroutefilename", "nproc", "hmmit", "atomtype",
                                           "selectatoms", "split"], ["mname", "atomnames", "allmoleculeroute", "splitmoleculeroute"])

    @staticmethod
    def getstype(rng):
        """Get a class for different methods.

        Following methonds are used to identify isomers:
        * SMILES (default)
        * VF2
        """
        if rng.SMILES:
            return _CollectSMILESPaths(rng)
        return _CollectMolPaths(rng)

    def collect(self):
        """Collect paths."""
        self.atomnames = self.atomname[self.atomtype]
        self._printmoleculename()
        atomeach, conflict = self._getatomeach()
        self.allmoleculeroute = self._printatomroute(atomeach)
        if self.split > 1:
            splittime = np.array_split(np.arange(self.step), self.split)
            self.splitmoleculeroute = list([self._printatomroute(
                atomeach[:, st], timeaxis=i) for i, st in enumerate(splittime)])
        self.returnkeys()
        ReactionsFinder(self.rng).findreactions(atomeach.T, conflict.T)

    @abstractmethod
    def _printmoleculename(self):
        pass

    def _getatomeach(self):
        """Values in atomeach starts from 1."""
        atomeach = np.zeros((self.N, self.step), dtype=int)
        conflict = np.zeros((self.N, self.step), dtype=int)
        with open(self.hmmfilename if self.runHMM else self.originfilename, 'rb') as fh, open(self.moleculetemp2filename, 'rb') as ft:
            for i, (linehz, linetz) in enumerate(tqdm(zip(fh, itertools.zip_longest(*[ft] * 3)),
                                                      total=self.hmmit, desc="Analyze atoms", unit="molecule"), start=1):
                lineh = bytestolist(linehz)
                atom = np.array(bytestolist(linetz[0]))
                index = np.where(lineh)[0]
                if index.size:
                    conflict[np.nonzero(atomeach[atom[:, None], index])] = 1
                    atomeach[atom[:, None], index] = i
        return atomeach, conflict

    def _getatomroute(self, item):
        i, (atomeachi, atomtypei) = item
        atomeachi = atomeachi[np.nonzero(atomeachi)[0]]
        route = atomeachi[np.concatenate([[0], np.nonzero(np.diff(atomeachi))[
                                         0]+1])] if atomeachi.size else np.zeros(0, dtype=int)
        moleculeroute = np.dstack((route[:-1], route[1:]))[
            0] if self.atomname[atomtypei] in self.selectatoms else np.zeros((0, 2), dtype=int)
        names = self.mname[route-1]
        routestr = f"Atom {i} {self.atomname[atomtypei]}: " + \
            " -> ".join(names)
        return moleculeroute, routestr

    def _printatomroute(self, atomeach, timeaxis=None):
        """For analysis without HMM, we may not need to use np.unique"""
        with WriteBuffer(open(self.atomroutefilename if timeaxis is None else f"{self.atomroutefilename}.{timeaxis}", 'w'), sep='\n') as f:
            allmoleculeroute = []
            if not self.runHMM:
                have_added = {}
            results = run_mp(self.nproc, func=self._getatomroute, l=zip(atomeach, self.atomtype), return_num=True, start=1, unordered=False,
                                total=self.N, desc="Collect reaction paths" if timeaxis is None else f"Collect reaction paths {timeaxis}", unit="atom")
            for ii, (moleculeroute, routestr) in enumerate(results):
                f.append(routestr)
                if moleculeroute.size > 0:
                    if not self.runHMM:
                        # check whether repeated or not if analyzing without HMM
                        for rr in moleculeroute:
                            tpr = tuple(rr)
                            if have_added.get(tpr, atomeach.shape[0]) >= ii:
                                have_added[tpr] = ii
                                allmoleculeroute.append(rr.reshape(1,2))
                    else:
                        allmoleculeroute.append(moleculeroute)
        allmoleculeroute = np.concatenate(allmoleculeroute) if allmoleculeroute else np.zeros((0, 2), dtype=int)
        if self.runHMM and allmoleculeroute.size:
            allmoleculeroute = np.unique(allmoleculeroute, axis=0)
        return allmoleculeroute

    def convertSMILES(self, atoms, bonds):
        """Convert atoms and bonds information to SMILES."""
        m = Chem.RWMol(Chem.MolFromSmiles(''))
        d = {}
        for name, number in zip(self.atomnames[atoms], atoms):
            d[number] = m.AddAtom(Chem.Atom(name))
        for atom1, atom2, level in bonds:
            m.AddBond(d[atom1], d[atom2], Chem.BondType(level))
        name = Chem.MolToSmiles(m)
        return name

    def _getatomsandbonds(self, line):
        atoms = np.array(bytestolist(line[0]), dtype=int)
        bonds = bytestolist(line[1])
        return atoms, bonds


class _CollectMolPaths(_CollectPaths):
    def _printmoleculename(self):
        mname = []
        d = defaultdict(list)
        em = iso.numerical_edge_match(['atom', 'level'], ["None", 1])
        with WriteBuffer(open(self.moleculefilename, 'w'), sep='\n') as fm, open(self.moleculetemp2filename, 'rb') as ft:
            for line in itertools.zip_longest(*[ft] * 3):
                atoms, bonds = self._getatomsandbonds(line)
                molecule = self._molecule(self, atoms, bonds)
                for isomer in d[str(molecule)]:
                    if isomer.isomorphic(molecule, em):
                        molecule.smiles = isomer.smiles
                        break
                else:
                    d[str(molecule)].append(molecule)
                mname.append(molecule.smiles)
                fm.append(listtostirng(
                    (molecule.smiles, atoms, bonds), sep=(' ', ';', ',')))
        self.mname = np.array(mname)

    class _molecule:
        def __init__(self, cmp, atoms, bonds):
            self.atoms = atoms
            self.bonds = bonds
            self._atomtypes = cmp.atomtype[atoms]
            self._atomnames = cmp.atomnames[atoms]
            self.graph = self._makemoleculegraph()
            counter = Counter(self._atomnames)
            self.name = "".join(map(lambda atomname: f"{atomname}{counter[atomname]}",
                                    cmp.atomname))
            self._smiles = None
            self._convertSMILES = cmp.convertSMILES

        def __str__(self):
            return self.name

        @property
        def smiles(self):
            """Return SMILES of a molecule."""
            if self._smiles is None:
                self._smiles = self._convertSMILES(self.atoms, self.bonds)
            return self._smiles

        @smiles.setter
        def smiles(self, value):
            self._smiles = value

        def _makemoleculegraph(self):
            graph = nx.Graph()
            for line in self.bonds:
                graph.add_edge(line[0], line[1], level=line[2])
            for atomnumber, atomtype in zip(self.atoms, self._atomtypes):
                graph.add_node(atomnumber, atom=atomtype)
            return graph

        def isomorphic(self, mol, em):
            """Return whether two molecules are isomorphic."""
            return nx.is_isomorphic(self.graph, mol.graph, em)


class _CollectSMILESPaths(_CollectPaths):
    def _printmoleculename(self):
        mname = []
        with WriteBuffer(open(self.moleculefilename, 'w'), sep='\n') as fm, open(self.moleculetemp2filename, 'rb') as ft:
            results = run_mp(self.nproc, func=self._calmoleculeSMILESname, l=ft, unordered=False,
                                nlines=3, total=self.hmmit, desc="Indentify isomers", unit="molecule")
            for name, atoms, bonds in results:
                mname.append(name)
                fm.append(listtostirng(
                    (name, atoms, bonds), sep=(' ', ';', ',')))
        self.mname = np.array(mname)

    def _calmoleculeSMILESname(self, item):
        line = item
        atoms, bonds = self._getatomsandbonds(line)
        name = self.convertSMILES(atoms, bonds)
        return name, atoms, bonds
