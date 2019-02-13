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
from itertools import groupby
from multiprocessing import Pool, Semaphore

import networkx as nx
import networkx.algorithms.isomorphism as iso
import numpy as np
from rdkit import Chem
from tqdm import tqdm


class _CollectPaths(metaclass=ABCMeta):
    def __init__(self, rng):
        self.rng = rng
        self.runHMM = rng.runHMM
        self._N = rng.N
        self._step = rng.step
        self.atomname = rng.atomname
        self.originfilename = rng.originfilename
        self.hmmfilename = rng.hmmfilename
        self.moleculefilename = rng.moleculefilename
        self.moleculetemp2filename = rng.moleculetemp2filename
        self.atomroutefilename = rng.atomroutefilename
        self.nproc = rng.nproc
        self._hmmit = rng.hmmit
        self.atomtype = rng.atomtype
        self.selectatoms = rng.selectatoms
        self._decompress = rng.decompress
        self._bytestolist = rng.bytestolist
        self._produce = rng.produce
        self._mname = None
        self.atomnames = None

    @staticmethod
    def getstype(SMILES):
        """Get a class for different methods.

        Following methonds are used to identify isomers:
        * SMILES (default)
        * VF2
        """
        if SMILES:
            return _CollectSMILESPaths
        return _CollectMolPaths

    def collect(self):
        """Collect paths."""
        self.atomnames = self.atomname[self.atomtype]
        self._printmoleculename()
        atomeach = self._getatomeach()
        self.rng.allmoleculeroute = self._printatomroute(atomeach)
        self.rng.mname = self._mname

    @abstractmethod
    def _printmoleculename(self):
        pass

    def _getatomeach(self):
        atomeach = np.zeros((self._N, self._step), dtype=np.int)
        with open(self.hmmfilename if self.runHMM else self.originfilename, 'rb') as fh, open(self.moleculetemp2filename, 'rb') as ft:
            for i, (linehz, linetz) in enumerate(zip(fh, itertools.zip_longest(*[ft] * 3)), start=1):
                lineh = self._bytestolist(linehz, nparray=True)
                atom = np.array(self._bytestolist(linetz[0]))
                index = np.where(lineh)[0]
                if index.size:
                    atomeach[atom[:, None], index] = i
        return atomeach

    def _getatomroute(self, item):
        (i, (atomeachi, atomtypei)), _ = item
        atomeachi = atomeachi[np.nonzero(atomeachi)[0]]
        route = atomeachi[np.nonzero(np.diff(atomeachi))[0]]
        moleculeroute = np.dstack((route[:-1], route[1:]))[
            0] if self.atomname[atomtypei] in self.selectatoms else np.zeros((0, 2), dtype=int)
        names = self._mname[route-1]
        routestr = "".join(
            [f"Atom {i} {self.atomname[atomtypei]}: ", " -> ".join(names)])
        return moleculeroute, routestr

    def _printatomroute(self, atomeach):
        with open(self.atomroutefilename, 'w') as f, Pool(self.nproc, maxtasksperchild=1000) as pool:
            allmoleculeroute = []
            semaphore = Semaphore(self.nproc*150)
            results = pool.imap(self._getatomroute, self._produce(
                semaphore, enumerate(zip(atomeach, self.atomtype), start=1), ()), 100)
            for moleculeroute, routestr in tqdm(
                    results, total=self._N, desc="Collect reaction paths",
                    unit="atom"):
                f.write("".join([routestr, '\n']))
                if moleculeroute.size > 0:
                    allmoleculeroute.append(moleculeroute)
                semaphore.release()
        pool.close()
        allmoleculeroute = np.unique(np.concatenate(allmoleculeroute), axis=0) if len(
        pool.join()
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
        atoms = np.array(self._bytestolist(line[0]), dtype=int)
        bonds = self._bytestolist(line[1])
        return atoms, bonds


class _CollectMolPaths(_CollectPaths):
    def _printmoleculename(self):
        mname = []
        d = defaultdict(list)
        em = iso.numerical_edge_match(['atom', 'level'], ["None", 1])
        buff = []
        with open(self.moleculefilename, 'w') as fm, open(self.moleculetemp2filename, 'rb') as ft:
            for line in itertools.zip_longest(*[ft] * 3):
                molecule = self._molecule(self, *self._getatomsandbonds(line))
                for isomer in d[str(molecule)]:
                    if isomer.isomorphic(molecule, em):
                        molecule.smiles = isomer.smiles
                        break
                else:
                    d[str(molecule)].append(molecule)
                mname.append(molecule.smiles)
                buff.append(' '.join((molecule.smiles, ",".join(map(str, atoms)), ";".join(
                    map(lambda x: ",".join(map(str, x)), bonds)))))
            fm.write('\n'.join(buff))
        self._mname = np.array(mname)

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
        buff = []
        with open(self.moleculefilename, 'w') as fm, open(self.moleculetemp2filename, 'rb') as ft, Pool(self.nproc, maxtasksperchild=1000) as pool:
            semaphore = Semaphore(self.nproc*150)
            results = pool.imap(self._calmoleculeSMILESname,
                                self._produce(semaphore, itertools.zip_longest(*[ft] * 3), None), 100)
            for name, atoms, bonds in tqdm(
                    results, total=self._hmmit, desc="Indentify isomers",
                    unit="molecule"):
                mname.append(name)
                buff.append(' '.join(
                    (name, ",".join(map(str, atoms)),
                     ";".join((map(lambda x: ",".join(map(str, x)), bonds))))))
                semaphore.release()
            fm.write('\n'.join(buff))
        pool.close()
        self._mname = np.array(mname)
        pool.join()

    def _calmoleculeSMILESname(self, item):
        line, _ = item
        atoms, bonds = self._getatomsandbonds(line)
        name = self.convertSMILES(atoms, bonds)
        return name, atoms, bonds
