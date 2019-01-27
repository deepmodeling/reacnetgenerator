''' Collect paths:
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
'''

from abc import ABCMeta, abstractmethod
from multiprocessing import Pool, Semaphore
from collections import Counter, defaultdict

import numpy as np
from tqdm import tqdm
import networkx as nx
import networkx.algorithms.isomorphism as iso
from rdkit import Chem


class _CollectPaths(metaclass=ABCMeta):
    def __init__(self, rng):
        self.rng = rng
        self.runHMM = rng.runHMM
        self._N = rng.N
        self._step = rng.step
        self._timestep = rng.timestep
        self.atomname = rng.atomname
        self.originfilename = rng.originfilename
        self.hmmfilename = rng.hmmfilename
        self.moleculefilename = rng.moleculefilename
        self.moleculetemp2filename = rng.moleculetemp2filename
        self.atomroutefilename = rng.atomroutefilename
        self.nproc = rng.nproc
        self._hmmit = rng.hmmit
        self._atomtype = rng.atomtype
        self.selectatoms = rng.selectatoms
        self._decompress = rng.decompress
        self._produce = rng.produce
        self._mname = None
        self._atomnames = None

    @staticmethod
    def getstype(SMILES):
        if SMILES:
            return _CollectSMILESPaths
        return _CollectMolPaths

    def collect(self):
        ''' Collect paths'''
        self._atomnames = self.atomname[self._atomtype-1]
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
            for i, (linehz, linetz) in enumerate(zip(fh, ft), start=1):
                lineh = self._decompress(linehz)
                linet = self._decompress(linetz)
                s = linet.split()
                key1 = np.array([int(x) for x in s[0].split(",")])
                index = np.array(
                    [j for j in range(len(lineh)) if lineh[j] == "1"])
                if index.size:
                    atomeach[key1[:, None]-1, index] = i
        return atomeach

    def _getatomroute(self, item):
        (i, (atomeachi, atomtypei)), _ = item
        routestrarr = []
        moleculeroute = []
        right = -1
        for j, atomeachij in enumerate(atomeachi.tolist()):
            if atomeachij > 0 and atomeachij != right:
                routestrarr.append(
                    f"{self._mname[atomeachij-1]} ({atomeachij} step { self._timestep[j]})")
                left, right = right, atomeachij
                if self.atomname[atomtypei-1] in self.selectatoms:
                    if left >= 0 and (left, right) not in moleculeroute:
                        moleculeroute.append((left, right))
        routestr = f"Atom {i} {self.atomname[atomtypei-1]}: " + \
            " -> ".join(routestrarr)
        return moleculeroute, routestr

    def _printatomroute(self, atomeach):
        with open(self.atomroutefilename, 'w') as f, Pool(self.nproc, maxtasksperchild=1000) as pool:
            allmoleculeroute = []
            semaphore = Semaphore(self.nproc*150)
            results = pool.imap(self._getatomroute, self._produce(
                semaphore, enumerate(zip(atomeach, self._atomtype), start=1), ()), 100)
            for route in tqdm(results, total=self._N, desc="Collect reaction paths", unit="atom"):
                moleculeroute, routestr = route
                print(routestr, file=f)
                for mroute in moleculeroute:
                    if mroute not in allmoleculeroute:
                        allmoleculeroute.append(mroute)
                semaphore.release()
        return allmoleculeroute

    def _convertSMILES(self, atoms, bonds):
        m = Chem.RWMol(Chem.MolFromSmiles(''))
        d = {}
        for name, number in zip(self._atomnames[atoms-1], atoms):
            d[number] = m.AddAtom(Chem.Atom(name))
        for atom1, atom2, level in bonds:
            m.AddBond(d[atom1], d[atom2], Chem.BondType(level))
        name = Chem.MolToSmiles(m)
        return name


class _CollectMolPaths(_CollectPaths):
    def _printmoleculename(self):
        mname = []
        d = defaultdict(list)
        em = iso.numerical_edge_match(['atom', 'level'], ["None", 1])
        with open(self.moleculefilename, 'w') as fm, open(self.moleculetemp2filename, 'rb') as ft:
            for line in ft:
                s = self._decompress(line).split()
                atoms = np.array([int(x) for x in s[0].split(",")])
                bonds = tuple(tuple(int(y) for y in x.split(","))
                              for x in s[1].split(";")) if len(s) == 3 else ()
                molecule = self._molecule(self, atoms, bonds)
                for isomer in d[str(molecule)]:
                    if isomer.isomorphic(molecule, em):
                        molecule.smiles = isomer.smiles
                        break
                else:
                    d[str(molecule)].append(molecule)
                mname.append(molecule.smiles)
                print(molecule.smiles, ",".join([str(x) for x in atoms]), ";".join(
                    [",".join([str(y) for y in x]) for x in bonds]), file=fm)
        self._mname = mname

    class _molecule:
        def __init__(self, cmp, atoms, bonds):
            self.atoms = atoms
            self.bonds = bonds
            self._atomtypes = cmp._atomtype[atoms-1]
            self._atomnames = cmp._atomnames[atoms-1]
            self.graph = self._makemoleculegraph()
            counter = Counter(self._atomnames)
            self.name = "".join(
                [f"{atomname}{counter[atomname]}" for atomname in cmp.atomname])
            self._smiles = None
            self._convertSMILES = cmp._convertSMILES

        def __str__(self):
            return self.name

        @property
        def smiles(self):
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
            return nx.is_isomorphic(self.graph, mol.graph, em)


class _CollectSMILESPaths(_CollectPaths):
    def _printmoleculename(self):
        mname = []
        with open(self.moleculefilename, 'w') as fm, open(self.moleculetemp2filename, 'rb') as ft, Pool(self.nproc, maxtasksperchild=1000) as pool:
            semaphore = Semaphore(self.nproc*150)
            results = pool.imap(self._calmoleculeSMILESname,
                                self._produce(semaphore, ft, None), 100)
            for name, atoms, bonds in tqdm(results, total=self._hmmit, desc="Indentify isomers", unit="molecule"):
                mname.append(name)
                fm.write(' '.join((name, ",".join((str(x) for x in atoms)), ";".join(
                    (",".join((str(y) for y in x)) for x in bonds)), '\n')))
                semaphore.release()
        self._mname = mname

    def _calmoleculeSMILESname(self, item):
        line, _ = item
        s = self._decompress(line).split()
        atoms = np.array([int(x) for x in s[0].split(",")])
        bonds = tuple(tuple(int(y) for y in x.split(","))
                      for x in s[1].split(";")) if len(s) == 3 else ()
        name = self._convertSMILES(atoms, bonds)
        return name, atoms, bonds
