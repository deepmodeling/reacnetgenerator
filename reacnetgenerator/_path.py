''' Collect paths '''

from abc import ABCMeta, abstractmethod
from multiprocessing import Pool, Semaphore, cpu_count

import numpy as np
from tqdm import tqdm
import networkx as nx
import networkx.algorithms.isomorphism as iso
from rdkit import Chem

class _CollectPaths(metaclass=ABCMeta):
    def __init__(self, rng):
        self.rng = rng
        self.runHMM = rng.runHMM
        self._N = rng._N
        self._step = rng._step
        self._timestep = rng._timestep
        self.atomname = rng.atomname
        self.originfilename = rng.originfilename
        self.hmmfilename = rng.hmmfilename
        self.moleculefilename = rng.moleculefilename
        self.moleculetemp2filename = rng.moleculetemp2filename
        self.atomroutefilename = rng.atomroutefilename
        self.nproc = rng.nproc
        self._hmmit = rng._hmmit
        self._atomtype = rng._atomtype
        self.selectatoms = rng.selectatoms
        self._decompress = rng._decompress
        self._produce = rng._produce
        self._mname = None

    @staticmethod
    def getstype(SMILES):
        if SMILES:
            return _CollectSMILESPaths
        return _CollectMolPaths

    def collect(self):
        ''' Collect paths'''
        self._printmoleculename()
        atomeach = self._getatomeach()
        self.rng.allmoleculeroute = self._printatomroute(atomeach)
        self.rng._mname = self._mname

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
            semaphore = Semaphore(self.nproc*15)
            results = pool.imap(self._getatomroute, self._produce(
                semaphore, enumerate(zip(atomeach, self._atomtype), start=1), ()), 10)
            for route in tqdm(results, total=self._N, desc="Collect reaction paths", unit="atom"):
                moleculeroute, routestr = route
                print(routestr, file=f)
                for mroute in moleculeroute:
                    if mroute not in allmoleculeroute:
                        allmoleculeroute.append(mroute)
                semaphore.release()
        return allmoleculeroute

class _CollectMolPaths(_CollectPaths):
    def __init__(self, rng):
        super(_CollectMolPaths, self).__init__(rng)
        self.moleculestructurefilename = rng.moleculestructurefilename

    def _printmoleculename(self):
        mname = []
        d = {}
        em = iso.numerical_edge_match(['atom', 'level'], ["None", 1])
        with open(self.moleculefilename, 'w') as fm, open(self.moleculetemp2filename, 'rb') as ft, open(self.moleculestructurefilename, 'w') as fs:
            for line in ft:
                s = self._decompress(line).split()
                atoms = np.array([int(x) for x in s[0].split(",")])
                bonds = np.array([tuple(int(y) for y in x.split(","))
                                  for x in s[1].split(";")] if len(s) == 3 else [])
                typenumber = np.zeros(len(self.atomname), dtype=np.int)
                atomtypes = []
                for atomnumber in atoms:
                    typenumber[self._atomtype[atomnumber-1]-1] += 1
                    atomtypes.append(
                        (atomnumber, self._atomtype[atomnumber-1]))
                G = self._makemoleculegraph(atomtypes, bonds)
                name = "".join([self.atomname[i]+(str(typenumber[i] if typenumber[i] > 1 else ""))
                                if typenumber[i] > 0 else "" for i in range(0, len(self.atomname))])
                if name in d:
                    for j in range(len(d[name])):
                        if nx.is_isomorphic(G, d[name][j], em):
                            if j > 0:
                                name += "_"+str(j+1)
                            break
                    else:
                        d[name].append(G)
                        name += "_"+str(len(d[name]))
                        print(self._getstructure(name, atoms, bonds), file=fs)
                else:
                    d[name] = [G]
                    print(self._getstructure(name, atoms, bonds), file=fs)
                mname.append(name)
                print(name, ",".join([str(x) for x in atoms]), ";".join(
                    [",".join([str(y) for y in x]) for x in bonds]), file=fm)
        self._mname = mname

    def _getstructure(self, name, atoms, bonds):
        index = {}
        for i, atom in enumerate(atoms, start=1):
            index[atom] = i
        return " ".join((name, ",".join([self.atomname[self._atomtype[x-1]-1] for x in atoms]), ";".join([",".join((str(index[x[0]]), str(index[x[1]]), str(x[2]))) for x in bonds])))

    @staticmethod
    def _makemoleculegraph(atoms, bonds):
        G = nx.Graph()
        for line in bonds:
            G.add_edge(line[0], line[1], level=line[2])
        for atom in atoms:
            atomnumber, atomtype = atom
            G.add_node(atomnumber, atom=atomtype)
        return G

class _CollectSMILESPaths(_CollectPaths):
    def _printmoleculename(self):
        mname = []
        with open(self.moleculefilename, 'w') as fm, open(self.moleculetemp2filename, 'rb') as ft, Pool(self.nproc, maxtasksperchild=1000) as pool:
            semaphore = Semaphore(self.nproc*15)
            results = pool.imap(self._calmoleculeSMILESname,
                                self._produce(semaphore, ft, ()), 10)
            for name, atoms, bonds in tqdm(results, total=self._hmmit, desc="Indentify isomers", unit="molecule"):
                mname.append(name)
                print(name, ",".join([str(x) for x in atoms]), ";".join(
                    [",".join([str(y) for y in x]) for x in bonds]), file=fm)
                semaphore.release()
        self._mname = mname
    
    def _calmoleculeSMILESname(self, item):
        line, _ = item
        s = self._decompress(line).split()
        atoms = [int(x) for x in s[0].split(",")]
        bonds = [tuple(int(y) for y in x.split(","))
                 for x in s[1].split(";")] if len(s) == 3 else []
        name = self._convertSMILES(atoms, bonds)
        return name, atoms, bonds

    def _convertSMILES(self, atoms, bonds):
        m = Chem.RWMol(Chem.MolFromSmiles(''))
        d = {}
        for atomnumber in atoms:
            d[atomnumber] = m.AddAtom(
                Chem.Atom(self.atomname[self._atomtype[atomnumber-1]-1]))
        for atom1, atom2, level in bonds:
            m.AddBond(d[atom1], d[atom2], Chem.BondType.DOUBLE if level == 2 else (
                Chem.BondType.TRIPLE if level == 3 else (Chem.BondType.AROMATIC if level == 12 else Chem.BondType.SINGLE)))
        name = Chem.MolToSmiles(m)
        return name