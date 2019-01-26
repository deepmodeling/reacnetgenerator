''' Detect molecules '''

import itertools
from collections import defaultdict
from multiprocessing import Pool, Semaphore
from enum import Enum, auto
import tempfile
from abc import ABCMeta, abstractmethod

import numpy as np
from ase import Atom, Atoms
from tqdm import tqdm
import openbabel


class InputFileType(Enum):
    LAMMPSBOND = auto()
    LAMMPSDUMP = auto()


class _Detect(metaclass=ABCMeta):
    '''Detect molecules'''

    def __init__(self, rng):
        self.rng = rng
        self.inputfilename = rng.inputfilename
        self.atomname = rng.atomname
        self.stepinterval = rng.stepinterval
        self.nproc = rng.nproc
        self._produce = rng._produce
        self._compress = rng._compress
        self._decompress = rng._decompress

        self._N = None
        self._atomtype = None
        self._step = None
        self.moleculetempfilename = None
        self._temp1it = None
        self._timestep = None

    @staticmethod
    def gettype(inputtype):
        if inputtype == InputFileType.LAMMPSBOND:
            detectclass = _DetectLAMMPSbond
        elif inputtype == InputFileType.LAMMPSDUMP:
            detectclass = _DetectLAMMPSdump
        else:
            raise RuntimeError("Wrong input file type")
        return detectclass

    def detect(self):
        ''' Detect molecules'''
        self._readinputfile()
        self.rng._N = self._N
        self.rng._atomtype = self._atomtype
        self.rng._step = self._step
        self.rng._timestep = self._timestep
        self.rng.moleculetempfilename = self.moleculetempfilename
        self.rng._temp1it = self._temp1it

    def _mo(self, i, bond, level, molecule, done, bondlist):
        ''' connect molecule '''
        molecule.append(i)
        done[i-1] = True
        for b, l in zip(bond[i-1], level[i-1]):
            if not done[b-1]:
                bondlist.append((i, b, l) if i < b else (b, i, l))
                molecule, done, bondlist = self._mo(
                    b, bond, level, molecule, done, bondlist)
        return molecule, done, bondlist

    def _readinputfile(self):
        d = defaultdict(list)
        timestep = {}
        with open(self.inputfilename) as f, Pool(self.nproc, maxtasksperchild=1000) as pool:
            _steplinenum = self._readNfunc(f)
            f.seek(0)
            semaphore = Semaphore(self.nproc*15)
            results = pool.imap_unordered(self._readstepfunc, self._produce(semaphore, enumerate(itertools.islice(
                itertools.zip_longest(*[f]*_steplinenum), 0, None, self.stepinterval)), None), 10)
            for molecules, (step, thetimestep) in tqdm(results, desc="Read bond information and Detect molecules", unit="timestep"):
                for molecule in molecules:
                    d[molecule].append(step)
                timestep[step] = thetimestep
                semaphore.release()
        self._writemoleculetempfile(d)
        self._timestep = timestep
        self._step = len(timestep)-1

    @abstractmethod
    def _readNfunc(self, f):
        pass

    @abstractmethod
    def _readstepfunc(self):
        pass

    def _connectmolecule(self, bond, level):
        molecules = []
        done = np.zeros(self._N, dtype=bool)
        for i in range(1, self._N+1):
            if not done[i-1]:
                mole, done, bondlist = self._mo(i, bond, level, [], done, [])
                moleculestr = ' '.join((",".join((str(x) for x in sorted(mole))), ";".join(
                    (",".join([str(y) for y in x]) for x in sorted(bondlist)))))
                molecules.append(self._compress(moleculestr))
        return molecules

    def _writemoleculetempfile(self, d):
        with tempfile.NamedTemporaryFile('wb', delete=False) as f:
            self.moleculetempfilename = f.name
            for key, value in d.items():
                f.write(self._compress(
                    ' '.join((self._decompress(key), ",".join((str(x) for x in value))))))
        self._temp1it = len(d)


class _DetectLAMMPSbond(_Detect):
    def _readNfunc(self, f):
        iscompleted = False
        for index, line in enumerate(f):
            if line.startswith("#"):
                if line.startswith("# Number of particles"):
                    if iscompleted:
                        stepbindex = index
                        break
                    else:
                        iscompleted = True
                        stepaindex = index
                    N = [int(s) for s in line.split() if s.isdigit()][0]
                    atomtype = np.zeros(N, dtype=np.int)
            else:
                s = line.split()
                atomtype[int(s[0])-1] = int(s[1])
        steplinenum = stepbindex-stepaindex
        self._N = N
        self._atomtype = atomtype
        return steplinenum

    def _readstepfunc(self, item):
        (step, lines), _ = item
        bond = [None for x in range(self._N)]
        level = [None for x in range(self._N)]
        for line in lines:
            if line:
                if line.startswith("#"):
                    if line.startswith("# Timestep"):
                        timestep = int(line.split()[-1])
                else:
                    s = line.split()
                    bond[int(s[0])-1] = [int(x) for x in s[3:3+int(s[2])]]
                    level[int(s[0])-1] = [max(1, round(float(x)))
                                          for x in s[4+int(s[2]):4+2*int(s[2])]]
        molecules = self._connectmolecule(bond, level)
        return molecules, (step, timestep)


class _DetectLAMMPSdump(_Detect):
    class LineType(Enum):
        TIMESTEP = auto()
        ATOMS = auto()
        NUMBER = auto()
        OTHER = auto()

    def _readNfunc(self, f):
        iscompleted = False
        for index, line in enumerate(f):
            if line.startswith("ITEM:"):
                linecontent = self.LineType.TIMESTEP if line.startswith("ITEM: TIMESTEP") else (self.LineType.ATOMS if line.startswith(
                    "ITEM: ATOMS") else (self.LineType.NUMBER if line.startswith("ITEM: NUMBER OF ATOMS") else self.LineType.OTHER))
            else:
                if linecontent == self.LineType.NUMBER:
                    if iscompleted:
                        stepbindex = index
                        break
                    else:
                        iscompleted = True
                        stepaindex = index
                    N = int(line.split()[0])
                    atomtype = np.zeros(N, dtype=np.int)
                elif linecontent == self.LineType.ATOMS:
                    s = line.split()
                    atomtype[int(s[0])-1] = int(s[1])
        steplinenum = stepbindex-stepaindex
        self._N = N
        self._atomtype = atomtype
        return steplinenum

    def _readstepfunc(self, item):
        (step, lines), _ = item
        step_atoms = []
        for line in lines:
            if line:
                if line.startswith("ITEM:"):
                    linecontent = self.LineType.TIMESTEP if line.startswith("ITEM: TIMESTEP") else (self.LineType.ATOMS if line.startswith(
                        "ITEM: ATOMS") else (self.LineType.NUMBER if line.startswith("ITEM: NUMBER OF ATOMS") else self.LineType.OTHER))
                else:
                    if linecontent == self.LineType.ATOMS:
                        s = line.split()
                        step_atoms.append(
                            (int(s[0]), Atom(self.atomname[int(s[1])-1], [float(x) for x in s[2:5]])))
                    elif linecontent == self.LineType.TIMESTEP:
                        timestep = step, int(line.split()[0])
        _, step_atoms = zip(*sorted(step_atoms, key=lambda a: a[0]))
        step_atoms = Atoms(step_atoms)
        bond, level = self._getbondfromcrd(step_atoms)
        molecules = self._connectmolecule(bond, level)
        return molecules, timestep

    @classmethod
    def _getbondfromcrd(cls, step_atoms):
        atomnumber = len(step_atoms)
        xyzstring = f"{atomnumber}\nReacNetGenerator\n"+"\n".join(
            [f'{s:2s} {x:22.15f} {y:22.15f} {z:22.15f}' for s, (x, y, z) in zip(step_atoms.get_chemical_symbols(), step_atoms.positions)])
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats('xyz', 'mol2')
        mol = openbabel.OBMol()
        conv.ReadString(mol, xyzstring)
        mol2string = conv.WriteString(mol)
        linecontent = -1
        bond = [[] for i in range(atomnumber)]
        bondlevel = [[] for i in range(atomnumber)]
        for line in mol2string.split('\n'):
            if line.startswith("@<TRIPOS>BOND"):
                linecontent = 0
            else:
                if linecontent == 0:
                    s = line.split()
                    if len(s) > 3:
                        bond[int(s[1])-1].append(int(s[2]))
                        bond[int(s[2])-1].append(int(s[1]))
                        level = 12 if s[3] == 'ar' else int(s[3])
                        bondlevel[int(s[1])-1].append(level)
                        bondlevel[int(s[2])-1].append(level)
        return bond, bondlevel
