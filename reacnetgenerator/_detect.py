"""Detect molecules.

There are two types of input files that could be imported by ReacNetGen,
the first part of necessary is the trajectory from reactive MD, the second
part can be the bond information normally given by simulation using ReaxFF.
In fact, atomic coordinates can be converted to the bond information with
the Open Babel software. As a results, ReacNetGen can both processes ReaxFF
trajectories, AIMD trajectories, and other kinds of reactive trajectories.
With the bond information, molecules can be detected from atoms by Depth-first
search at every timestep. By using this way, all molecules in a given
trajectory can be acquired. Molecules consisting of same atoms and bonds can
be considered as the same molecule.

Reference:
[1] O’Boyle, N. M.; Banck, M.; James, C. A.; Morley, C.; Vandermeersch, T.;
Hutchison, G. Open Babel: An open chemical toolbox. J. Cheminf. 2011, 3(1),
33-47.
[2] Tarjan, R. Depth-first search and linear graph algorithms. SIAM J. Comput.
1972, 1 (2), 146-160.
"""

import itertools
import tempfile
from abc import ABCMeta, abstractmethod
from collections import defaultdict
from enum import Enum, auto
from multiprocessing import Pool, Semaphore

import numpy as np
import openbabel
from ase import Atom, Atoms
from tqdm import tqdm

from .dps import dps


class InputFileType(Enum):
    """Enum for input file types.

    Now ReacNetGen support the following files:
    LAMMPS bond files: http://lammps.sandia.gov/doc/fix_reax_bonds.html
    LAMMPS dump files: https://lammps.sandia.gov/doc/dump.html
    """

    LAMMPSBOND = auto()
    LAMMPSDUMP = auto()


class _Detect(metaclass=ABCMeta):
    """Detect molecules."""

    def __init__(self, rng):
        self.rng = rng
        self.inputfilenames = rng.inputfilenames
        self.atomname = rng.atomname
        self.stepinterval = rng.stepinterval
        self.nproc = rng.nproc
        self._produce = rng.produce
        self._compress = rng.compress
        self._listtobytes = rng.listtobytes

        self._N = None
        self._atomtype = None
        self._step = None
        self.moleculetempfilename = None
        self._temp1it = None
        self._timestep = None

    @staticmethod
    def gettype(inputtype):
        """Get the class for the input file type."""
        if inputtype == InputFileType.LAMMPSBOND:
            detectclass = _DetectLAMMPSbond
        elif inputtype == InputFileType.LAMMPSDUMP:
            detectclass = _DetectLAMMPSdump
        else:
            raise RuntimeError("Wrong input file type")
        return detectclass

    def detect(self):
        """Detect molecules."""
        self._readinputfile()
        self.rng.N = self._N
        self.rng.atomtype = self._atomtype
        self.rng.step = self._step
        self.rng.timestep = self._timestep
        self.rng.moleculetempfilename = self.moleculetempfilename
        self.rng.temp1it = self._temp1it

    def _readinputfile(self):
        d = defaultdict(list)
        timestep = {}
        with Pool(self.nproc, maxtasksperchild=1000) as pool:
            semaphore = Semaphore(self.nproc*150)
            for i, inputfilename in enumerate(self.inputfilenames):
                with open(inputfilename) as f:
                    if i==0:
                        _steplinenum = self._readNfunc(f)
                        f.seek(0)
                    results = pool.imap_unordered(
                        self._readstepfunc, self._produce(
                            semaphore,
                            enumerate(
                                itertools.islice(
                                    itertools.zip_longest(*[f] * _steplinenum),
                                    0, None, self.stepinterval), len(d)),
                            None),
                        100)
                    for molecules, (step, thetimestep) in tqdm(results,
                                                            desc="Read bond information and Detect molecules",
                                                            unit="timestep"):
                        for molecule in molecules:
                            d[molecule].append(step)
                        timestep[step] = thetimestep
                        semaphore.release()
            self._temp1it = len(d)
            values_c = list(tqdm(pool.imap(self._compressvalue,
                                                    d.values(),
                                                    100), desc="Save molecules", unit="molecule", total=self._temp1it))
        pool.close()
        self._writemoleculetempfile((d.keys(), values_c))
        self._timestep = timestep
        self._step = len(timestep)
        pool.join()

    def _compressvalue(self, x):
        return self._listtobytes(np.array(x))

    @abstractmethod
    def _readNfunc(self, f):
        pass

    @abstractmethod
    def _readstepfunc(self, item):
        pass

    def _connectmolecule(self, bond, level):
        return list([b' '.join((self._listtobytes(mol),
                                self._listtobytes(bondlist))) for mol, bondlist in zip(*dps(bond, level))])

    def _writemoleculetempfile(self, d):
        buff = []
        with tempfile.NamedTemporaryFile('wb', delete=False) as f:
            self.moleculetempfilename = f.name
            for mol in zip(*d):
                buff.extend(mol)
                if len(buff) > 30*self.nproc:
                    f.write(b''.join(buff))
                    buff[:] = []
            if buff:
                f.write(b''.join(buff))


class _DetectLAMMPSbond(_Detect):
    def _readNfunc(self, f):
        iscompleted = False
        for index, line in enumerate(f):
            if line[0] == '#':
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
                atomtype[int(s[0])-1] = int(s[1])-1
        steplinenum = stepbindex-stepaindex
        self._N = N
        self._atomtype = atomtype
        return steplinenum

    def _readstepfunc(self, item):
        (step, lines), _ = item
        bond = [None]*self._N
        level = [None]*self._N
        for line in lines:
            if line:
                if line[0] == "#":
                    if line.startswith("# Timestep"):
                        timestep = int(line.split()[-1])
                else:
                    s = line.split()
                    s0 = int(s[0])-1
                    s2 = int(s[2])
                    bond[s0] = map(lambda x: int(x)-1, s[3:3+s2])
                    level[s0] = map(lambda x: max(
                        1, round(float(x))), s[4+s2:4+2*s2])
        molecules = self._connectmolecule(bond, level)
        return molecules, (step, timestep)


class _DetectLAMMPSdump(_Detect):
    class LineType(Enum):
        """Line type in the LAMMPS dump files."""

        TIMESTEP = auto()
        ATOMS = auto()
        NUMBER = auto()
        OTHER = auto()

        @classmethod
        def linecontent(cls, line):
            """Return line content."""
            if line.startswith("ITEM: TIMESTEP"):
                return cls.TIMESTEP
            if line.startswith("ITEM: ATOMS"):
                return cls.ATOMS
            if line.startswith("ITEM: NUMBER OF ATOMS"):
                return cls.NUMBER
            return cls.OTHER

    def _readNfunc(self, f):
        iscompleted = False
        for index, line in enumerate(f):
            if line.startswith("ITEM:"):
                linecontent = self.LineType.linecontent(line)
            else:
                if linecontent == self.LineType.NUMBER:
                    if iscompleted:
                        stepbindex = index
                        break
                    else:
                        iscompleted = True
                        stepaindex = index
                    N = int(line.split()[0])
                    atomtype = np.zeros(N, dtype=int)
                elif linecontent == self.LineType.ATOMS:
                    s = line.split()
                    atomtype[int(s[0])-1] = int(s[1])-1
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
                    linecontent = self.LineType.linecontent(line)
                else:
                    if linecontent == self.LineType.ATOMS:
                        s = line.split()
                        step_atoms.append(
                            (int(s[0]),
                             Atom(
                                 self.atomname[int(s[1]) - 1],
                                 tuple(map(float, s[2: 5])))))
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
        xyzstring = ''.join((f"{atomnumber}\nReacNetGenerator\n", "\n".join(
            [f'{s:2s} {x:22.15f} {y:22.15f} {z:22.15f}'
             for s, (x, y, z)
             in
             zip(
                 step_atoms.get_chemical_symbols(),
                 step_atoms.positions)])))
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
                        bond[int(s[1])-1].append(int(s[2])-1)
                        bond[int(s[2])-1].append(int(s[1])-1)
                        level = 12 if s[3] == 'ar' else int(s[3])
                        bondlevel[int(s[1])-1].append(level)
                        bondlevel[int(s[2])-1].append(level)
        return bond, bondlevel
