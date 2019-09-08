# cython: language_level=3
# cython: linetrace=True
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
import fileinput
from abc import ABCMeta, abstractmethod
from collections import defaultdict
from enum import Enum, auto
from multiprocessing import Pool, Semaphore

import numpy as np
import openbabel
from scipy.spatial import cKDTree
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
        self.pbc = rng.pbc

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
        pool = Pool(self.nproc, maxtasksperchild=1000)
        try:
            semaphore = Semaphore(self.nproc*150)
            with fileinput.input(files=self.inputfilenames) as f:
                _steplinenum = self._readNfunc(f)
            with fileinput.input(files=self.inputfilenames) as f:
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
        finally:
            pool.close()
            pool.join()
        self._writemoleculetempfile((d.keys(), values_c))
        self._timestep = timestep
        self._step = len(timestep)

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
                        steplinenum = index-stepaindex
                        break
                    else:
                        iscompleted = True
                        stepaindex = index
                    N = [int(s) for s in line.split() if s.isdigit()][0]
                    atomtype = np.zeros(N, dtype=np.int)
            else:
                s = line.split()
                atomtype[int(s[0])-1] = int(s[1])-1
        else:
            steplinenum = index + 1
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
        BOX = auto()
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
            if line.startswith("ITEM: BOX"):
                return cls.BOX
            return cls.OTHER

    def _readNfunc(self, f):
        iscompleted = False
        for index, line in enumerate(f):
            if line.startswith("ITEM:"):
                linecontent = self.LineType.linecontent(line)
                if linecontent == self.LineType.ATOMS:
                    keys = line.split()
                    self.id_idx = keys.index('id') - 2
                    self.tidx = keys.index('type') - 2
                    self.xidx = keys.index('x') - 2
                    self.yidx = keys.index('y') - 2
                    self.zidx = keys.index('z') - 2
            else:
                if linecontent == self.LineType.NUMBER:
                    if iscompleted:
                        steplinenum = index-stepaindex
                        break
                    else:
                        iscompleted = True
                        stepaindex = index
                    N = int(line.split()[0])
                    atomtype = np.zeros(N, dtype=int)
                elif linecontent == self.LineType.ATOMS:
                    s = line.split()
                    atomtype[int(s[0])-1] = int(s[1])-1
        else:
            steplinenum = index + 1
        self._N = N
        self._atomtype = atomtype
        return steplinenum

    def _readstepfunc(self, item):
        (step, lines), _ = item
        step_atoms = []
        boxsize = []
        for line in lines:
            if line:
                if line.startswith("ITEM:"):
                    linecontent = self.LineType.linecontent(line)
                else:
                    if linecontent == self.LineType.ATOMS:
                        s = line.split()
                        step_atoms.append(
                            (int(s[self.id_idx]),
                             Atom(
                                 self.atomname[int(s[self.tidx]) - 1],
                                 (float(s[self.xidx]), float(s[self.yidx]), float(s[self.zidx])))))
                    elif linecontent == self.LineType.TIMESTEP:
                        timestep = step, int(line.split()[0])
                    elif linecontent == self.LineType.BOX:
                        s = line.split()
                        boxsize.append(float(s[1])-float(s[0]))
        _, step_atoms = zip(*sorted(step_atoms, key=lambda a: a[0]))
        step_atoms = Atoms(step_atoms)
        bond, level = self._getbondfromcrd(step_atoms, boxsize)
        molecules = self._connectmolecule(bond, level)
        return molecules, timestep

    def _getbondfromcrd(self, step_atoms, cell):
        atomnumber = len(step_atoms)
        ghosts = {}
        if self.pbc:
            # Apply period boundry conditions
            step_atoms.set_pbc(True)
            step_atoms.set_cell(cell)
            # add ghost atoms
            repeated_atoms = step_atoms.repeat(2)[atomnumber:]
            tree = cKDTree(step_atoms.get_positions())
            d = tree.query(repeated_atoms.get_positions(), k=1)[0]
            nearest = d<5
            ghost_atoms = repeated_atoms[nearest]
            realnumber = np.where(nearest)[0] % atomnumber
            step_atoms += ghost_atoms
        xyzstring = ''.join((f"{len(step_atoms)}\nReacNetGenerator\n", "\n".join(
            [f'{s:2s} {x:22.15f} {y:22.15f} {z:22.15f}'
             for s, (x, y, z) in zip(
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
                        s1 = int(s[1])-1
                        s2 = int(s[2])-1
                        if s1 >= atomnumber and s2 >= atomnumber:
                            # duplicated
                            continue
                        elif s1 >= atomnumber:
                            s1 = realnumber[s1-atomnumber]
                        elif s2 >= atomnumber:
                            s2 = realnumber[s2-atomnumber]
                        bond[s1].append(s2)
                        bond[s2].append(s1)
                        level = 12 if s[3] == 'ar' else (1 if s[3] == 'am' else int(s[3]))
                        bondlevel[s1].append(level)
                        bondlevel[s2].append(level)
        return bond, bondlevel
