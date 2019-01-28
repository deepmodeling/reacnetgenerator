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
        self.inputfilename = rng.inputfilename
        self.atomname = rng.atomname
        self.stepinterval = rng.stepinterval
        self.nproc = rng.nproc
        self._produce = rng.produce
        self._compress = rng.compress
        self._decompress = rng.decompress

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

    def _mo(self, i, bond, level, molecule, done, bondlist):
        """Connect molecule with Depth-First Search."""
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
            semaphore = Semaphore(self.nproc*150)
            results = pool.imap_unordered(
                self._readstepfunc, self._produce(
                    semaphore,
                    enumerate(
                        itertools.islice(
                            itertools.zip_longest(*[f] * _steplinenum),
                            0, None, self.stepinterval)),
                    None),
                100)
            for molecules, (step, thetimestep) in tqdm(results,
                                                       desc="Read bond information and Detect molecules",
                                                       unit="timestep"):
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
    def _readstepfunc(self, item):
        pass

    def _connectmolecule(self, bond, level):
        molecules = []
        done = np.zeros(self._N, dtype=bool)
        for i in range(1, self._N+1):
            if not done[i-1]:
                mole, done, bondlist = self._mo(i, bond, level, [], done, [])
                moleculestr = ' '.join(
                    (",".join((str(x) for x in sorted(mole))),
                     ";".join(
                         (",".join([str(y) for y in x])
                          for x in sorted(bondlist)))))
                molecules.append(self._compress(moleculestr))
        return molecules

    def _writemoleculetempfile(self, d):
        with tempfile.NamedTemporaryFile('wb', delete=False) as f:
            self.moleculetempfilename = f.name
            for key, value in d.items():
                f.write(
                    self._compress(
                        ' '.join(
                            (self._decompress(key),
                             ",".join((str(x) for x in value))))))
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
            elif line.startswith("ITEM: ATOMS"):
                return cls.ATOMS
            elif line.startswith("ITEM: NUMBER OF ATOMS"):
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
                    linecontent = self.LineType.linecontent(line)
                else:
                    if linecontent == self.LineType.ATOMS:
                        s = line.split()
                        step_atoms.append(
                            (int(s[0]),
                             Atom(
                                 self.atomname[int(s[1]) - 1],
                                 [float(x) for x in s[2: 5]])))
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
                        bond[int(s[1])-1].append(int(s[2]))
                        bond[int(s[2])-1].append(int(s[1]))
                        level = 12 if s[3] == 'ar' else int(s[3])
                        bondlevel[int(s[1])-1].append(level)
                        bondlevel[int(s[2])-1].append(level)
        return bond, bondlevel
