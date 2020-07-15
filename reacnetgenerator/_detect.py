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

import tempfile
import fileinput
import operator
from abc import ABCMeta, abstractmethod
from collections import defaultdict
from enum import Enum, auto

import numpy as np
try:
    from openbabel import openbabel
except ImportError:  # pragma: no cover
    import openbabel
from scipy.spatial import cKDTree
from ase import Atom, Atoms

from .dps import dps
from .utils import WriteBuffer, listtobytes, run_mp, SharedRNGData


class _Detect(SharedRNGData, metaclass=ABCMeta):
    """Detect molecules."""
    subclasses = {}

    def __init__(self, rng):
        SharedRNGData.__init__(self, rng, ['inputfilename', 'atomname', 'stepinterval', 'nproc', 'pbc',
                                           'cell'],
                               ['N', 'atomtype', 'step', 'timestep', 'temp1it', 'moleculetempfilename'])

    @classmethod
    def gettype(cls, rng):
        """Get the class for the input file type.

        Now ReacNetGen support the following files:
        LAMMPS bond files: http://lammps.sandia.gov/doc/fix_reax_bonds.html
        LAMMPS dump files: https://lammps.sandia.gov/doc/dump.html
        """
        if rng.inputfiletype not in cls.subclasses:
            raise ValueError(
                f"Unsupported input file type {rng.inputfiletype}")
        return cls.subclasses[rng.inputfiletype](rng)

    @classmethod
    def register_subclass(cls, message_type):
        def decorator(subclass):
            cls.subclasses[message_type] = subclass
            return subclass
        return decorator

    def detect(self):
        """Detect molecules."""
        self._readinputfile()
        self.returnkeys()

    def _readinputfile(self):
        d = defaultdict(list)
        timestep = {}
        with fileinput.input(files=self.inputfilename) as f:
            _steplinenum = self._readNfunc(f)
        with fileinput.input(files=self.inputfilename) as f:
            results = run_mp(self.nproc, func=self._readstepfunc, l=f, nlines=_steplinenum, return_num=True, interval=self.stepinterval,
                             desc="Read bond information and Detect molecules", unit="timestep")
            for molecules, (step, thetimestep) in results:
                for molecule in molecules:
                    d[molecule].append(step)
                timestep[step] = thetimestep
        self.temp1it = len(d)
        values_c = list(run_mp(self.nproc, func=self._compressvalue, l=d.values(
        ), unordered=False, desc="Save molecules", unit="molecule", total=self.temp1it))
        self._writemoleculetempfile((d.keys(), values_c))
        self.timestep = timestep
        self.step = len(timestep)

    def _compressvalue(self, x):
        return listtobytes(np.array(x))

    @abstractmethod
    def _readNfunc(self, f):
        pass

    @abstractmethod
    def _readstepfunc(self, item):
        pass

    def _connectmolecule(self, bond, level):
        return list([b''.join((listtobytes(mol),
                               listtobytes(bondlist))) for mol, bondlist in zip(*dps(bond, level))])

    def _writemoleculetempfile(self, d):
        with WriteBuffer(tempfile.NamedTemporaryFile('wb', delete=False)) as f:
            self.moleculetempfilename = f.name
            for mol in zip(*d):
                f.extend(mol)


@_Detect.register_subclass("lammpsbondfile")
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
        self.N = N
        self.atomtype = atomtype
        return steplinenum

    def _readstepfunc(self, item):
        step, lines = item
        bond = [None]*self.N
        level = [None]*self.N
        for line in lines:
            if line:
                if line[0] == "#":
                    if line.startswith("# Timestep"):
                        timestep = int(line.split()[-1])
                else:
                    s = line.split()
                    s0 = int(s[0])-1
                    s2 = int(s[2])
                    bond[s0] = map(self._get_idx, s[3:3+s2])
                    level[s0] = map(self._get_bo, s[4+s2:4+2*s2])
        molecules = self._connectmolecule(bond, level)
        return molecules, (step, timestep)

    @staticmethod
    def _get_idx(x):
        return int(x) - 1

    @staticmethod
    def _get_bo(x):
        return max(1, round(float(x)))


class _DetectCrd(_Detect):
    def _getbondfromcrd(self, step_atoms, cell):
        atomnumber = len(step_atoms)
        if self.pbc:
            # Apply period boundry conditions
            step_atoms.set_pbc(True)
            step_atoms.set_cell(cell)
            # add ghost atoms
            repeated_atoms = step_atoms.repeat(2)[atomnumber:]
            tree = cKDTree(step_atoms.get_positions())
            d = tree.query(repeated_atoms.get_positions(), k=1)[0]
            nearest = d < 5
            ghost_atoms = repeated_atoms[nearest]
            realnumber = np.where(nearest)[0] % atomnumber
            step_atoms += ghost_atoms
        # Use openbabel to connect atoms
        mol = openbabel.OBMol()
        mol.BeginModify()
        for idx, (num, position) in enumerate(zip(step_atoms.get_atomic_numbers(), step_atoms.positions)):
            a = mol.NewAtom(idx)
            a.SetAtomicNum(int(num))
            a.SetVector(*position)
        mol.ConnectTheDots()
        mol.PerceiveBondOrders()
        mol.EndModify()
        bond = [[] for i in range(atomnumber)]
        bondlevel = [[] for i in range(atomnumber)]
        for b in openbabel.OBMolBondIter(mol):
            s1 = b.GetBeginAtom().GetId()
            s2 = b.GetEndAtom().GetId()
            if s1 >= atomnumber and s2 >= atomnumber:
                # duplicated
                continue
            elif s1 >= atomnumber:
                s1 = realnumber[s1-atomnumber]
            elif s2 >= atomnumber:
                s2 = realnumber[s2-atomnumber]
            level = b.GetBondOrder()
            if level == 5:
                # aromatic, 5 in openbabel but 12 in rdkit
                level = 12
            bond[s1].append(s2)
            bond[s2].append(s1)
            bondlevel[s1].append(level)
            bondlevel[s2].append(level)
        return bond, bondlevel


@_Detect.register_subclass("lammpsdumpfile")
class _DetectLAMMPSdump(_DetectCrd):
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
        self.N = N
        self.atomtype = atomtype
        return steplinenum

    def _readstepfunc(self, item):
        step, lines = item
        step_atoms = []
        ss = []
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
                        ss.append(list(map(float, s)))
        ss = np.array(ss)
        if ss.shape[1]>2:
            xy = ss[0][2]
            xz = ss[1][2]
            yz = ss[2][2]
        else:
            xy, xz, yz = 0., 0., 0.
        xlo = ss[0][0] - min(0., xy, xz, xy+xz)
        xhi = ss[0][1] - max(0., xy, xz, xy+xz)
        ylo = ss[1][0] - min(0., yz)
        yhi = ss[1][1] - max(0., yz)
        zlo = ss[2][0]
        zhi = ss[2][1] 
        boxsize = np.array([[xhi-xlo, 0., 0.],
                   [xy, yhi-ylo, 0.],
                   [xz, yz, zhi-zlo]])
        _, step_atoms = zip(*sorted(step_atoms, key=operator.itemgetter(0)))
        step_atoms = Atoms(step_atoms)
        bond, level = self._getbondfromcrd(step_atoms, boxsize)
        molecules = self._connectmolecule(bond, level)
        return molecules, timestep


@_Detect.register_subclass("xyz")
class _Detectxyz(_DetectCrd):
    """xyz file. Two frames are connected. Cell information must be inputed by user."""

    def _readNfunc(self, f):
        atomname_dict = dict(zip(self.atomname.tolist(), range(self.atomname.size)))
        for index, line in enumerate(f):
            s = line.split()
            if index == 0:
                N = int(line.strip())
                atomtype = np.zeros(N, dtype=int)
            elif index > N+1:
                break
            elif index > 1:
                atomtype[index-2] = atomname_dict[s[0]]
        self.N = N
        self.atomtype = atomtype
        steplinenum = N + 2
        return steplinenum

    def _readstepfunc(self, item):
        step, lines = item
        timestep = step, step
        step_atoms = []
        boxsize = self.cell
        if self.pbc and boxsize is None:
            raise RuntimeError("No cell information is given.")
        for index, line in enumerate(lines):
            s = line.split()
            if index > 1:
                step_atoms.append((index-1, Atom(s[0] ,tuple((float(x) for x in s[1:4])))))
        _, step_atoms = zip(*sorted(step_atoms, key=operator.itemgetter(0)))
        step_atoms = Atoms(step_atoms)
        bond, level = self._getbondfromcrd(step_atoms, boxsize)
        molecules = self._connectmolecule(bond, level)
        return molecules, timestep
