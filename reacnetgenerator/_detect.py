# SPDX-License-Identifier: LGPL-3.0-or-later
# cython: language_level=3
# cython: linetrace=True
"""Detect molecules.

There are two types of input files that could be imported by ReacNetGenerator,
the first part of necessary is the trajectory from reactive MD, the second
part can be the bond information normally given by simulation using ReaxFF.
In fact, atomic coordinates can be converted to the bond information with
the Open Babel software.[1]_ As a results, ReacNetGenerator can both processes ReaxFF
trajectories, AIMD trajectories, and other kinds of reactive trajectories.
With the bond information, molecules can be detected from atoms by Depth-first
search at every timestep.[2]_ By using this way, all molecules in a given
trajectory can be acquired. Molecules consisting of same atoms and bonds can
be considered as the same molecule.

References
----------
.. [1] O'Boyle, N. M.; Banck, M.; James, C. A.; Morley, C.; Vandermeersch, T.;
   Hutchison, G. Open Babel: An open chemical toolbox. J. Cheminf. 2011, 3(1),
   33-47.
.. [2] Tarjan, R. Depth-first search and linear graph algorithms. SIAM J. Comput.
   1972, 1 (2), 146-160.
"""

import fileinput
import operator
import tempfile
from abc import ABCMeta, abstractmethod
from collections import defaultdict
from enum import Enum, auto
from typing import ClassVar, List, Optional, Tuple

import numpy as np

try:
    from openbabel import __version__ as obversion
    from openbabel import openbabel
except ImportError:  # pragma: no cover
    raise ImportError("Open Babel 3.1.0 is required.")
from ase import Atom, Atoms
from packaging import version

from ._logging import logger
from .dps import dps  # type:ignore
from .utils import SharedRNGData, WriteBuffer, listtobytes, run_mp

# ASE imports
try:
    from ase.neighborlist import natural_cutoffs, neighbor_list

    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False

# Scipy imports for connected components
try:
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

if version.parse(obversion) < version.parse("3.1.0"):  # pragma: no cover
    raise ImportError("Open Babel 3.1.0 is required.")

openbabel.obErrorLog.StopLogging()


class _Detect(SharedRNGData, metaclass=ABCMeta):
    """Detect molecules.

    Parameters
    ----------
    rng: reacnetgenerator.ReacNetGenerator
        The ReacNetGenerator class.
    """

    subclasses: ClassVar[dict] = {}

    # type hints
    # TODO: we need a better way to communicate with parameters
    inputfilename: list
    atomname: np.ndarray
    stepinterval: int
    nproc: int
    pbc: bool
    cell: np.ndarray

    def __init__(self, rng):
        SharedRNGData.__init__(
            self,
            rng,
            ["inputfilename", "atomname", "stepinterval", "nproc", "pbc", "cell"],
            ["N", "atomtype", "step", "timestep", "temp1it", "moleculetempfilename"],
        )

    @classmethod
    def gettype(cls, rng):
        """Get the class for the input file type.

        Now ReacNetGen support the following files:
            - lammpsbondfile: LAMMPS bond files, see http://lammps.sandia.gov/doc/fix_reax_bonds.html
            - lammpsdumpfile: LAMMPS dump files, see https://lammps.sandia.gov/doc/dump.html
            - xyz: xyz files
            - extxyz: extended xyz files

        Parameters
        ----------
        rng : reacnetgenerator.ReacNetGenerator
            The ReacNetGenerator class.

        Returns
        -------
        cls.subclasses[rng.inputfiletype](rng): class
            The _Detect class for specific file type.
        """
        if rng.inputfiletype not in cls.subclasses:
            raise ValueError(f"Unsupported input file type {rng.inputfiletype}")
        return cls.subclasses[rng.inputfiletype](rng)

    @classmethod
    def register_subclass(cls, message_type):
        """Register the file type, used as a decorator.

        Parameters
        ----------
        message_type : str
            The type name to register, such as `xyz`.

        Returns
        -------
        decorator: function
            The decorator that used for a subclass.

        Examples
        --------
        @_Detect.register_subclass("lammpsbondfile")
        class _DetectLAMMPSbond(_Detect):
        """

        def decorator(subclass):
            cls.subclasses[message_type] = subclass
            return subclass

        return decorator

    def detect(self):
        """Detect molecules."""
        self._readinputfile()
        self.returnkeys()

    def _readinputfile(self):
        """Read the input file."""
        d = defaultdict(list)
        timestep = {}
        with fileinput.input(files=self.inputfilename) as f:
            _steplinenum = self._readNfunc(f)
        with fileinput.input(files=self.inputfilename) as f:
            results = run_mp(
                self.nproc,
                func=self._readstepfunc,
                l=f,
                nlines=_steplinenum,
                return_num=True,
                interval=self.stepinterval,
                desc="Read bond information and Detect molecules",
                unit="timestep",
            )
            for molecules, (step, thetimestep) in results:
                for molecule in molecules:
                    d[molecule].append(step)
                timestep[step] = thetimestep
        self.temp1it = len(d)
        values_c = list(
            run_mp(
                self.nproc,
                func=self._compressvalue,
                l=d.values(),
                unordered=False,
                desc="Save molecules",
                unit="molecule",
                total=self.temp1it,
            )
        )
        self._writemoleculetempfile((d.keys(), values_c))
        self.timestep = timestep
        self.step = len(timestep)

    def _compressvalue(self, x):
        return listtobytes(np.array(x))

    @abstractmethod
    def _readNfunc(self, f) -> int:
        pass

    @abstractmethod
    def _readstepfunc(self, item) -> Tuple[List[bytes], Tuple[int, int]]:
        pass

    def _connectmolecule(self, bond, level):
        # Use scipy for connected components if available (better for large periodic systems)
        if SCIPY_AVAILABLE:
            # Build adjacency list into a sparse matrix for connected_components
            n_atoms = len(bond)
            row_ind = []
            col_ind = []

            for i, neighbors in enumerate(bond):
                if neighbors:
                    for j in neighbors:
                        row_ind.append(i)
                        col_ind.append(j)

            # Create CSR matrix
            data = [1] * len(row_ind)
            graph = csr_matrix((data, (row_ind, col_ind)), shape=(n_atoms, n_atoms))

            # Find connected components
            _n_components, labels = connected_components(graph, directed=False)

            # Group atoms by component
            components = defaultdict(list)
            for atom_idx, comp_label in enumerate(labels):
                components[comp_label].append(atom_idx)

            # For each component, create the bond list
            results = []
            for comp_atoms in components.values():
                # Extract bonds within this component
                comp_bonds = []
                for atom_idx in comp_atoms:
                    if bond[atom_idx]:
                        for neighbor_idx in bond[atom_idx]:
                            # Only add bond if it's within the same component and in increasing order
                            if neighbor_idx in comp_atoms and atom_idx < neighbor_idx:
                                bond_idx = bond[atom_idx].index(neighbor_idx)
                                bond_level = (
                                    level[atom_idx][bond_idx] if level[atom_idx] else 1
                                )
                                comp_bonds.append((atom_idx, neighbor_idx, bond_level))

                # Create the required format
                results.append(
                    b"".join(
                        (
                            listtobytes(comp_atoms),
                            listtobytes([(int(b[0]), int(b[1])) for b in comp_bonds]),
                            listtobytes([int(b[2]) for b in comp_bonds]),
                        )
                    )
                )
            return results
        else:
            # Fallback to the existing DPS implementation
            mols, bondlists = dps(bond, level)
            return [
                b"".join(
                    (
                        listtobytes(mol),
                        listtobytes([(int(i[0]), int(i[1])) for i in bondlist]),
                        listtobytes([int(i[2]) for i in bondlist]),
                    )
                )
                for mol, bondlist in zip(mols, bondlists)
            ]

    def _writemoleculetempfile(self, d):
        with WriteBuffer(tempfile.NamedTemporaryFile("wb", delete=False)) as f:
            self.moleculetempfilename = f.name
            for mol in zip(*d):
                f.extend(mol)


@_Detect.register_subclass("bond")
@_Detect.register_subclass("lammpsbondfile")
class _DetectLAMMPSbond(_Detect):
    def _readNfunc(self, f) -> int:
        iscompleted = False
        N = None
        atomtype = None
        stepaindex = None
        index = -1
        for index, line in enumerate(f):
            if line[0] == "#":
                if line.startswith("# Number of particles"):
                    if iscompleted:
                        assert stepaindex is not None
                        steplinenum = index - stepaindex
                        break
                    else:
                        iscompleted = True
                        stepaindex = index
                    N = next(int(s) for s in line.split() if s.isdigit())
                    atomtype = np.zeros(N, dtype=int)
            else:
                s = line.split()
                assert atomtype is not None
                atomtype[int(s[0]) - 1] = int(s[1]) - 1
        else:
            steplinenum = index + 1
        assert N is not None and atomtype is not None
        self.N = N
        self.atomtype = atomtype
        return steplinenum

    def _readstepfunc(self, item) -> Tuple[List[bytes], Tuple[int, int]]:
        step, lines = item
        bond: list[Optional[Tuple[int, ...]]] = [None] * self.N
        level: list[Optional[Tuple[int, ...]]] = [None] * self.N
        timestep = None
        for line in lines:
            if line:
                if line[0] == "#":
                    if line.startswith("# Timestep"):
                        timestep = int(line.split()[-1])
                else:
                    s = line.split()
                    s0 = int(s[0]) - 1
                    s2 = int(s[2])
                    bond[s0] = tuple(map(self._get_idx, s[3 : 3 + s2]))
                    level[s0] = tuple(map(self._get_bo, s[4 + s2 : 4 + 2 * s2]))
        molecules = self._connectmolecule(bond, level)
        assert timestep is not None
        return molecules, (step, timestep)

    @staticmethod
    def _get_idx(x):
        return int(x) - 1

    @staticmethod
    def _get_bo(x):
        return max(1, round(float(x)))


class _DetectCrd(_Detect):
    def _parse_custom_cutoffs(self, cutoff_str):
        """Parse the custom cutoff string into a dictionary.

        Parameters
        ----------
        cutoff_str : str
            Custom cutoffs in the format "El1-El2:dist,El3-El4:dist"

        Returns
        -------
        dict
            Dictionary with frozenset pairs as keys and cutoffs as values
        """
        if not cutoff_str:
            return {}

        result = {}
        pairs = cutoff_str.split(",")
        for pair in pairs:
            pair = pair.strip()  # Clean input to allow spaces
            if not pair:  # Skip empty segments
                continue

            if ":" not in pair:
                raise ValueError(
                    f"Invalid custom cutoff format '{pair}'. Expected 'Element1-Element2:distance'. "
                    f"Example: 'Al-O:2.5,C-H:1.1'"
                )

            atom_pair_str, dist_str = pair.split(":", 1)
            atom_pair_str = atom_pair_str.strip()
            dist_str = dist_str.strip()

            if "-" not in atom_pair_str:
                raise ValueError(
                    f"Invalid custom cutoff format '{pair}'. Expected 'Element1-Element2:distance'. "
                    f"Example: 'Al-O:2.5,C-H:1.1'"
                )

            atom1, atom2 = atom_pair_str.split("-", 1)
            atom1 = atom1.strip()
            atom2 = atom2.strip()

            try:
                distance = float(dist_str)
            except ValueError:
                raise ValueError(
                    f"Invalid distance value '{dist_str}' in '{pair}'. Expected a number. "
                    f"Example: 'Al-O:2.5,C-H:1.1'"
                )

            result[frozenset({atom1, atom2})] = distance
        return result

    def _getbondfromase(self, step_atoms: Atoms, cell: np.ndarray):
        """Perceive bonds using ASE neighbor list with custom cutoffs.

        Parameters
        ----------
        step_atoms : ase.Atoms
            Atoms in a step.
        cell : np.ndarray
            Cell in the shape (3, 3).

        Returns
        -------
        list[list[int]]
            Connected atoms for each atom.
        list[list[int]]
            Bond levels (all set to 1 since ASE doesn't detect bond order).
        """
        if not ASE_AVAILABLE:
            raise ImportError("ASE is not available. Please install ase package.")

        # Step 1: Preparation
        # Get raw covalent radii
        radii = natural_cutoffs(step_atoms, mult=1.0)
        symbols = step_atoms.get_chemical_symbols()

        custom_cutoffs = {}
        if self.rng.custom_cutoffs:
            custom_cutoffs = self._parse_custom_cutoffs(self.rng.custom_cutoffs)

        step_atoms.set_cell(cell)
        step_atoms.set_pbc(self.pbc)  # Ensure neighbor_list respects PBC
        if self.pbc:
            step_atoms.wrap()

        # Step 2: Calculate Global Search Cutoff (Absolute Angstroms)
        # Determine max_radius = max(radii) (handle empty case)
        if len(radii) > 0:
            max_radius = max(radii)
        else:
            max_radius = 0.0

        # base_limit = 2 * max_radius * self.rng.ase_cutoff_mult
        base_limit = 2 * max_radius * self.rng.ase_cutoff_mult

        # max_custom = max(custom_cutoffs.values()) if exists, else 0
        max_custom = 0.0
        if custom_cutoffs:
            max_custom = max(custom_cutoffs.values())

        # Final Cutoff: search_cutoff = max(base_limit, max_custom) + 0.45
        # Aligned with OpenBabel's convention.
        search_cutoff = max(base_limit, max_custom) + 0.45

        # Step 3: Neighbor Search
        i_list, j_list, d_list = neighbor_list("ijd", step_atoms, cutoff=search_cutoff)

        atomnumber = len(step_atoms)
        bond = [[] for _ in range(atomnumber)]
        bondlevel = [[] for _ in range(atomnumber)]

        # Step 4: Filtering & Deduplication
        # Iterate through results. Skip if i >= j
        for i, j, d in zip(i_list, j_list, d_list):
            if i >= j:
                continue

            # Threshold:
            # If pair (sym_i, sym_j) in custom_cutoffs: threshold = custom_value
            elem_i = symbols[i]
            elem_j = symbols[j]
            pair_key = frozenset({elem_i, elem_j})

            if pair_key in custom_cutoffs:
                threshold = custom_cutoffs[pair_key]
            else:
                # Else: threshold = (radii[i] + radii[j]) * self.rng.ase_cutoff_mult
                threshold = (radii[i] + radii[j]) * self.rng.ase_cutoff_mult

            # Bond Decision:
            # If d < threshold:
            if d < threshold:
                # Fix Duplicate Bonds: Check if j not in bond[i] before appending
                # This prevents RDKit "bond already exists" errors caused by PBC images
                if j not in bond[i]:
                    bond[i].append(j)
                    bond[j].append(i)
                    # Set bondlevel to 1
                    bondlevel[i].append(1)
                    bondlevel[j].append(1)

        return bond, bondlevel

    def _getbondfromcrd(
        self, step_atoms: Atoms, cell: np.ndarray
    ) -> Tuple[List[List[int]], List[List[int]]]:
        """Perceive bonds from atomic coordinates.

        Parameters
        ----------
        step_atoms : ase.Atoms
            Atoms in a step.
        cell : np.ndarray
            Cell in the shape (3, 3).

        Returns
        -------
        list[list[int]]
            Connected atoms for each atom.
        list[list[int]]
            Bond orders for each atom. 12 is an aromatic bond.
        """
        # Check if ASE mode is enabled
        if getattr(self.rng, "use_ase", False):
            return self._getbondfromase(step_atoms, cell)

        # Otherwise, use the original OpenBabel logic
        atomnumber = len(step_atoms)
        # Use openbabel to connect atoms
        mol = openbabel.OBMol()
        mol.BeginModify()
        for idx, (num, position) in enumerate(
            zip(step_atoms.get_atomic_numbers(), step_atoms.positions)
        ):
            a = mol.NewAtom(idx)
            a.SetAtomicNum(int(num))
            a.SetVector(*position)
        # Apply period boundary conditions
        # openbabel#1853, supported in v3.1.0
        if self.pbc:
            uc = openbabel.OBUnitCell()
            uc.SetData(
                openbabel.vector3(cell[0][0], cell[0][1], cell[0][2]),
                openbabel.vector3(cell[1][0], cell[1][1], cell[1][2]),
                openbabel.vector3(cell[2][0], cell[2][1], cell[2][2]),
            )
            mol.CloneData(uc)
            mol.SetPeriodicMol()
        mol.ConnectTheDots()
        mol.PerceiveBondOrders()
        mol.EndModify()
        bond = [[] for i in range(atomnumber)]
        bondlevel = [[] for i in range(atomnumber)]
        for b in openbabel.OBMolBondIter(mol):
            s1 = b.GetBeginAtom().GetId()
            s2 = b.GetEndAtom().GetId()
            level = b.GetBondOrder()
            if level == 5:
                # aromatic, 5 in openbabel but 12 in rdkit
                level = 12
            bond[s1].append(s2)
            bond[s2].append(s1)
            bondlevel[s1].append(level)
            bondlevel[s2].append(level)
        return bond, bondlevel


@_Detect.register_subclass("dump")
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
        logger.info("Reading from LAMMPS dump file")
        iscompleted = False
        N = None
        linecontent = None
        atomtype = None
        index = -1
        stepaindex = None
        for index, line in enumerate(f):
            if line.startswith("ITEM:"):
                linecontent = self.LineType.linecontent(line)
                if linecontent == self.LineType.ATOMS:
                    keys = line.split()
                    self.id_idx = keys.index("id") - 2
                    self.tidx = keys.index("type") - 2
                    self.xidx = keys.index("x") - 2
                    self.yidx = keys.index("y") - 2
                    self.zidx = keys.index("z") - 2
            else:
                assert linecontent is not None
                if linecontent == self.LineType.NUMBER:
                    if iscompleted:
                        assert stepaindex is not None
                        steplinenum = index - stepaindex
                        break
                    else:
                        iscompleted = True
                        stepaindex = index
                    N = int(line.split()[0])
                    atomtype = np.zeros(N, dtype=int)
                    logger.info(f"Found N = {N}, initializing atomtype array")
                elif linecontent == self.LineType.ATOMS:
                    s = line.split()
                    assert atomtype is not None
                    atomtype[int(s[self.id_idx]) - 1] = int(s[self.tidx]) - 1
        else:
            steplinenum = index + 1
        logger.info(f"Finished reading N, total lines processed: {index + 1}")
        assert N is not None and atomtype is not None
        self.N = N
        self.atomtype = atomtype
        return steplinenum

    def _readstepfunc(self, item) -> Tuple[List[bytes], Tuple[int, int]]:
        step, lines = item
        step_atoms = []
        ss = []
        linecontent = None
        timestep = None
        for line in lines:
            if line:
                if line.startswith("ITEM:"):
                    linecontent = self.LineType.linecontent(line)
                else:
                    if linecontent is None:
                        raise ValueError("LAMMPS dump file format error")
                    elif linecontent == self.LineType.ATOMS:
                        s = line.split()
                        step_atoms.append(
                            (
                                int(s[self.id_idx]),
                                Atom(
                                    self.atomname[int(s[self.tidx]) - 1],
                                    (
                                        float(s[self.xidx]),
                                        float(s[self.yidx]),
                                        float(s[self.zidx]),
                                    ),
                                ),
                            )
                        )
                    elif linecontent == self.LineType.TIMESTEP:
                        timestep = step, int(line.split()[0])
                    elif linecontent == self.LineType.BOX:
                        s = line.split()
                        ss.append(list(map(float, s)))
        assert timestep is not None
        ss = np.array(ss)
        if ss.shape[1] > 2:
            xy = ss[0][2]
            xz = ss[1][2]
            yz = ss[2][2]
        else:
            xy, xz, yz = 0.0, 0.0, 0.0
        xlo = ss[0][0] - min(0.0, xy, xz, xy + xz)
        xhi = ss[0][1] - max(0.0, xy, xz, xy + xz)
        ylo = ss[1][0] - min(0.0, yz)
        yhi = ss[1][1] - max(0.0, yz)
        zlo = ss[2][0]
        zhi = ss[2][1]
        boxsize = np.array(
            [
                [xhi - xlo, 0.0, 0.0],
                [xy, yhi - ylo, 0.0],
                [xz, yz, zhi - zlo],
            ]
        )
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
        N = None
        atomtype = None
        for index, line in enumerate(f):
            s = line.split()
            if index == 0:
                N = int(line.strip())
                atomtype = np.zeros(N, dtype=int)
            elif N is None:
                assert N is not None  # pragma: no cover
            elif index > N + 1:
                break
            elif index > 1:
                assert atomtype is not None
                atomtype[index - 2] = atomname_dict[s[0]]
        assert N is not None and atomtype is not None
        self.N = N
        self.atomtype = atomtype
        steplinenum = N + 2
        return steplinenum

    def _readstepfunc(self, item) -> Tuple[List[bytes], Tuple[int, int]]:
        step, lines = item
        timestep = step, step
        step_atoms = []
        boxsize = self.cell
        if self.pbc and boxsize is None:
            raise RuntimeError("No cell information is given.")
        for index, line in enumerate(lines):
            s = line.split()
            if index > 1:
                step_atoms.append(
                    (index - 1, Atom(s[0], tuple(float(x) for x in s[1:4])))
                )
        _, step_atoms = zip(*sorted(step_atoms, key=operator.itemgetter(0)))
        step_atoms = Atoms(step_atoms)
        bond, level = self._getbondfromcrd(step_atoms, boxsize)
        molecules = self._connectmolecule(bond, level)
        return molecules, timestep


@_Detect.register_subclass("extxyz")
class _Detectextxyz(_DetectCrd):
    """extxyz file. xyz with extended metadata support like cell, force, etc."""

    def _readNfunc(self, f):
        atomname_dict = dict(zip(self.atomname.tolist(), range(self.atomname.size)))
        N = None
        atomtype = None
        for index, line in enumerate(f):
            if index == 0:
                N = int(line.strip())
                atomtype = np.zeros(N, dtype=int)
            elif index == 1:
                # Extract lattice (cell) from comment line
                if "Lattice=" in line:
                    lattice_str = line.split("Lattice=")[1].split('"')[1]
                    lattice_floats = list(map(float, lattice_str.split()))
                    self.cell = np.array(lattice_floats).reshape((3, 3))
                else:
                    raise RuntimeError("Missing Lattice= information in extxyz.")
            elif (N is not None) and (index > N + 1):
                break
            elif index > 1:
                s = line.split()
                assert atomtype is not None
                atomtype[index - 2] = atomname_dict[s[0]]
        assert N is not None and atomtype is not None
        self.N = N
        self.atomtype = atomtype
        steplinenum = N + 2
        return steplinenum

    def _readstepfunc(self, item) -> Tuple[List[bytes], Tuple[int, int]]:
        step, lines = item
        step_atoms = []
        timestep = step, step  # Use step as timestep fallback
        boxsize = self.cell
        if self.pbc and boxsize is None:
            raise RuntimeError("No cell information is given in extxyz.")
        for index, line in enumerate(lines):
            if index > 1:
                s = line.split()
                step_atoms.append(
                    (index - 1, Atom(s[0], tuple(float(x) for x in s[1:4])))
                )
        _, step_atoms = zip(*sorted(step_atoms, key=operator.itemgetter(0)))
        step_atoms = Atoms(step_atoms)
        bond, level = self._getbondfromcrd(step_atoms, boxsize)
        molecules = self._connectmolecule(bond, level)
        return molecules, timestep
