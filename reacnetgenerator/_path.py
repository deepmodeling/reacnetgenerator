# SPDX-License-Identifier: LGPL-3.0-or-later
# cython: language_level=3
# cython: linetrace=True
"""Collect paths.

To produce a reaction network, every molecule (species) should be treated as a
node in the network. Therefore, all detected species are indexed by canonical
SMILES to guarantee its uniqueness. Isomers are also identified according to
SMILES codes._[1] The VF2 algorithm can be also used to identify isomers, which is
an option in ReacNetGenerator._[2] After filtering out noise, the reaction path of atoms
and the number of intermolecular reactions can be calculated.

References
----------
.. [1] Landrum, G. RDKit: Open-Source Cheminformatics Software 2016.
.. [2] Cordella, L. P.; Foggia, P.; Sansone, C.; Vento, M. A (Sub)Graph
   Isomorphism Algorith for Matching Large Graphs. IEEE Trans. Pattern Analysis
   and Machine Intelligence 2004, 26, 1367-1372.
"""

import itertools
import re
from abc import ABCMeta, abstractmethod
from collections import Counter, defaultdict

import networkx as nx
import networkx.algorithms.isomorphism as iso
import numpy as np
from rdkit import Chem
from tqdm.auto import tqdm

from ._reaction import ReactionsFinder
from .utils import (
    SharedRNGData,
    WriteBuffer,
    bytestolist,
    get_timestep_value,
    listtostirng,
    read_compressed_block,
    run_mp,
)


class _CollectPaths(SharedRNGData, metaclass=ABCMeta):
    runHMM: bool
    N: int
    step: int
    atomname: np.ndarray
    originfilename: str
    hmmfilename: str
    moleculefilename: str
    moleculetemp2filename: str
    atomroutefilename: str
    nproc: int
    hmmit: int
    atomtype: np.ndarray
    selectatoms: list
    split: int
    miso: int
    timestep: dict
    printmoleculetime: bool
    moleculeframes: list
    moleculetimesteps: list
    mname: np.ndarray

    def __init__(self, rng):
        SharedRNGData.__init__(
            self,
            rng,
            [
                "runHMM",
                "N",
                "step",
                "atomname",
                "originfilename",
                "hmmfilename",
                "moleculefilename",
                "moleculetemp2filename",
                "atomroutefilename",
                "nproc",
                "hmmit",
                "atomtype",
                "selectatoms",
                "split",
                "miso",
                "timestep",
                "printmoleculetime",
                "moleculeframes",
                "moleculetimesteps",
            ],
            ["mname", "atomnames", "allmoleculeroute", "splitmoleculeroute"],
        )

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
            self.splitmoleculeroute = [
                self._printatomroute(atomeach[:, st], timeaxis=i)
                for i, st in enumerate(splittime)
            ]
        self.returnkeys()
        ReactionsFinder(self.rng).findreactions(atomeach.T, conflict.T)

    @abstractmethod
    def _printmoleculename(self):
        pass

    def _getatomeach(self):
        """Values in atomeach starts from 1."""
        atomeach = np.zeros((self.N, self.step), dtype=int)
        conflict = np.zeros((self.N, self.step), dtype=int)
        with open(
            self.hmmfilename if self.runHMM else self.originfilename, "rb"
        ) as fh, open(self.moleculetemp2filename, "rb") as ft:
            for i, (linehz, linetz) in enumerate(
                tqdm(
                    zip(
                        read_compressed_block(fh),
                        itertools.zip_longest(*[read_compressed_block(ft)] * 4),
                    ),
                    total=self.hmmit,
                    desc="Analyze atoms",
                    unit="molecule",
                    disable=None,
                ),
                start=1,
            ):
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
        if atomeachi.size:
            time = np.concatenate(
                [
                    np.zeros((1,), dtype=int),
                    np.nonzero(np.diff(atomeachi))[0] + 1,
                ]
            )
            route = atomeachi[time]
        else:
            time = np.zeros(0, dtype=int)
            route = np.zeros(0, dtype=int)
        moleculeroute = (
            np.dstack((route[:-1], route[1:]))[0]
            if self.atomname[atomtypei] in self.selectatoms
            else np.zeros((0, 2), dtype=int)
        )
        names = self.mname[route - 1]
        # Atom {idx}: {time} {SMILES} -> {time} {SMILES} -> ...
        routestr = f"Atom {i} {self.atomname[atomtypei]}: " + " -> ".join(
            [f"{tt} {name}" for tt, name in zip(time, names)]
        )
        return moleculeroute, routestr

    def _printatomroute(self, atomeach, timeaxis=None):
        """For analysis without HMM, we may not need to use np.unique."""
        with WriteBuffer(
            open(
                (
                    self.atomroutefilename
                    if timeaxis is None
                    else f"{self.atomroutefilename}.{timeaxis}"
                ),
                "w",
            ),
            sep="\n",
        ) as f:
            allmoleculeroute = []
            if not self.runHMM:
                have_added = {}
            else:
                have_added = None
            results = run_mp(
                self.nproc,
                func=self._getatomroute,
                l=zip(atomeach, self.atomtype),
                return_num=True,
                start=1,
                unordered=False,
                total=self.N,
                desc=(
                    "Collect reaction paths"
                    if timeaxis is None
                    else f"Collect reaction paths {timeaxis}"
                ),
                unit="atom",
            )
            for ii, (moleculeroute, routestr) in enumerate(results):
                f.append(routestr)
                if moleculeroute.size > 0:
                    if not self.runHMM:
                        # check whether repeated or not if analyzing without HMM
                        for rr in moleculeroute:
                            tpr = tuple(rr)
                            assert have_added is not None
                            if have_added.get(tpr, atomeach.shape[0]) >= ii:
                                have_added[tpr] = ii
                                allmoleculeroute.append(rr.reshape(1, 2))
                    else:
                        allmoleculeroute.append(moleculeroute)
        allmoleculeroute = (
            np.concatenate(allmoleculeroute)
            if allmoleculeroute
            else np.zeros((0, 2), dtype=int)
        )
        if self.runHMM and allmoleculeroute.size:
            allmoleculeroute = np.unique(allmoleculeroute, axis=0)
        return allmoleculeroute

    def _re(self, smi):
        """If you use RDkit to convert a methyl radical to SMILES, you will get something
        like [H]C([H])[H]. However, OpenBabel will consider it as a methane molecule. So,
        you have to use [H][C]([H])[H], if you need to process some radicals.

        Examples
        --------
        >>> self._re('C')
        [C]
        >>> self._re('[C]')
        [C]
        >>> self._re('[CH]')
        [CH]
        >>> self._re('Na')
        [Na]
        >>> self._re('[H]c(Cl)C([H])Cl')
        [H][c]([Cl])[C]([H])[Cl]
        """
        if "_unknownSMILES" in smi:
            # not SMILES
            return smi
        Satom = sorted(self.atomname, key=lambda x: len(x), reverse=True)
        elements = "|".join(
            [
                ((an.upper() + "|" + an.lower()) if len(an) == 1 else an)
                for an in Satom
                if an != "H"
            ]
        )
        smi = re.sub(r"(?<!\[)(" + elements + r")(?!H|\])", r"[\1]", smi)
        return smi.replace("[HH]", "[H]")

    def convertSMILES(self, atoms, bonds):
        """Convert atoms and bonds information to SMILES.

        Raises
        ------
        ValueError
            (RDKit error) Maximum BFS search size exceeded.
        """
        m = Chem.RWMol(Chem.MolFromSmiles(""))
        d = {}
        for name, number in zip(self.atomnames[atoms], atoms):
            d[number] = m.AddAtom(Chem.Atom(name))
        for atom1, atom2, level in bonds:
            m.AddBond(d[atom1], d[atom2], Chem.BondType(level))
        # https://github.com/rdkit/rdkit/discussions/6613#discussioncomment-6688021
        for a in m.GetAtoms():
            a.SetNoImplicit(True)
        name = Chem.MolToSmiles(m)
        return self._re(name)

    def _getatomsandbonds(self, line):
        atoms = np.array(bytestolist(line[0]), dtype=int)
        pairs = bytestolist(line[1])
        levels = bytestolist(line[2])
        bonds = [[*pair, level] for pair, level in zip(pairs, levels)]
        return atoms, bonds

    def _getmoleculeframes(self, line):
        return np.array(bytestolist(line[-1]), dtype=int)

    def _needmoleculeframes(self):
        return (
            self.printmoleculetime
            or self._hasmoleculefilter(self.moleculeframes)
            or self._hasmoleculefilter(self.moleculetimesteps)
        )

    def _getmoleculeframesandtimesteps(self, line, need_timesteps=True):
        if not self._needmoleculeframes():
            return None, None
        frames = self._getmoleculeframes(line)
        timesteps = (
            self._getmoleculetimesteps(frames)
            if need_timesteps
            and (
                self.printmoleculetime
                or self._hasmoleculefilter(self.moleculetimesteps)
            )
            else None
        )
        return frames, timesteps

    def _getmoleculetimesteps(self, frames):
        return [get_timestep_value(self.timestep[int(frame)]) for frame in frames]

    @staticmethod
    def _hasmoleculefilter(values):
        return values is not None and len(values) > 0

    def _shouldprintmolecule(self, frames, timesteps=None):
        has_frame_filter = self._hasmoleculefilter(self.moleculeframes)
        has_timestep_filter = self._hasmoleculefilter(self.moleculetimesteps)
        if not has_frame_filter and not has_timestep_filter:
            return True
        assert frames is not None
        if has_timestep_filter:
            if timesteps is None:
                timesteps = self._getmoleculetimesteps(frames)
            timestep_filter = set(self.moleculetimesteps)
        if has_frame_filter:
            frame_filter = set(self.moleculeframes)
        if has_frame_filter and has_timestep_filter:
            return any(
                int(frame) in frame_filter and timestep in timestep_filter
                for frame, timestep in zip(frames, timesteps)
            )
        if has_frame_filter:
            return bool(set(map(int, frames)).intersection(frame_filter))
        return bool(set(timesteps).intersection(timestep_filter))

    def _formatmoleculename(self, name, atoms, bonds, frames=None, timesteps=None):
        if self.printmoleculetime:
            assert frames is not None
            if timesteps is None:
                timesteps = self._getmoleculetimesteps(frames)
            return listtostirng(
                (name, atoms, bonds, frames, timesteps), sep=(" ", ";", ",")
            )
        return listtostirng((name, atoms, bonds), sep=(" ", ";", ","))


class _CollectMolPaths(_CollectPaths):
    """VF2 is used to identify isomers.

    If SMILES is failed to generate, fallback to the name like CxHyOz.
    """

    def _printmoleculename(self):
        mname = []
        d = defaultdict(list)
        em = iso.numerical_edge_match(["atom", "level"], ["None", 1])
        # idx for unknown SMILES
        self.n_unknown = 0
        with WriteBuffer(open(self.moleculefilename, "w"), sep="\n") as fm, open(
            self.moleculetemp2filename, "rb"
        ) as ft:
            for line in tqdm(
                itertools.zip_longest(*[read_compressed_block(ft)] * 4),
                total=self.hmmit,
                desc="Indentify isomers",
                unit="molecule",
                disable=None,
            ):
                atoms, bonds = self._getatomsandbonds(line)
                frames, timesteps = self._getmoleculeframesandtimesteps(line)
                molecule = Molecule(self, atoms, bonds)
                for isomer in d[str(molecule)]:
                    if isomer.isomorphic(molecule, em):
                        molecule.smiles = isomer.smiles
                        break
                else:
                    d[str(molecule)].append(molecule)
                mname.append(molecule.smiles)
                if self._shouldprintmolecule(frames, timesteps):
                    fm.append(
                        self._formatmoleculename(
                            molecule.smiles, atoms, bonds, frames, timesteps
                        )
                    )
        self.mname = np.array(mname)


class _CollectSMILESPaths(_CollectPaths):
    def _printmoleculename(self):
        mname = []
        d = defaultdict(list)
        name_mapping = {}
        name_mapping_graph = defaultdict(dict)
        em = iso.numerical_edge_match(["atom", "level"], ["None", 1])
        self.n_unknown = 0
        with WriteBuffer(open(self.moleculefilename, "w"), sep="\n") as fm, open(
            self.moleculetemp2filename, "rb"
        ) as ft:
            results = run_mp(
                self.nproc,
                func=self._calmoleculeSMILESname,
                l=read_compressed_block(ft),
                unordered=False,
                nlines=4,
                total=self.hmmit,
                desc="Indentify isomers",
                unit="molecule",
            )
            for name, atoms, bonds, frames in results:
                timesteps = (
                    self._getmoleculetimesteps(frames)
                    if frames is not None
                    and (self.printmoleculetime or self.moleculetimesteps is not None)
                    else None
                )
                if name is None:
                    # SMILES failed, fallback to VF2 identify isomers
                    molecule = Molecule(self, atoms, bonds)

                    # directly raise ValueError to save time
                    def _raise_anyway(*args, **kwargs):
                        raise ValueError("Maximum BFS search size exceeded.")

                    molecule._convertSMILES = _raise_anyway
                    for isomer in d[str(molecule)]:
                        if isomer.isomorphic(molecule, em):
                            molecule.smiles = isomer.smiles
                            break
                    else:
                        d[str(molecule)].append(molecule)
                    name = molecule.smiles
                if self.miso > 0:
                    if name in name_mapping:
                        name = name_mapping[name]
                    else:
                        # check if the name is isomorphic to the previous molecules
                        molecule = Molecule(self, atoms, bonds)
                        # the formula should be the same
                        mng = name_mapping_graph[molecule.name]
                        for isomer, mol in mng.items():
                            if mol.isomorphic(molecule, em):
                                # use the previous SMILES
                                name_mapping[name] = isomer
                                name = isomer
                                break
                        else:
                            mng[name] = molecule
                            name_mapping[name] = name
                mname.append(name)
                if self._shouldprintmolecule(frames, timesteps):
                    fm.append(
                        self._formatmoleculename(name, atoms, bonds, frames, timesteps)
                    )
        self.mname = np.array(mname)

    def _calmoleculeSMILESname(self, item):
        line = item
        atoms, bonds = self._getatomsandbonds(line)
        frames, _ = self._getmoleculeframesandtimesteps(line, need_timesteps=False)
        try:
            name = self.convertSMILES(atoms, bonds)
        except ValueError:
            # fallback to VF2
            name = None
        return name, atoms, bonds, frames


class Molecule:
    """A molecule class for isomer identification."""

    def __init__(self, cmp, atoms, bonds):
        self.cmp = cmp
        self.atoms = atoms
        self.bonds = bonds
        self._atomtypes = cmp.atomtype[atoms]
        self._atomnames = cmp.atomnames[atoms]
        self._miso = cmp.miso
        self.graph = self._makemoleculegraph()
        counter = Counter(self._atomnames)
        self.name = "".join(
            f"{atomname}{counter[atomname]}" for atomname in cmp.atomname
        )
        self._smiles = None
        self._convertSMILES = cmp.convertSMILES

    def __str__(self):
        return self.name

    @property
    def smiles(self):
        """Return SMILES of a molecule."""
        if self._smiles is None:
            try:
                self._smiles = self._convertSMILES(self.atoms, self.bonds)
            except ValueError:
                # when RDKit error: Maximum BFS search size exceeded
                # fallback to the name of the molecule
                # blank should be avoided
                self._smiles = self.name + f"_unknownSMILES_{self.cmp.n_unknown}"
                self.cmp.n_unknown += 1
        return self._smiles

    @smiles.setter
    def smiles(self, value):
        self._smiles = value

    def _makemoleculegraph(self):
        graph = nx.Graph()
        for line in self.bonds:
            if self._miso == 0:
                # normal mode
                graph.add_edge(line[0], line[1], level=line[2])
            elif self._miso == 1:
                # merge the isomers with same atoms and same bond-network but different bond orders
                graph.add_edge(line[0], line[1], level=1)
            elif self._miso == 2:
                # merge the isomers with same atoms with different bond-network
                pass
            else:
                raise ValueError(f"Unknown isomer identification method: {self._miso}.")
        for atomnumber, atomtype in zip(self.atoms, self._atomtypes):
            graph.add_node(atomnumber, atom=atomtype)
        return graph

    def isomorphic(self, mol, em):
        """Return whether two molecules are isomorphic."""
        return nx.is_isomorphic(self.graph, mol.graph, em)
