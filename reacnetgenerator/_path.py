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

import csv
import heapq
import itertools
import os
import re
import shutil
import tempfile
from abc import ABCMeta, abstractmethod
from collections import Counter, defaultdict
from contextlib import ExitStack, closing
from multiprocessing.util import Finalize

import networkx as nx
import networkx.algorithms.isomorphism as iso
import numpy as np
from rdkit import Chem
from tqdm.auto import tqdm

from ._reaction import ReactionsFinder
from ._step3state import (
    _STEP3_SCAN_ROWS,
    _AtomFrameReader,
    _AtomFrameStore,
    _close_mapping,
    _MoleculeNameBuilder,
    _MoleculeNameTable,
)
from .utils import (
    SharedRNGData,
    WriteBuffer,
    bytestolist,
    get_timestep_value,
    listtostirng,
    read_compressed_block,
    run_mp,
)

_ROUTE_WORKER_STATE = None


def _route_changes(timeline, active_transitions=None, block_rows=_STEP3_SCAN_ROWS):
    """Compress a route in bounded blocks while marking raw-frame changes.

    Legacy route positions count only nonzero observations. Reaction indices
    instead refer to adjacent original frames, including transitions to/from 0.
    Keep these two coordinate systems separate across block boundaries.
    """
    if block_rows <= 0:
        raise ValueError("Scan block size must be positive")
    times, routes = [], []
    nonzero_count = 0
    previous = 0
    for start in range(0, len(timeline), block_rows):
        stop = min(start + block_rows, len(timeline))
        if active_transitions is not None:
            first = max(1, start)
            changed = timeline[first:stop] != timeline[first - 1 : stop - 1]
            active_transitions[first - 1 : stop - 1][changed] = 1
        block = timeline[start:stop]
        values = block[block != 0]
        if not len(values):
            continue
        changed = np.empty(len(values), dtype=np.bool_)
        changed[0] = values[0] != previous
        changed[1:] = values[1:] != values[:-1]
        indices = np.flatnonzero(changed)
        if len(indices):
            times.append(nonzero_count + indices)
            routes.append(values[indices])
        nonzero_count += len(values)
        previous = values[-1]
    if not routes:
        return np.zeros(0, dtype=int), np.zeros(0, dtype=int)
    return np.concatenate(times), np.concatenate(routes)


def _atom_route_result(atom_index, timeline, atom_name, selected, names, active=None):
    """Format one atom's route with the historical nonzero-position semantics."""
    time, route = _route_changes(timeline, active)
    molecule_route = (
        np.column_stack((route[:-1], route[1:]))
        if selected
        else np.zeros((0, 2), dtype=int)
    )
    route_string = f"Atom {atom_index + 1} {atom_name}: " + " -> ".join(
        f"{tt} {name}" for tt, name in zip(time, names[route - 1])
    )
    return molecule_route, route_string


def _initialize_route_worker(
    reader_args, atomtype, atomname, selectatoms, names, frame_range, active_path
):
    """Attach mappings once; no trajectory arrays travel with atom tasks."""
    global _ROUTE_WORKER_STATE
    reader = _AtomFrameReader(*reader_args)
    Finalize(None, reader.close, exitpriority=10)
    active = None
    if active_path is not None:
        active = np.memmap(
            active_path,
            mode="r+",
            dtype=np.uint8,
            shape=(reader.atomeach.shape[1] - 1,),
        )
        Finalize(None, _close_mapping, args=(active,), exitpriority=10)
    _ROUTE_WORKER_STATE = (
        reader,
        atomtype,
        atomname,
        selectatoms,
        names,
        frame_range,
        active,
    )


def _get_atom_route_by_index(atom_index):
    """Process a zero-based atom index using this worker's attached matrices."""
    assert _ROUTE_WORKER_STATE is not None
    reader, atomtype, atomname, selectatoms, names, frame_range, active = (
        _ROUTE_WORKER_STATE
    )
    start, stop = frame_range
    name = atomname[atomtype[atom_index]]
    return _atom_route_result(
        atom_index,
        reader.atomeach[atom_index, start:stop],
        name,
        name in selectatoms,
        names,
        active,
    )


class _MoleculeTimelineSpool:
    def __init__(self, filename, buffer_rows=10000, max_open_chunks=64):
        self.filename = filename
        self.buffer_rows = max(1, int(buffer_rows))
        self.max_open_chunks = max(2, int(max_open_chunks))
        output_dir = os.path.dirname(os.path.abspath(filename))
        self.tempdir = tempfile.mkdtemp(prefix=".molecule-timeline-", dir=output_dir)
        self.paths = []
        self.buffer = []
        self.n_chunks = 0

    def append(self, row):
        self.buffer.append(row)
        if len(self.buffer) >= self.buffer_rows:
            self.flush()

    def extend(self, rows):
        for row in rows:
            self.append(row)

    def flush(self):
        if not self.buffer:
            return
        self.buffer.sort(key=self._sortkey)
        path = self._newpath()
        with open(path, "w", newline="") as f:
            csv.writer(f).writerows(self.buffer)
        self.paths.append(path)
        self.buffer = []

    def write(self):
        self.flush()
        self.paths = self._mergechunks(self.paths)
        with open(self.filename, "w", newline="") as timeline_file:
            timeline_writer = csv.writer(timeline_file)
            timeline_writer.writerow(["Timestep", "Species", "AtomIDs", "BondIDs"])
            timeline_writer.writerows(self._itermergedrows(self.paths))

    def close(self):
        shutil.rmtree(self.tempdir, ignore_errors=True)

    @staticmethod
    def _sortkey(row):
        return int(row[0])

    def _newpath(self):
        path = os.path.join(self.tempdir, f"{self.n_chunks}.csv")
        self.n_chunks += 1
        return path

    def _mergechunks(self, paths):
        while len(paths) > self.max_open_chunks:
            merged_paths = []
            for i in range(0, len(paths), self.max_open_chunks):
                group = paths[i : i + self.max_open_chunks]
                if len(group) == 1:
                    merged_paths.append(group[0])
                    continue
                path = self._newpath()
                with open(path, "w", newline="") as f:
                    csv.writer(f).writerows(self._itermergedrows(group))
                for old_path in group:
                    os.remove(old_path)
                merged_paths.append(path)
            paths = merged_paths
        return paths

    def _itermergedrows(self, paths):
        if not paths:
            return
        with ExitStack() as stack:
            readers = []
            for path in paths:
                readers.append(csv.reader(stack.enter_context(open(path, newline=""))))
            yield from heapq.merge(*readers, key=self._sortkey)


class _CollectPaths(SharedRNGData, metaclass=ABCMeta):
    runHMM: bool
    N: int
    step: int
    atomname: np.ndarray
    originfilename: str
    hmmfilename: str
    moleculefilename: str
    moleculetimelinefilename: str
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
    mname: _MoleculeNameTable
    _moleculetimelinebufferrows: int = 10000

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
                "moleculetimelinefilename",
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
        self._moleculeframefilter = self._getmoleculefilterset(self.moleculeframes)
        self._moleculetimestepfilter = self._getmoleculefilterset(
            self.moleculetimesteps
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
        with self._getatomeach() as matrix_store:
            atomeach = matrix_store.atomeach
            matrix_store.prepare_active_transitions()
            self.allmoleculeroute = self._printatomroute(matrix_store)
            if self.split > 1:
                # Integer ranges match array_split, including empty trailing
                # splits, without allocating frame indices or copying matrices.
                size, remainder = divmod(self.step, self.split)
                start = 0
                self.splitmoleculeroute = []
                for i in range(self.split):
                    stop = start + size + (i < remainder)
                    self.splitmoleculeroute.append(
                        self._printatomroute(
                            matrix_store, timeaxis=i, frame_range=(start, stop)
                        )
                    )
                    start = stop
            self.returnkeys()
            matrix_store.flush()
            ReactionsFinder(self.rng).findreactions(
                atomeach.T,
                matrix_store.conflict.T,
                matrix_store=matrix_store,
            )

    @abstractmethod
    def _printmoleculename(self):
        pass

    def _getatomeach(self):
        """Build compact atom-frame matrices; molecule IDs start from 1."""
        try:
            store = _AtomFrameStore(
                (self.N, self.step),
                self.hmmit,
                directory=os.path.dirname(os.path.abspath(self.atomroutefilename)),
            )
        except PermissionError:
            # Writable outputs such as /dev/null can have unwritable parents.
            store = _AtomFrameStore((self.N, self.step), self.hmmit)
        try:
            with (
                open(
                    self.hmmfilename if self.runHMM else self.originfilename,
                    "rb",
                ) as fh,
                open(self.moleculetemp2filename, "rb") as ft,
            ):
                for molecule_id, (linehz, linetz) in enumerate(
                    tqdm(
                        zip(
                            read_compressed_block(fh),
                            itertools.zip_longest(*[read_compressed_block(ft)] * 4),
                            strict=True,
                        ),
                        total=self.hmmit,
                        desc="Analyze atoms",
                        unit="molecule",
                        disable=None,
                    ),
                    start=1,
                ):
                    if molecule_id > self.hmmit:
                        raise RuntimeError(
                            "More molecule records than the declared count"
                        )
                    lineh = bytestolist(linehz)
                    atoms = np.asarray(bytestolist(linetz[0]), dtype=np.int64)
                    frames = np.flatnonzero(lineh)
                    if frames.size:
                        overlap = np.not_equal(
                            store.atomeach[atoms[:, None], frames],
                            0,
                        )
                        store.conflict.mark(atoms, frames, overlap)
                        store.atomeach[atoms[:, None], frames] = molecule_id
            store.flush()
            return store
        except BaseException:
            store.close()
            raise

    def _getatomroute(self, item):
        """Adapt a legacy one-based atom task to the shared route formatter."""
        i, (atomeachi, atomtypei) = item
        name = self.atomname[atomtypei]
        return _atom_route_result(
            i - 1, atomeachi, name, name in self.selectatoms, self.mname
        )

    def _printatomroute(self, matrix_store, timeaxis=None, frame_range=None):
        """For analysis without HMM, we may not need to use np.unique."""
        if frame_range is None:
            frame_range = (0, self.step)
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
                func=_get_atom_route_by_index,
                l=range(self.N),
                initializer=_initialize_route_worker,
                initargs=(
                    matrix_store.reader_args,
                    self.atomtype,
                    self.atomname,
                    self.selectatoms,
                    self.mname,
                    frame_range,
                    matrix_store.active_transition_path if timeaxis is None else None,
                ),
                unordered=False,
                chunksize=1,
                max_inflight=max(2, 2 * self.nproc),
                disk_ordered=True,
                total=self.N,
                desc=(
                    "Collect reaction paths"
                    if timeaxis is None
                    else f"Collect reaction paths {timeaxis}"
                ),
                unit="atom",
            )
            # Finish or terminate workers before the parent releases mappings,
            # including when writing a route or consuming a result fails.
            with closing(results):
                for ii, (moleculeroute, routestr) in enumerate(results):
                    f.append(routestr)
                    if moleculeroute.size > 0:
                        if not self.runHMM:
                            # Deduplicate atom routes exactly as in the legacy path.
                            for rr in moleculeroute:
                                tpr = tuple(rr)
                                assert have_added is not None
                                if have_added.get(tpr, self.N) >= ii:
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

    def _needmoleculetimeline(self):
        return (
            self.printmoleculetime
            or self._moleculeframefilter is not None
            or self._moleculetimestepfilter is not None
        )

    def _getmoleculeframesandtimesteps(self, line, need_timesteps=True):
        if not self._needmoleculetimeline():
            return None, None
        frames = self._getmoleculeframes(line)
        timesteps = (
            self._getmoleculetimesteps(frames)
            if need_timesteps
            and (self.printmoleculetime or self._moleculetimestepfilter is not None)
            else None
        )
        return frames, timesteps

    def _getmoleculetimesteps(self, frames):
        return [get_timestep_value(self.timestep[int(frame)]) for frame in frames]

    @staticmethod
    def _hasmoleculefilter(values):
        return values is not None and len(values) > 0

    @classmethod
    def _getmoleculefilterset(cls, values):
        return set(values) if cls._hasmoleculefilter(values) else None

    def _shouldprintmoleculetimelinerow(self, frame, timestep):
        if (
            self._moleculeframefilter is not None
            and int(frame) not in self._moleculeframefilter
        ):
            return False
        return not (
            self._moleculetimestepfilter is not None
            and int(timestep) not in self._moleculetimestepfilter
        )

    def _formatmoleculename(self, name, atoms, bonds):
        return listtostirng((name, atoms, bonds), sep=(" ", ";", ","))

    @staticmethod
    def _formatmoleculeatomids(atoms):
        return ";".join(str(atom) for atom in atoms)

    @staticmethod
    def _formatmoleculebondids(bonds):
        return ";".join("-".join(str(item) for item in bond) for bond in bonds)

    def _getmoleculetimelinerows(self, name, atoms, bonds, frames, timesteps):
        assert frames is not None
        atom_ids = self._formatmoleculeatomids(atoms)
        bond_ids = self._formatmoleculebondids(bonds)
        if timesteps is None:
            frame_timestep_pairs = (
                (frame, get_timestep_value(self.timestep[int(frame)]))
                for frame in frames
            )
        else:
            frame_timestep_pairs = zip(frames, timesteps)
        for frame, timestep in frame_timestep_pairs:
            timestep = int(timestep)
            if self._shouldprintmoleculetimelinerow(frame, timestep):
                yield (timestep, name, atom_ids, bond_ids)

    def _writemoleculetimeline(self, rows):
        if not self._needmoleculetimeline():
            return
        timeline = self._openmoleculetimelinespool()
        try:
            timeline.extend(rows)
            timeline.write()
        finally:
            timeline.close()

    def _openmoleculetimelinespool(self):
        return _MoleculeTimelineSpool(
            self.moleculetimelinefilename,
            buffer_rows=self._moleculetimelinebufferrows,
        )


class _CollectMolPaths(_CollectPaths):
    """VF2 is used to identify isomers.

    If SMILES is failed to generate, fallback to the name like CxHyOz.
    """

    def _printmoleculename(self):
        mname = _MoleculeNameBuilder(self.hmmit)
        d = defaultdict(list)
        em = iso.numerical_edge_match(["atom", "level"], ["None", 1])
        # idx for unknown SMILES
        self.n_unknown = 0
        timeline = (
            self._openmoleculetimelinespool() if self._needmoleculetimeline() else None
        )
        try:
            with (
                WriteBuffer(open(self.moleculefilename, "w"), sep="\n") as fm,
                open(self.moleculetemp2filename, "rb") as ft,
            ):
                for line in tqdm(
                    itertools.zip_longest(*[read_compressed_block(ft)] * 4),
                    total=self.hmmit,
                    desc="Indentify isomers",
                    unit="molecule",
                    disable=None,
                ):
                    atoms, bonds = self._getatomsandbonds(line)
                    frames, timesteps = self._getmoleculeframesandtimesteps(
                        line, need_timesteps=False
                    )
                    molecule = Molecule(self, atoms, bonds)
                    for isomer in d[str(molecule)]:
                        if isomer.isomorphic(molecule, em):
                            molecule.smiles = isomer.smiles
                            break
                    else:
                        d[str(molecule)].append(molecule)
                    mname.append(molecule.smiles)
                    fm.append(self._formatmoleculename(molecule.smiles, atoms, bonds))
                    if timeline is not None:
                        timeline.extend(
                            self._getmoleculetimelinerows(
                                molecule.smiles,
                                atoms,
                                bonds,
                                frames,
                                timesteps,
                            )
                        )
            if timeline is not None:
                timeline.write()
        finally:
            if timeline is not None:
                timeline.close()
        self.mname = mname.finish()


class _CollectSMILESPaths(_CollectPaths):
    def _printmoleculename(self):
        mname = _MoleculeNameBuilder(self.hmmit)
        d = defaultdict(list)
        name_mapping = {}
        name_mapping_graph = defaultdict(dict)
        em = iso.numerical_edge_match(["atom", "level"], ["None", 1])
        self.n_unknown = 0
        timeline = (
            self._openmoleculetimelinespool() if self._needmoleculetimeline() else None
        )
        try:
            with (
                WriteBuffer(open(self.moleculefilename, "w"), sep="\n") as fm,
                open(self.moleculetemp2filename, "rb") as ft,
            ):
                results = run_mp(
                    self.nproc,
                    func=self._calmoleculeSMILESname,
                    l=read_compressed_block(ft),
                    unordered=False,
                    nlines=4,
                    chunksize=1,
                    max_inflight=max(2, 2 * self.nproc),
                    disk_ordered=True,
                    total=self.hmmit,
                    desc="Indentify isomers",
                    unit="molecule",
                )
                for name, atoms, bonds, frames in results:
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
                    fm.append(self._formatmoleculename(name, atoms, bonds))
                    if timeline is not None:
                        timeline.extend(
                            self._getmoleculetimelinerows(
                                name, atoms, bonds, frames, None
                            )
                        )
            if timeline is not None:
                timeline.write()
        finally:
            if timeline is not None:
                timeline.close()
        self.mname = mname.finish()

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
