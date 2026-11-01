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
from contextlib import ExitStack, closing, nullcontext
from multiprocessing.util import Finalize

import networkx as nx
import networkx.algorithms.isomorphism as iso
import numpy as np
from rdkit import Chem
from tqdm.auto import tqdm

from ._moleculepresence import _iter_molecule_presence
from ._reaction import ReactionsFinder
from ._smilesworker import _run_smiles_batches, _should_batch_smiles
from ._species import species_molecule
from ._step3state import (
    _STEP3_SCAN_ROWS,
    _AtomFrameReader,
    _AtomFrameStore,
    _close_mapping,
    _MoleculeNameBuilder,
    _MoleculeNameTable,
)
from ._timedoutput import _TimedOutputWriter
from .utils import (
    SharedRNGData,
    WriteBuffer,
    bytestolist,
    get_timestep_value,
    listtostirng,
    read_compressed_block,
    run_mp,
)

_ROUTE_COMPACTION_MIN_ROWS = 256
_ROUTE_COMPACTION_SAMPLE_ROWS = 1024
_ROUTE_COMPACTION_MAX_PAIRS = 65536
_ROUTE_WRITE_BUFFER_BYTES = 1024 * 1024

_ROUTE_WORKER_STATE = None
# Starting workers and attaching three mappings costs more than scanning a
# single block in the parent. Larger workloads retain the indexed worker path.
_STEP3_SERIAL_ATOM_FRAMES = _STEP3_SCAN_ROWS
_ATOM_MATCH = iso.categorical_node_match("atom", None)
_BOND_MATCH = iso.categorical_edge_match("level", None)
_MATRIX_WRITE_CELLS = 1024 * 1024


def _write_molecule_presence(store, atoms, presence, molecule_id):
    """Bound matrix gathers in both dimensions, preserving overlap conflicts."""
    atoms = np.asarray(atoms, dtype=np.int64)
    # Duplicate atom IDs within a record must not create a self-conflict when
    # they happen to straddle two write batches.
    if np.any(atoms[1:] <= atoms[:-1]):
        atoms = np.unique(atoms)
    for frames in presence.selections():
        frame_count = (
            frames.stop - frames.start if isinstance(frames, slice) else len(frames)
        )
        atoms_per_batch = max(1, _MATRIX_WRITE_CELLS // frame_count)
        for start in range(0, len(atoms), atoms_per_batch):
            batch = atoms[start : start + atoms_per_batch]
            index = (
                (batch, frames) if isinstance(frames, slice) else np.ix_(batch, frames)
            )
            overlap = np.not_equal(store.atomeach[index], 0)
            if np.any(overlap):
                # The packed conflict API takes indices. Materialize them only
                # for an overlapping bounded span, never for a full trajectory.
                conflict_frames = (
                    np.arange(frames.start, frames.stop)
                    if isinstance(frames, slice)
                    else frames
                )
                store.conflict.mark(batch, conflict_frames, overlap)
            store.atomeach[index] = molecule_id


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


def _count_route_pairs(pairs):
    """Count a bounded block, packing the common narrow unsigned IDs exactly."""
    if pairs.dtype.kind == "u" and pairs.dtype.itemsize <= 4:
        # Shared matrices use unsigned IDs. Packing two uint32 fields makes
        # repetitive routes one integer sort instead of a two-field row sort.
        keys = pairs[:, 0].astype(np.uint64)
        np.left_shift(keys, np.uint64(32), out=keys)
        np.bitwise_or(keys, pairs[:, 1], out=keys)
        unique, counts = np.unique(keys, return_counts=True)
        values = np.column_stack(
            (unique >> np.uint64(32), unique & np.uint64(2**32 - 1))
        )
        return values.astype(pairs.dtype), counts
    # Wider and signed arrays retain both original fields without truncation.
    return np.unique(pairs, axis=0, return_counts=True)


def _compact_route_pairs(pairs):
    """Count repetitive pairs with bounded scratch space and a unique-pair cap.

    Short or apparently unique routes keep their original array. The sample
    only selects an optimization; either representation has identical counts.
    Stop aggregating if the dictionary would grow beyond one scan block.
    """
    if len(pairs) < _ROUTE_COMPACTION_MIN_ROWS:
        return pairs, None
    stride = max(1, len(pairs) // _ROUTE_COMPACTION_SAMPLE_ROWS)
    sample = pairs[::stride][:_ROUTE_COMPACTION_SAMPLE_ROWS]
    if len(_count_route_pairs(sample)[0]) == len(sample):
        return pairs, None
    counts = Counter()
    for start in range(0, len(pairs), _STEP3_SCAN_ROWS):
        unique, frequencies = _count_route_pairs(
            pairs[start : start + _STEP3_SCAN_ROWS]
        )
        for pair, frequency in zip(unique, frequencies, strict=True):
            counts[int(pair[0]), int(pair[1])] += int(frequency)
            if len(counts) > _ROUTE_COMPACTION_MAX_PAIRS:
                return pairs, None
    ordered = sorted(counts)
    return (
        np.asarray(ordered, dtype=pairs.dtype).reshape((-1, 2)),
        np.fromiter((counts[pair] for pair in ordered), dtype=np.int64),
    )


class _SpeciesRouteAccumulator:
    """Reduce ordered atom results without retaining event-level pair arrays."""

    def __init__(self, names, run_hmm):
        self.names = names
        self.run_hmm = run_hmm
        self.counts = Counter()
        # Python integer keys avoid fixed-width overflow and tuple/NumPy scalar
        # overhead. Ownership must survive across blocks and across atoms.
        self._span = len(names) + 1
        self._owners = {}

    def add(self, atom_index, pairs, frequencies):
        """Keep every occurrence on the first atom, or one globally with HMM."""
        if frequencies is not None and len(frequencies) != len(pairs):
            raise ValueError("Route pair counts do not match route pairs")
        for start in range(0, len(pairs), _STEP3_SCAN_ROWS):
            block = pairs[start : start + _STEP3_SCAN_ROWS]
            selected = []
            for index, pair in enumerate(block):
                key = int(pair[0]) * self._span + int(pair[1])
                owner = self._owners.get(key)
                if owner is None:
                    self._owners[key] = atom_index
                elif self.run_hmm or owner != atom_index:
                    continue
                selected.append(index)
            if not selected:
                continue
            # Resolve names in bounded NumPy batches. Scalar name lookups for
            # every unique event would penalize routes that cannot be compacted.
            species = self.names[block[selected] - 1]
            different = species[:, 0] != species[:, 1]
            if frequencies is None or self.run_hmm:
                self.counts.update(map(tuple, species[different]))
            else:
                weights = frequencies[start : start + len(block)][selected]
                for pair, count in zip(
                    species[different], weights[different], strict=True
                ):
                    self.counts[tuple(pair)] += int(count)


def _compact_atom_route_result(*args):
    """Share the same counted representation between local scans and workers."""
    pairs, text = _atom_route_result(*args)
    pairs, frequencies = _compact_route_pairs(pairs)
    return pairs, frequencies, text


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
    return _compact_atom_route_result(
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
        context = (
            _TimedOutputWriter(self.rng) if self.rng.timed_output else nullcontext()
        )
        with context as writer:
            self._collect(writer)

    def _collect(self, writer):
        """Keep the optional writer alive until all reaction results are consumed."""
        self.atomnames = self.atomname[self.atomtype]
        self._printmoleculename()
        if writer is not None:
            writer.write_molecules(self.mname)
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
                timed_writer=writer,
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
            with closing(_iter_molecule_presence(self)) as records:
                for molecule_id, (record, presence) in enumerate(
                    tqdm(
                        records,
                        total=self.hmmit,
                        desc="Analyze atoms",
                        unit="molecule",
                        disable=None,
                    ),
                    start=1,
                ):
                    _write_molecule_presence(
                        store, bytestolist(record[0]), presence, molecule_id
                    )
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

    def _iter_local_atom_routes(self, matrix_store, frame_range, mark_active):
        """Read small routes directly without changing parent mapping ownership."""
        start, stop = frame_range
        active = matrix_store.active_transitions if mark_active else None
        for atom_index in range(self.N):
            name = self.atomname[self.atomtype[atom_index]]
            yield _compact_atom_route_result(
                atom_index,
                matrix_store.atomeach[atom_index, start:stop],
                name,
                name in self.selectatoms,
                self.mname,
                active,
            )

    def _printatomroute(self, matrix_store, timeaxis=None, frame_range=None):
        """Write atom routes and reduce their pairs to species-level counts."""
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
            maxbytes=_ROUTE_WRITE_BUFFER_BYTES,
        ) as f:
            routes = _SpeciesRouteAccumulator(self.mname, self.runHMM)
            scan_size = self.N * (frame_range[1] - frame_range[0])
            if scan_size <= _STEP3_SERIAL_ATOM_FRAMES:
                results = self._iter_local_atom_routes(
                    matrix_store, frame_range, mark_active=timeaxis is None
                )
            else:
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
                        matrix_store.active_transition_path
                        if timeaxis is None
                        else None,
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
                for ii, (pairs, frequencies, routestr) in enumerate(results):
                    f.append(routestr)
                    routes.add(ii, pairs, frequencies)
        return routes.counts

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
        m = species_molecule(self.atomnames[atoms], atoms, bonds)
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
                        if isomer.isomorphic(molecule):
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


def _isomer_candidate_key(graph):
    """Bucket candidates by a necessary invariant of the policy-adjusted graph.

    Isomorphisms preserve node labels and multisets of labelled neighbors.
    Equal keys are only candidates: regular graphs and hash collisions still
    require the unchanged VF2 comparison. The process-local key is never a
    species identity. Retain one integer per bucket instead of graph-sized
    signature objects; temporary storage is proportional to one molecule.
    """
    neighborhoods = Counter(
        (
            data.get("atom"),
            frozenset(
                Counter(
                    (graph.nodes[neighbor].get("atom"), edge.get("level"))
                    for neighbor, edge in graph[node].items()
                ).items()
            ),
        )
        for node, data in graph.nodes(data=True)
    )
    return hash(frozenset(neighborhoods.items()))


class _IsomerCandidateBuckets:
    """Index formulas with multiple structural groups, keeping cheap cases local."""

    def __init__(self):
        """Keep initial representatives without computing graph keys."""
        self.first = {}
        self.indexed = None

    def candidates(self, molecule):
        """Return possible matches in their original insertion order."""
        if self.indexed is None:
            # Composition-only graphs cannot be separated by neighborhoods.
            # Repeated variants of one topology still need its VF2 comparison;
            # indexing them cannot reduce the number of candidate groups.
            if molecule._miso == 2 or len(self.first) < 2:
                return self.first
            self.indexed = defaultdict(dict)
            for name, first in self.first.items():
                self.indexed[_isomer_candidate_key(first.graph)][name] = first
            self.first.clear()
        return self.indexed[_isomer_candidate_key(molecule.graph)]


class _CollectSMILESPaths(_CollectPaths):
    def _printmoleculename(self):
        if self.miso > 0:
            self._printmergedmoleculename()
            return
        mname = _MoleculeNameBuilder(self.hmmit)
        d = defaultdict(list)
        self.n_unknown = 0
        timeline = (
            self._openmoleculetimelinespool() if self._needmoleculetimeline() else None
        )
        try:
            with (
                WriteBuffer(open(self.moleculefilename, "w"), sep="\n") as fm,
                open(self.moleculetemp2filename, "rb") as ft,
            ):
                with closing(
                    self._getSMILESresults(ft, self._calmoleculeSMILESname)
                ) as results:
                    for name, atoms, bonds, frames in results:
                        name = self._resolveSMILESname(name, atoms, bonds, d)
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

    def _printmergedmoleculename(self):
        """Select representatives before writing any merged-name output.

        Frequency is the number of detected molecule-frame occurrences stored
        before HMM filtering. Keep only compact name IDs and one graph per
        candidate in memory, then stream the original records a second time so
        every output uses the final representative.
        """
        source_names = _MoleculeNameBuilder(self.hmmit)
        name_frequency = Counter()
        name_group = {}
        name_mapping_graph = defaultdict(_IsomerCandidateBuckets)
        d = defaultdict(list)
        self.n_unknown = 0
        with (
            open(self.moleculetemp2filename, "rb") as ft,
            closing(
                self._getSMILESresults(ft, self._calmoleculeSMILESfrequency)
            ) as results,
        ):
            for name, atoms, bonds, frequency in results:
                name = str(self._resolveSMILESname(name, atoms, bonds, d))
                source_names.append(name)
                name_frequency[name] += frequency
                if name in name_group:
                    continue
                molecule = Molecule(self, atoms, bonds)
                # The graph already applies miso's bond/connectivity policy.
                # A bucket only excludes impossible matches; VF2 decides all
                # remaining candidates in their original insertion order.
                groups = name_mapping_graph[molecule.name].candidates(molecule)
                for group_name, group_molecule in groups.items():
                    if group_molecule.isomorphic(molecule):
                        name_group[name] = group_name
                        break
                else:
                    groups[name] = molecule
                    name_group[name] = name
        source_names = source_names.finish()

        group_members = defaultdict(list)
        for name, group_name in name_group.items():
            group_members[group_name].append(name)
        representative = {}
        for members in group_members.values():
            selected = min(members, key=lambda name: (-name_frequency[name], name))
            representative.update(dict.fromkeys(members, selected))

        mname = _MoleculeNameBuilder(self.hmmit)
        timeline = (
            self._openmoleculetimelinespool() if self._needmoleculetimeline() else None
        )
        try:
            with (
                WriteBuffer(open(self.moleculefilename, "w"), sep="\n") as fm,
                open(self.moleculetemp2filename, "rb") as ft,
            ):
                records = itertools.zip_longest(*[read_compressed_block(ft)] * 4)
                for source_name, line in zip(source_names, records, strict=True):
                    name = representative[str(source_name)]
                    atoms, bonds = self._getatomsandbonds(line)
                    frames, _ = self._getmoleculeframesandtimesteps(
                        line, need_timesteps=False
                    )
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

    def _getSMILESresults(self, molecule_file, func):
        """Convert four-block records lazily, preserving molecule order."""
        if self.nproc == 1:
            # Keep conversion in the caller without accumulating decoded graphs.
            records = itertools.zip_longest(*[read_compressed_block(molecule_file)] * 4)
            return tqdm(
                map(func, records),
                total=self.hmmit,
                desc="Indentify isomers",
                unit="molecule",
                disable=None,
            )
        if (
            _should_batch_smiles(self.hmmit, self.nproc)
            and type(self) is _CollectSMILESPaths
            and getattr(func, "__func__", None)
            in (
                _CollectSMILESPaths._calmoleculeSMILESname,
                _CollectSMILESPaths._calmoleculeSMILESfrequency,
            )
            and not any(callable(value) for value in vars(self).values())
        ):
            # Custom collectors/callables keep their bound-method behavior.
            return _run_smiles_batches(
                self, molecule_file, func.__name__ == "_calmoleculeSMILESfrequency"
            )
        return run_mp(
            self.nproc,
            func=func,
            l=read_compressed_block(molecule_file),
            unordered=False,
            nlines=4,
            chunksize=1,
            max_inflight=max(2, 2 * self.nproc),
            disk_ordered=True,
            total=self.hmmit,
            desc="Indentify isomers",
            unit="molecule",
        )

    def _resolveSMILESname(self, name, atoms, bonds, molecules):
        """Return a full-structure name, using VF2 when SMILES conversion fails."""
        if name is not None:
            return name
        # A canonical SMILES identifies the structure before ``miso`` grouping.
        # Preserve the same boundary in the fallback so distinct candidates do
        # not pool their frequencies before a representative is selected.
        molecule = Molecule(self, atoms, bonds, miso=0)

        # Avoid repeating the failed and potentially expensive RDKit call.
        def _raise_anyway(*args, **kwargs):
            raise ValueError("Maximum BFS search size exceeded.")

        molecule._convertSMILES = _raise_anyway
        atom_match = iso.categorical_node_match("atom", None)
        bond_match = iso.categorical_edge_match("level", None)
        for isomer in molecules[str(molecule)]:
            if nx.is_isomorphic(
                isomer.graph,
                molecule.graph,
                node_match=atom_match,
                edge_match=bond_match,
            ):
                molecule.smiles = isomer.smiles
                break
        else:
            molecules[str(molecule)].append(molecule)
        return molecule.smiles

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

    def _calmoleculeSMILESfrequency(self, item):
        """Return a candidate name and its pre-HMM occurrence count."""
        line = item
        atoms, bonds = self._getatomsandbonds(line)
        frequency = len(bytestolist(line[-1]))
        try:
            name = self.convertSMILES(atoms, bonds)
        except ValueError:
            name = None
        return name, atoms, bonds, frequency


class Molecule:
    """A molecule graph represented at a selected isomer-identity level."""

    def __init__(self, cmp, atoms, bonds, *, miso=None):
        self.cmp = cmp
        self.atoms = atoms
        self.bonds = bonds
        self._atomtypes = cmp.atomtype[atoms]
        self._atomnames = cmp.atomnames[atoms]
        self._miso = cmp.miso if miso is None else miso
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
        """Build a VF2 graph according to the selected isomer equivalence mode.

        Normal mode shares aromatic bond normalization with SMILES conversion
        and retains raw bond levels if RDKit cannot construct the molecule.
        The broader modes discard bond orders or connectivity, respectively.
        Trajectory atom indices and the original bond records are preserved.
        """
        if self._miso not in (0, 1, 2):
            raise ValueError(f"Unknown isomer identification method: {self._miso}.")
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
        for atomnumber, atomtype in zip(self.atoms, self._atomtypes):
            graph.add_node(atomnumber, atom=atomtype)
        if self._miso == 0:
            try:
                molecule = species_molecule(self._atomnames, self.atoms, self.bonds)
            except (ValueError, RuntimeError):
                # VF2 also supports graphs RDKit cannot construct. Keep the
                # existing raw-graph fallback when even construction fails.
                return graph
            for bond in molecule.GetBonds():
                atom1 = self.atoms[bond.GetBeginAtomIdx()]
                atom2 = self.atoms[bond.GetEndAtomIdx()]
                graph[atom1][atom2]["level"] = int(bond.GetBondType())
        return graph

    def isomorphic(self, mol):
        """Return whether two molecules are isomorphic."""
        return nx.is_isomorphic(
            self.graph,
            mol.graph,
            node_match=_ATOM_MATCH,
            edge_match=_BOND_MATCH,
        )
