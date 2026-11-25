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
import itertools
import operator
import os
import pickle
import tempfile
from abc import ABCMeta, abstractmethod
from array import array
from collections import Counter, defaultdict
from contextlib import closing, suppress
from enum import Enum, auto
from multiprocessing import get_start_method
from typing import ClassVar

import numpy as np

try:
    from openbabel import __version__ as obversion
    from openbabel import openbabel
except ImportError:  # pragma: no cover
    raise ImportError("Open Babel 3.1.0 is required.")
from ase import Atoms
from ase.data import atomic_numbers, covalent_radii
from ase.neighborlist import natural_cutoffs, neighbor_list
from packaging import version
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree  # ty: ignore[unresolved-import]
from tqdm.auto import tqdm

from .dps import dps  # type:ignore
from .utils import SharedRNGData, WriteBuffer, listtobytes, run_mp

if version.parse(obversion) < version.parse("3.1.0"):  # pragma: no cover
    raise ImportError("Open Babel 3.1.0 is required.")

openbabel.obErrorLog.StopLogging()

_OPENBABEL_PERIODIC_NEIGHBOR_MIN_ATOMS = 2048
_OPENBABEL_COVALENT_RADIUS_TOLERANCE = 0.45
_OPENBABEL_MIN_BOND_DISTANCE = 0.4
_OPENBABEL_MIN_BOND_ANGLE = 45.0
_OPENBABEL_DISTANCE_AMBIGUITY = 1e-7
_HYSTERESIS_KDTREE_MARGIN = 1e-9
_DETECT_HYSTERESIS_RECORD_CACHE_MAX = 16384
_DETECT_WRITE_BYTES = 1024 * 1024
_DETECT_COMPRESS_MIN_RECORDS = 24576
_DETECT_COMPRESS_MAX_FRAME_BYTES = 1024
_DETECT_COMPRESS_RAW_BYTES = 64 * 1024
_DETECT_COMPRESS_INPUT_BYTES = 256 * 1024
_DETECT_COMPRESS_RESULT_BYTES = 1024 * 1024
_DETECT_COMPRESS_MAX_RECORDS = 2048
_DETECT_COMPRESS_MAX_INFLIGHT = 4


def _compression_batch_size(values, nproc):
    """Select bounded batches only for many short records under fork.

    Spawn startup and small/dense workloads retain the serial path. A single
    large record also selects that path, so no oversized task enters the pool.
    The largest raw record determines a conservative byte-bounded batch size;
    the separate item cap bounds array/pickle metadata even for empty records.
    """
    if (
        nproc < 2
        or len(values) < _DETECT_COMPRESS_MIN_RECORDS
        or get_start_method() != "fork"
    ):
        return 0
    largest = max((len(v) * v.itemsize for v in values), default=0)
    if largest > _DETECT_COMPRESS_MAX_FRAME_BYTES or not values:
        return 0
    return min(
        _DETECT_COMPRESS_MAX_RECORDS,
        _DETECT_COMPRESS_RAW_BYTES // max(1, largest),
    )


def _frame_batches(values, batch_size):
    """Serialize one planned batch, checking wire bytes before task submission.

    Inputs contain trusted, local unsigned arrays. Their raw storage and record
    count are bounded by the plan; this additional check covers pickle metadata
    rather than treating raw array bytes as the complete IPC payload.
    """
    source = iter(values)
    while batch := list(itertools.islice(source, batch_size)):
        payload = pickle.dumps(batch, protocol=pickle.HIGHEST_PROTOCOL)
        if len(payload) > _DETECT_COMPRESS_INPUT_BYTES:
            raise ValueError("Detection compression input exceeds its byte budget")
        yield payload


def _compress_frame_batch(payload):
    """Compress a bounded batch and cap serialized results before queueing.

    The internal byte envelope keeps both directions of run_mp IPC measurable.
    Each contained record retains the existing independent compressed format.
    """
    if len(payload) > _DETECT_COMPRESS_INPUT_BYTES:
        raise ValueError("Detection compression input exceeds its byte budget")
    # Exact histories often repeat across distinct molecule definitions. Reuse
    # only their immutable compressed bytes, retaining one output per record.
    # This cache dies with the batch: <=64 KiB of raw keys, <=2048 entries and
    # <=1 MiB of unique compressed values, plus bounded container overhead.
    cache = {}
    blocks = []
    stored_bytes = 0
    for value in pickle.loads(payload):
        key = (value.typecode, value.tobytes())
        block = cache.get(key)
        if block is None:
            block = listtobytes(np.asarray(value))
            stored_bytes += len(block)
            if stored_bytes > _DETECT_COMPRESS_RESULT_BYTES:
                raise ValueError("Detection compression result exceeds its byte budget")
            cache[key] = block
        blocks.append(block)
    result = pickle.dumps(blocks, protocol=pickle.HIGHEST_PROTOCOL)
    if len(result) > _DETECT_COMPRESS_RESULT_BYTES:
        raise ValueError("Detection compression result exceeds its byte budget")
    return result


def _widen_frame_indices(values, frame):
    """Append an overflowing frame without losing earlier order or duplicates.

    Array widths are platform-dependent. Choose by actual itemsize and skip
    intermediate widths that cannot hold a sparse jump in the frame index.
    The original array stays intact if no unsigned representation is available.
    """
    frame = operator.index(frame)
    for typecode in ("H", "I", "Q"):
        width = array(typecode).itemsize
        if width > values.itemsize and 0 <= frame < 1 << (8 * width):
            widened = array(typecode, values)
            widened.append(frame)
            return widened
    raise OverflowError(f"Frame index {frame} exceeds available unsigned storage")


class _Detect(SharedRNGData, metaclass=ABCMeta):
    """Detect molecules.

    Parameters
    ----------
    rng: reacnetgenerator.ReacNetGenerator
        The ReacNetGenerator class.
    """

    subclasses: ClassVar[dict] = {}
    # Set only for ordered, opt-in detection. Edges use stable zero-based atom IDs.
    _active_bonds: set[tuple[int, int]] | None = None

    # type hints
    # TODO: we need a better way to communicate with parameters
    inputfilename: list
    atomname: np.ndarray
    stepinterval: int
    nproc: int
    pbc: bool
    cell: np.ndarray
    max_component_atoms: int
    max_component_fraction: float
    N: int
    atomtype: np.ndarray

    def __init__(self, rng):
        SharedRNGData.__init__(
            self,
            rng,
            [
                "inputfilename",
                "atomname",
                "stepinterval",
                "nproc",
                "pbc",
                "cell",
                "max_component_atoms",
                "max_component_fraction",
            ],
            ["N", "atomtype", "step", "timestep", "temp1it", "moleculetempfilename"],
        )
        self._prev_graph_pair: tuple[list, list] | None = None
        self._prev_molecules: list[bytes] | None = None
        self._record_bytes_cache: (
            dict[tuple[tuple[int, ...], tuple[tuple[int, int, int], ...]], bytes] | None
        ) = None

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

    def _clear_hysteresis_state(self) -> None:
        """Reset transient per-run hysteresis caches."""
        self._active_bonds = None
        self._prev_graph_pair = None
        self._prev_molecules = None
        self._record_bytes_cache = None
        if hasattr(self, "_active_edge_keys"):
            self._active_edge_keys = None
        if hasattr(self, "_prev_hysteresis_bonds"):
            self._prev_hysteresis_bonds = None
        if hasattr(self, "_prev_retained_orders"):
            self._prev_retained_orders = None
        if hasattr(self, "_prev_bond_graph"):
            self._prev_bond_graph = None

    def _readinputfile(self):
        """Read the input file."""
        # Detection results may arrive out of frame order. Keep that order and
        # duplicate observations while avoiding Python references per visit.
        d = defaultdict(lambda: array("H"))
        timestep = {}
        with fileinput.input(files=self.inputfilename) as f:
            _steplinenum = self._readNfunc(f)
        hysteresis = getattr(self.rng, "bond_hysteresis", None) is not None
        self._clear_hysteresis_state()
        if hysteresis:
            self._active_bonds = set()
            self._record_bytes_cache = {}
        source_stream = None
        try:
            if hysteresis:
                # A retained pair depends on the preceding analyzed frame. Keep
                # this pass ordered across input files; later stages can use workers.
                use_chunk_stream = (
                    type(self)._source_lines is _Detect._source_lines
                    and "_source_lines" not in self.__dict__
                )
                if use_chunk_stream:
                    source_stream = self._source_frames(_steplinenum)
                    sampled = itertools.islice(
                        source_stream, 0, None, self.stepinterval
                    )
                else:
                    source_stream = self._source_lines(_steplinenum)
                    source_frames = itertools.zip_longest(
                        *[iter(source_stream)] * _steplinenum
                    )
                    sampled = itertools.islice(
                        source_frames, 0, None, self.stepinterval
                    )
                results = (self._readstepfunc(item) for item in enumerate(sampled))
            else:
                source_stream = (
                    self._source_lines(_steplinenum)
                    if self.rng.timed_output
                    else fileinput.FileInput(files=self.inputfilename)
                )
                results = run_mp(
                    self.nproc,
                    func=self._readstepfunc,
                    l=source_stream,
                    max_inflight=max(1, self.nproc * 2),
                    nlines=_steplinenum,
                    return_num=True,
                    interval=self.stepinterval,
                    desc="Read bond information and Detect molecules",
                    unit="timestep",
                )
            last_molecules = None
            last_arrays = None
            with closing(results):
                for molecules, (step, thetimestep) in results:
                    if (
                        hysteresis
                        and molecules is last_molecules
                        and last_arrays is not None
                    ):
                        appended = 0
                        try:
                            for frames in last_arrays:
                                frames.append(step)
                                appended += 1
                        except OverflowError:
                            for frames in last_arrays[:appended]:
                                frames.pop()
                            for molecule in molecules:
                                frames = d[molecule]
                                try:
                                    frames.append(step)
                                except OverflowError:
                                    d[molecule] = _widen_frame_indices(frames, step)
                            last_arrays = [d[molecule] for molecule in molecules]
                    else:
                        for molecule in molecules:
                            frames = d[molecule]
                            try:
                                frames.append(step)
                            except OverflowError:
                                d[molecule] = _widen_frame_indices(frames, step)
                        if hysteresis:
                            last_molecules = molecules
                            last_arrays = [d[molecule] for molecule in molecules]
                    timestep[step] = thetimestep
        finally:
            if source_stream is not None:
                source_stream.close()
            self._clear_hysteresis_state()
        self.temp1it = len(d)
        # Keep the cheap serial path, but amortize parallel-compression IPC for
        # many short fork records. Explicit closure also cancels suspended worker
        # results when compression or the bounded writer fails.
        with closing(self._compressed_values(d.values())) as compressed:
            with tqdm(
                compressed,
                desc="Save molecules",
                unit="molecule",
                total=self.temp1it,
                disable=None,
            ) as values_c:
                self._writemoleculetempfile((d.keys(), values_c))
        self.timestep = timestep
        self.step = len(timestep)

    @staticmethod
    def _open_source_file(filename):
        """Open a trajectory input file directly unless stdin is requested."""
        if str(filename) == "-":
            return fileinput.FileInput(files=(filename,))
        return open(filename)

    def _source_frames(self, lines_per_frame):
        """Stream fixed-size frame line tuples while recording source frame counts."""
        counts = []
        for source_index, filename in enumerate(self.inputfilename):
            frame_count = 0
            with self._open_source_file(filename) as source:
                while chunk := tuple(itertools.islice(source, lines_per_frame)):
                    if len(chunk) != lines_per_frame:
                        if source_index + 1 == len(self.inputfilename):
                            raise ValueError("Incomplete final source frame")
                        raise ValueError("Incomplete source frame at file boundary")
                    frame_count += 1
                    yield chunk
            counts.append(frame_count)
        self.rng.source_frame_counts = counts

    def _source_lines(self, lines_per_frame):
        """Record per-source frame counts during the existing streaming read.

        Counts are indexed by input occurrence, so repeating the same filename
        remains distinguishable. Reject partial frames at file boundaries rather
        than silently inventing a source mapping for a cross-file frame.
        """
        counts = []
        for source_index, filename in enumerate(self.inputfilename):
            lines = 0
            with self._open_source_file(filename) as source:
                for line in source:
                    lines += 1
                    yield line
            if lines % lines_per_frame:
                if source_index + 1 == len(self.inputfilename):
                    raise ValueError("Incomplete final source frame")
                raise ValueError("Incomplete source frame at file boundary")
            counts.append(lines // lines_per_frame)
        if len(counts) != len(self.inputfilename):
            raise ValueError("Empty input sources cannot be mapped to frames")
        self.rng.source_frame_counts = counts

    def _compressed_values(self, values):
        """Yield frames in input order with a fixed byte-bounded task window.

        At most four submitted-but-unconsumed batches carry 256 KiB input and
        1 MiB result envelopes each. Ordered run_mp spools early results to disk
        and does not replenish the window past its oldest unconsumed batch.
        Only the current decoded batch and the 1 MiB writer buffer are needed
        for delivery. Resident detection records and process overhead are separate.
        """
        batch_size = _compression_batch_size(values, self.nproc)
        if not batch_size:
            yield from map(self._compressvalue, values)
            return
        with closing(
            run_mp(
                min(self.nproc, _DETECT_COMPRESS_MAX_INFLIGHT),
                func=_compress_frame_batch,
                l=_frame_batches(values, batch_size),
                max_inflight=_DETECT_COMPRESS_MAX_INFLIGHT,
                disk_ordered=True,
                unordered=False,
                total=(len(values) + batch_size - 1) // batch_size,
                bar=False,
                desc="Save molecules",
            )
        ) as results:
            for payload in results:
                yield from pickle.loads(payload)

    @staticmethod
    def _compressvalue(x):
        """Serialize the compact frame buffer without an intermediate array copy."""
        return listtobytes(np.asarray(x))

    @abstractmethod
    def _readNfunc(self, f) -> int:
        pass

    @abstractmethod
    def _readstepfunc(self, item) -> tuple[list[bytes], tuple[int, int]]:
        pass

    def _validate_component_sizes(
        self, mols: list[list[int]], frame=None, timestep=None
    ):
        if self.max_component_atoms == 0:
            return
        limit = max(
            self.max_component_atoms,
            int(np.ceil(self.N * self.max_component_fraction)),
        )
        if not mols:
            return
        largest = max(mols, key=list.__len__)
        if len(largest) <= limit:
            return

        composition = Counter(
            str(self.atomname[self.atomtype[atom]]) for atom in largest
        )
        composition_text = ", ".join(
            f"{element}:{count}" for element, count in sorted(composition.items())
        )
        location = f"frame {frame}"
        if timestep is not None:
            location += f" (timestep {timestep})"

        if self._active_bonds is not None:
            backend = "bond hysteresis"
            advice = (
                f"Inspect the configured bond thresholds: {self.rng.bond_hysteresis}."
            )
        elif isinstance(self, _DetectCrd):
            if self.use_ase:
                backend = "ASE"
                cutoffs = self.custom_cutoffs or "none"
                advice = (
                    "ASE natural cutoffs can include non-covalent contacts. "
                    "Inspect the connectivity and use "
                    "--ase-pair-cutoffs El1-El2:0 to exclude an unintended "
                    f"element pair. Current ASE pair cutoffs: {cutoffs}."
                )
            else:
                backend = "Open Babel"
                advice = (
                    "For coordinate trajectories, consider --use-ase together "
                    "with --ase-pair-cutoffs El1-El2:0 to control connectivity."
                )
        else:
            backend = "explicit bond input"
            advice = (
                "For explicit bond trajectories, revise the upstream bond-order "
                "cutoff or bonding rule."
            )

        raise RuntimeError(
            f"Oversized connected component detected at {location}: "
            f"{len(largest)} of {self.N} atoms exceed the configured limit of "
            f"{limit} (max_component_atoms={self.max_component_atoms}, "
            f"max_component_fraction={self.max_component_fraction}). "
            f"Composition: {composition_text}. Bond detection backend: {backend}. "
            f"{advice} If this component is intentional, increase the component "
            "limits or set --max-component-atoms 0 to disable this guard."
        )

    def _connectmolecule(self, bond, level, frame=None, timestep=None):
        if self._active_bonds is not None:
            if (
                self._prev_graph_pair is not None
                and bond is self._prev_graph_pair[0]
                and level is self._prev_graph_pair[1]
                and self._prev_molecules is not None
            ):
                return self._prev_molecules
            mols, bondlists = dps(bond, level)
            self._validate_component_sizes(mols, frame=frame, timestep=timestep)
            cache = self._record_bytes_cache
            if cache is None:
                cache = {}
                self._record_bytes_cache = cache
            records = []
            for mol, bondlist in zip(mols, bondlists):
                key = (tuple(mol), tuple(bondlist))
                encoded = cache.get(key)
                if encoded is None:
                    encoded = b"".join(
                        (
                            listtobytes(mol),
                            listtobytes([(int(i[0]), int(i[1])) for i in bondlist]),
                            listtobytes([int(i[2]) for i in bondlist]),
                        )
                    )
                    if len(cache) >= _DETECT_HYSTERESIS_RECORD_CACHE_MAX:
                        cache.pop(next(iter(cache)))
                    cache[key] = encoded
                records.append(encoded)
            self._prev_graph_pair = (bond, level)
            self._prev_molecules = records
            return records
        mols, bondlists = dps(bond, level)
        self._validate_component_sizes(mols, frame=frame, timestep=timestep)
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
        """Stream aligned definitions/frames, publishing only a complete file.

        The outer file context also closes the handle if a buffer flush fails.
        Remove partial output after closing it, including on Windows, while
        leaving any previously published filename untouched.
        """
        handle = tempfile.NamedTemporaryFile("wb", delete=False)
        try:
            with handle, WriteBuffer(handle, maxbytes=_DETECT_WRITE_BYTES) as output:
                for mol in zip(*d, strict=True):
                    output.extend(mol)
        except BaseException:
            with suppress(FileNotFoundError):
                os.unlink(handle.name)
            raise
        self.moleculetempfilename = handle.name


def _bond_lists(n_atoms, edges):
    """Expand undirected, unique (pair, order) entries for molecule detection."""
    bonds = [[] for _ in range(n_atoms)]
    orders = [[] for _ in range(n_atoms)]
    for (i, j), order in edges:
        bonds[i].append(j)
        bonds[j].append(i)
        orders[i].append(order)
        orders[j].append(order)
    return bonds, orders


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

    def _readstepfunc(self, item) -> tuple[list[bytes], tuple[int, int]]:
        step, lines = item
        bond: list[tuple[int, ...] | None] = [None] * self.N
        level: list[tuple[int, ...] | None] = [None] * self.N
        timestep = None
        active_bonds = self._active_bonds
        reported = {} if active_bonds is not None else None
        for line in lines:
            if line:
                if line[0] == "#":
                    if line.startswith("# Timestep"):
                        timestep = int(line.split()[-1])
                else:
                    s = line.split()
                    s0 = int(s[0]) - 1
                    s2 = int(s[2])
                    if reported is None:
                        bond[s0] = tuple(map(self._get_idx, s[3 : 3 + s2]))
                        level[s0] = tuple(map(self._get_bo, s[4 + s2 : 4 + 2 * s2]))
                    else:
                        neighbors = tuple(map(self._get_idx, s[3 : 3 + s2]))
                        orders = s[4 + s2 : 4 + 2 * s2]
                        for j, raw in zip(neighbors, orders, strict=True):
                            if not (0 <= s0 < self.N and 0 <= j < self.N and s0 != j):
                                raise ValueError(
                                    "Invalid reported ReaxFF bond endpoints"
                                )
                            order = float(raw)
                            if not np.isfinite(order) or order < 0:
                                raise ValueError("Invalid reported ReaxFF bond order")
                            pair = (min(s0, j), max(s0, j))
                            reported[pair] = max(order, reported.get(pair, 0.0))
        assert timestep is not None
        if reported is not None:
            assert active_bonds is not None
            thresholds = self.rng.bond_hysteresis
            assert thresholds is not None
            form_threshold = thresholds["form"]
            break_threshold = thresholds["break"]
            retained = {
                pair: max(1, round(order))
                for pair, order in reported.items()
                if order >= form_threshold
                or (pair in active_bonds and order > break_threshold)
            }
            self._active_bonds = set(retained)
            prev_retained = getattr(self, "_prev_retained_orders", None)
            prev_graph = getattr(self, "_prev_bond_graph", None)
            if prev_retained == retained and prev_graph is not None:
                bond, level = prev_graph
            else:
                bond, level = _bond_lists(self.N, retained.items())
                self._prev_retained_orders = retained
                self._prev_bond_graph = (bond, level)
        molecules = self._connectmolecule(bond, level, frame=step, timestep=timestep)
        return molecules, (step, timestep)

    @staticmethod
    def _get_idx(x):
        return int(x) - 1

    @staticmethod
    def _get_bo(x):
        return max(1, round(float(x)))


class _DetectCrd(_Detect):
    use_ase: bool
    ase_cutoff_mult: float
    custom_cutoffs: str | None

    def __init__(self, rng):
        SharedRNGData.__init__(
            self,
            rng,
            [
                "inputfilename",
                "atomname",
                "stepinterval",
                "nproc",
                "pbc",
                "cell",
                "use_ase",
                "ase_cutoff_mult",
                "custom_cutoffs",
                "max_component_atoms",
                "max_component_fraction",
            ],
            ["N", "atomtype", "step", "timestep", "temp1it", "moleculetempfilename"],
        )
        self._prev_graph_pair: tuple[list, list] | None = None
        self._prev_molecules: list[bytes] | None = None
        self._record_bytes_cache: (
            dict[tuple[tuple[int, ...], tuple[tuple[int, int, int], ...]], bytes] | None
        ) = None
        self._parsed_custom_cutoffs = self._parse_custom_cutoffs(self.custom_cutoffs)
        self._pair_break = None
        self._atom_outer = None
        self._active_edge_keys: np.ndarray | None = None
        self._prev_hysteresis_bonds: tuple[list[list[int]], list[list[int]]] | None = (
            None
        )
        self._cached_hysteresis_numbers: np.ndarray | None = None
        self._cached_hysteresis_radii: np.ndarray | None = None
        self._cached_hysteresis_max_cutoff: float = 0.0
        self._numbers_by_type = Atoms(symbols=list(self.atomname)).numbers
        thresholds = getattr(rng, "bond_hysteresis", None)
        if thresholds and "pair_break" in thresholds:
            # Index by atomic number in either endpoint order, independently of
            # the dump type mapping and the direction returned by ASE.
            n = len(covalent_radii)
            self._pair_break = np.full((n, n), thresholds["break"], dtype=float)
            for pair, value in thresholds["pair_break"].items():
                i, j = (atomic_numbers[symbol] for symbol in pair.split("-"))
                self._pair_break[i, j] = self._pair_break[j, i] = value
            # Each element's search radius only needs to cover its widest pair.
            self._atom_outer = self._pair_break.max(axis=1)

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
        step_atoms.set_cell(cell)
        step_atoms.set_pbc(self.pbc)
        if self.pbc:
            step_atoms.wrap()

        # Build per-atom natural cutoff radii (with global multiplier)
        radii = natural_cutoffs(step_atoms, mult=self.ase_cutoff_mult)
        symbols = step_atoms.get_chemical_symbols()

        # Build a per-symbol radius lookup (same symbol always has same radius)
        radii_by_symbol: dict = {}
        for sym, r in zip(symbols, radii):
            radii_by_symbol[sym] = r

        # Build complete per-pair cutoff dict for all unique element pairs.
        # ASE neighbor_list dict format: {('El1', 'El2'): distance}
        # Pairs NOT in the dict get cutoff=0, so we must cover ALL pairs.
        unique_symbols = sorted(set(symbols))
        cutoff_dict: dict = {}
        for s1 in unique_symbols:
            for s2 in unique_symbols:
                pair = frozenset({s1, s2})
                if pair in self._parsed_custom_cutoffs:
                    cutoff_dict[(s1, s2)] = self._parsed_custom_cutoffs[pair]
                else:
                    cutoff_dict[(s1, s2)] = radii_by_symbol[s1] + radii_by_symbol[s2]

        i_list, j_list = neighbor_list("ij", step_atoms, cutoff=cutoff_dict)

        atomnumber = len(step_atoms)
        bond: list = [[] for _ in range(atomnumber)]
        bondlevel: list = [[] for _ in range(atomnumber)]

        for i, j in zip(i_list, j_list):
            if i < j and j not in bond[i]:
                bond[i].append(j)
                bond[j].append(i)
                bondlevel[i].append(1)
                bondlevel[j].append(1)

        return bond, bondlevel

    def _new_openbabel_molecule(
        self,
        step_atoms: Atoms,
        cell: np.ndarray,
    ):
        """Create an Open Babel molecule while preserving input atom IDs."""
        mol = openbabel.OBMol()
        mol.BeginModify()
        for idx, (num, position) in enumerate(
            zip(step_atoms.get_atomic_numbers(), step_atoms.positions, strict=True)
        ):
            atom = mol.NewAtom(idx)
            atom.SetAtomicNum(int(num))
            atom.SetVector(*position)
        if self.pbc:
            unit_cell = openbabel.OBUnitCell()
            unit_cell.SetData(
                openbabel.vector3(*cell[0]),
                openbabel.vector3(*cell[1]),
                openbabel.vector3(*cell[2]),
            )
            mol.CloneData(unit_cell)
            mol.SetPeriodicMol()
        return mol

    @staticmethod
    def _getbondlistsfromopenbabel(mol, atomnumber):
        """Convert one perceived Open Babel molecule to adjacency lists."""
        bond = [[] for _ in range(atomnumber)]
        bondlevel = [[] for _ in range(atomnumber)]
        for obbond in openbabel.OBMolBondIter(mol):
            atom1 = obbond.GetBeginAtom().GetId()
            atom2 = obbond.GetEndAtom().GetId()
            level = obbond.GetBondOrder()
            if level == 5:
                # aromatic, 5 in openbabel but 12 in rdkit
                level = 12
            bond[atom1].append(atom2)
            bond[atom2].append(atom1)
            bondlevel[atom1].append(level)
            bondlevel[atom2].append(level)
        return bond, bondlevel

    def _getbondfromopenbabel(
        self,
        step_atoms: Atoms,
        cell: np.ndarray,
    ) -> tuple[list[list[int]], list[list[int]]]:
        """Run the original Open Babel coordinate-bond perception path."""
        mol = self._new_openbabel_molecule(step_atoms, cell)
        mol.ConnectTheDots()
        mol.PerceiveBondOrders()
        mol.EndModify()
        return self._getbondlistsfromopenbabel(mol, len(step_atoms))

    @staticmethod
    def _orthorhombic_cell_lengths(cell):
        """Return positive axis-aligned cell lengths when safely supported."""
        cell = np.asarray(cell, dtype=np.float64)
        if cell.shape != (3, 3) or not np.all(np.isfinite(cell)):
            return None
        lengths = np.diag(cell)
        if np.any(lengths <= 0.0) or not np.allclose(
            cell,
            np.diag(lengths),
            rtol=0.0,
            atol=1e-12,
        ):
            return None
        return lengths

    def _getperiodicbondcandidates(
        self,
        step_atoms: Atoms,
        cell: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray] | None:
        """Return Open Babel-ordered periodic bond candidates when safe."""
        atomnumber = len(step_atoms)
        positions = np.asarray(step_atoms.positions, dtype=np.float64)
        atomic_numbers = np.asarray(step_atoms.get_atomic_numbers())
        if (
            positions.shape != (atomnumber, 3)
            or not np.all(np.isfinite(positions))
            or atomic_numbers.shape != (atomnumber,)
        ):
            return None
        radii = np.fromiter(
            (
                openbabel.GetCovalentRad(int(atomic_number))
                for atomic_number in atomic_numbers
            ),
            dtype=np.float64,
            count=atomnumber,
        )
        if not np.all(np.isfinite(radii)) or np.any(radii <= 0.0):
            return None

        cell_lengths = self._orthorhombic_cell_lengths(cell)
        maximum_pair_cutoff = (
            2.0 * float(np.max(radii)) + _OPENBABEL_COVALENT_RADIUS_TOLERANCE
        )
        if cell_lengths is None:
            return None
        if np.min(cell_lengths) <= 2.0 * (
            maximum_pair_cutoff + _OPENBABEL_DISTANCE_AMBIGUITY
        ):
            # Avoid duplicate lattice images and unbounded neighbor-list growth
            # in cells narrower than twice the largest possible bond distance.
            return None

        z_order = np.argsort(positions[:, 2], kind="stable")
        sorted_z = positions[z_order, 2]
        has_z_ties = np.any(sorted_z[1:] == sorted_z[:-1])

        wrapped_positions = np.mod(positions, cell_lengths)
        # Floating-point modulo can round a tiny negative coordinate to the
        # excluded upper edge of cKDTree's periodic box.
        if np.any(wrapped_positions >= cell_lengths):
            return None
        candidate_pairs = cKDTree(
            wrapped_positions,
            boxsize=cell_lengths,
        ).query_pairs(
            maximum_pair_cutoff + _OPENBABEL_DISTANCE_AMBIGUITY,
            output_type="ndarray",
        )
        atom_indices = candidate_pairs[:, 0]
        neighbor_indices = candidate_pairs[:, 1]
        distance_vectors = (
            wrapped_positions[neighbor_indices] - wrapped_positions[atom_indices]
        )
        distance_vectors -= np.rint(distance_vectors / cell_lengths) * cell_lengths
        distances = np.sqrt(np.einsum("ij,ij->i", distance_vectors, distance_vectors))

        pair_cutoffs = (
            radii[atom_indices]
            + radii[neighbor_indices]
            + _OPENBABEL_COVALENT_RADIUS_TOLERANCE
        )
        boundary_ambiguous = (
            np.abs(distances - pair_cutoffs) <= _OPENBABEL_DISTANCE_AMBIGUITY
        ) | (
            np.abs(distances - _OPENBABEL_MIN_BOND_DISTANCE)
            <= _OPENBABEL_DISTANCE_AMBIGUITY
        )
        if np.any(boundary_ambiguous):
            return None

        accepted = (distances <= pair_cutoffs) & (
            distances >= _OPENBABEL_MIN_BOND_DISTANCE
        )
        atom_indices = atom_indices[accepted]
        neighbor_indices = neighbor_indices[accepted]

        if has_z_ties:
            # Open Babel leaves equal-z ordering to the C++ standard library.
            # Ties in different candidate components cannot change the relative
            # insertion order inside either component. Connectivity cleanup
            # only removes candidate edges, so those components cannot join.
            # Keep the reference fallback for ties connected even indirectly.
            graph = csr_matrix(
                (
                    np.ones(len(atom_indices), dtype=np.bool_),
                    (atom_indices, neighbor_indices),
                ),
                shape=(atomnumber, atomnumber),
            )
            _, labels = connected_components(graph, directed=False)
            component_order = np.lexsort((positions[:, 2], labels))
            component_labels = labels[component_order]
            component_z = positions[component_order, 2]
            if np.any(
                (component_labels[1:] == component_labels[:-1])
                & (component_z[1:] == component_z[:-1])
            ):
                return None

        z_rank = np.empty(atomnumber, dtype=np.intp)
        z_rank[z_order] = np.arange(atomnumber, dtype=np.intp)
        reverse_pair = z_rank[atom_indices] > z_rank[neighbor_indices]
        first_atoms = np.where(reverse_pair, neighbor_indices, atom_indices)
        second_atoms = np.where(reverse_pair, atom_indices, neighbor_indices)
        insertion_order = np.lexsort((z_rank[second_atoms], z_rank[first_atoms]))
        return first_atoms[insertion_order], second_atoms[insertion_order]

    @staticmethod
    def _add_openbabel_candidate_bonds(mol, first_atoms, second_atoms):
        """Insert candidates and preserve Open Babel's connectivity cleanup."""
        for atom_index, neighbor_index in zip(first_atoms, second_atoms, strict=True):
            atom = mol.GetAtom(int(atom_index) + 1)
            neighbor = mol.GetAtom(int(neighbor_index) + 1)
            if (
                openbabel.GetMaxBonds(atom.GetAtomicNum()) == 0
                or openbabel.GetMaxBonds(neighbor.GetAtomicNum()) == 0
            ):
                continue
            if atom.GetAtomicNum() == 15 and atom.GetExplicitValence() == 5:
                if neighbor.GetAtomicNum() not in (9, 17):
                    continue
            if neighbor.GetAtomicNum() == 15 and neighbor.GetExplicitValence() == 5:
                if atom.GetAtomicNum() not in (9, 17):
                    continue
            mol.AddBond(int(atom_index) + 1, int(neighbor_index) + 1, 1)

        for atom in openbabel.OBMolAtomIter(mol):
            while (
                atom.GetExplicitValence() > openbabel.GetMaxBonds(atom.GetAtomicNum())
                or atom.SmallestBondAngle() < _OPENBABEL_MIN_BOND_ANGLE
            ):
                bonds = list(openbabel.OBAtomBondIter(atom))
                if not bonds:
                    break
                if atom.GetAtomicNum() == 1:
                    hydrogen_bond = next(
                        (
                            bond
                            for bond in bonds
                            if bond.GetNbrAtom(atom).GetAtomicNum() == 1
                        ),
                        None,
                    )
                    if hydrogen_bond is not None:
                        mol.DeleteBond(hydrogen_bond)
                        continue
                longest_bond = bonds[0]
                longest_length = longest_bond.GetLength()
                for bond in bonds[1:]:
                    length = bond.GetLength()
                    if length > longest_length:
                        longest_bond = bond
                        longest_length = length
                mol.DeleteBond(longest_bond)

    def _getbondfromperiodicneighborlist(
        self,
        step_atoms: Atoms,
        cell: np.ndarray,
    ) -> tuple[list[list[int]], list[list[int]]] | None:
        """Accelerate periodic Open Babel connectivity with bounded candidates.

        A periodic cKDTree supplies nearby atom pairs for axis-aligned cells,
        then Open Babel retains responsibility for insertion rules,
        valence/angle cleanup, and bond-order perception. Unsupported or
        numerically ambiguous geometries return ``None`` so the caller can use
        ``ConnectTheDots`` without changing its historical result.
        """
        atomnumber = len(step_atoms)
        if atomnumber < 2:
            return ([[] for _ in range(atomnumber)], [[] for _ in range(atomnumber)])

        candidates = self._getperiodicbondcandidates(step_atoms, cell)
        if candidates is None:
            return None

        mol = self._new_openbabel_molecule(step_atoms, cell)
        self._add_openbabel_candidate_bonds(mol, *candidates)

        mol.PerceiveBondOrders()
        mol.EndModify()
        return self._getbondlistsfromopenbabel(mol, atomnumber)

    def _hysteresis_radii_and_cutoff(
        self, numbers: np.ndarray
    ) -> tuple[np.ndarray, float] | None:
        """Return per-atom covalent radii and tightest search cutoff."""
        if (
            self._cached_hysteresis_numbers is not None
            and np.array_equal(numbers, self._cached_hysteresis_numbers)
            and self._cached_hysteresis_radii is not None
        ):
            return self._cached_hysteresis_radii, self._cached_hysteresis_max_cutoff
        if (
            numbers.ndim != 1
            or not len(numbers)
            or np.any((numbers <= 0) | (numbers >= len(covalent_radii)))
        ):
            return None
        radii = covalent_radii[numbers]
        unique_nums = np.unique(numbers)
        unique_r = covalent_radii[unique_nums]
        if self._pair_break is None:
            thresholds = self.rng.bond_hysteresis
            assert thresholds is not None
            max_cutoff = float(2.0 * np.max(unique_r) * thresholds["break"])
        else:
            limits_u = self._pair_break[np.ix_(unique_nums, unique_nums)]
            max_cutoff = float(
                np.max((unique_r[:, None] + unique_r[None, :]) * limits_u)
            )
        self._cached_hysteresis_numbers = numbers.copy()
        self._cached_hysteresis_radii = radii
        self._cached_hysteresis_max_cutoff = max_cutoff
        return radii, max_cutoff

    def _get_kdtree_hysteresis_pairs(
        self,
        numbers: np.ndarray,
        positions: np.ndarray,
        cell: np.ndarray,
        max_cutoff: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
        """Return candidate pairs (i < j) and exact distances via cKDTree when safe."""
        n_atoms = len(numbers)
        if n_atoms < 2:
            return (
                np.empty(0, dtype=np.intp),
                np.empty(0, dtype=np.intp),
                np.empty(0, dtype=np.float64),
            )
        if positions.shape != (n_atoms, 3) or not np.all(np.isfinite(positions)):
            return None
        search_radius = max_cutoff + _HYSTERESIS_KDTREE_MARGIN
        if self.pbc:
            cell_lengths = self._orthorhombic_cell_lengths(cell)
            if cell_lengths is None or np.min(cell_lengths) <= 2.0 * (
                max_cutoff + _OPENBABEL_DISTANCE_AMBIGUITY
            ):
                return None
            wrapped = np.mod(positions, cell_lengths)
            # Negative tiny coordinates can modulo-round to the exact box edge,
            # which cKDTree rejects as out-of-bounds.
            wrapped[wrapped >= cell_lengths] = 0.0
            pairs = cKDTree(wrapped, boxsize=cell_lengths).query_pairs(
                search_radius, output_type="ndarray"
            )
            if not len(pairs):
                return (
                    np.empty(0, dtype=np.intp),
                    np.empty(0, dtype=np.intp),
                    np.empty(0, dtype=np.float64),
                )
            first = pairs[:, 0]
            second = pairs[:, 1]
            dv = positions[second] - positions[first]
            dv -= np.rint(dv / cell_lengths) * cell_lengths
        else:
            pairs = cKDTree(positions).query_pairs(search_radius, output_type="ndarray")
            if not len(pairs):
                return (
                    np.empty(0, dtype=np.intp),
                    np.empty(0, dtype=np.intp),
                    np.empty(0, dtype=np.float64),
                )
            first = pairs[:, 0]
            second = pairs[:, 1]
            dv = positions[second] - positions[first]
        distances = np.sqrt(np.einsum("ij,ij->i", dv, dv))
        return first, second, distances

    def _get_hysteresis_bonds(self, atoms, cell):
        """Use the closest periodic image and assign orders on the retained graph."""
        previous = self._active_bonds
        assert previous is not None
        thresholds = self.rng.bond_hysteresis
        assert thresholds is not None
        form_mult = thresholds["form"]
        break_mult = thresholds["break"]
        n_atoms = len(atoms)
        numbers = np.asarray(atoms.numbers)
        positions = np.asarray(atoms.positions, dtype=np.float64)

        if self._active_edge_keys is None:
            prev_keys = (
                np.array(
                    sorted(min(i, j) * n_atoms + max(i, j) for i, j in previous),
                    dtype=np.int64,
                )
                if previous
                else np.empty(0, dtype=np.int64)
            )
        else:
            prev_keys = self._active_edge_keys

        radii_and_cutoff = self._hysteresis_radii_and_cutoff(numbers)
        kdtree_result = (
            self._get_kdtree_hysteresis_pairs(
                numbers, positions, cell, radii_and_cutoff[1]
            )
            if radii_and_cutoff is not None
            else None
        )
        if kdtree_result is not None:
            assert radii_and_cutoff is not None
            radii = radii_and_cutoff[0]
            first, second, distances = kdtree_result
        else:
            atoms.set_cell(cell)
            atoms.set_pbc(self.pbc)
            radii = covalent_radii[numbers]
            # ASE returns both directions and can return multiple periodic images.
            # Existence of any qualifying image is equivalent to using the closest.
            # Widen search radii only for elements participating in wider overrides.
            cutoffs = (
                radii * break_mult
                if self._atom_outer is None
                else radii * self._atom_outer[numbers]
            )
            first, second, distances = neighbor_list("ijd", atoms, cutoffs)
            mask = first < second
            first = first[mask]
            second = second[mask]
            distances = distances[mask]

        if len(first):
            pair_radii = radii[first] + radii[second]
            form_mask = distances <= pair_radii * form_mult
            keys_form = first[form_mask].astype(np.int64) * n_atoms + second[form_mask]
            if prev_keys.size:
                limits = (
                    break_mult
                    if self._pair_break is None
                    else self._pair_break[numbers[first], numbers[second]]
                )
                hold_mask = (~form_mask) & (distances < pair_radii * limits)
                if np.any(hold_mask):
                    keys_cand = (
                        first[hold_mask].astype(np.int64) * n_atoms + second[hold_mask]
                    )
                    pos_idx = np.searchsorted(prev_keys, keys_cand)
                    in_prev = (pos_idx < prev_keys.size) & (
                        prev_keys[np.minimum(pos_idx, prev_keys.size - 1)] == keys_cand
                    )
                    if np.any(in_prev):
                        retained_keys = np.unique(
                            np.concatenate((keys_form, keys_cand[in_prev]))
                        )
                    else:
                        retained_keys = np.unique(keys_form)
                else:
                    retained_keys = np.unique(keys_form)
            else:
                retained_keys = np.unique(keys_form)
        else:
            retained_keys = np.empty(0, dtype=np.int64)

        unchanged = self._active_edge_keys is not None and np.array_equal(
            retained_keys, prev_keys
        )
        self._active_edge_keys = retained_keys

        if self.use_ase:
            if unchanged and self._prev_hysteresis_bonds is not None:
                return self._prev_hysteresis_bonds
            u_arr, v_arr = np.divmod(retained_keys, n_atoms)
            bonds = [[] for _ in range(n_atoms)]
            orders = [[] for _ in range(n_atoms)]
            for i, j in zip(u_arr.tolist(), v_arr.tolist(), strict=True):
                bonds[i].append(j)
                bonds[j].append(i)
                orders[i].append(1)
                orders[j].append(1)
            result = (bonds, orders)
            self._prev_hysteresis_bonds = result
            return result

        u_arr, v_arr = np.divmod(retained_keys, n_atoms)
        pairs = list(zip(u_arr.tolist(), v_arr.tolist(), strict=True))
        mol = openbabel.OBMol()
        mol.BeginModify()
        mol.ReserveAtoms(n_atoms)
        nums_list = numbers.tolist()
        pos_list = positions.tolist()
        new_atom = mol.NewAtom
        for idx in range(n_atoms):
            atom = new_atom(idx)
            atom.SetAtomicNum(int(nums_list[idx]))
            px, py, pz = pos_list[idx]
            atom.SetVector(px, py, pz)
        if self.pbc:
            unit_cell = openbabel.OBUnitCell()
            unit_cell.SetData(
                openbabel.vector3(*cell[0]),
                openbabel.vector3(*cell[1]),
                openbabel.vector3(*cell[2]),
            )
            mol.CloneData(unit_cell)
            mol.SetPeriodicMol()
        add_bond = mol.AddBond
        for i, j in pairs:
            add_bond(i + 1, j + 1, 1)
        mol.PerceiveBondOrders()
        mol.EndModify()
        assert mol.NumBonds() == len(pairs)
        bond = [[] for _ in range(n_atoms)]
        bondlevel = [[] for _ in range(n_atoms)]
        get_bond = mol.GetBond
        for bond_idx, (u, v) in enumerate(pairs):
            raw_level = get_bond(bond_idx).GetBondOrder()
            level = 12 if raw_level == 5 else raw_level
            bond[u].append(v)
            bond[v].append(u)
            bondlevel[u].append(level)
            bondlevel[v].append(level)
        return bond, bondlevel

    def _getbondfromcrd(
        self, step_atoms: Atoms, cell: np.ndarray
    ) -> tuple[list[list[int]], list[list[int]]]:
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
        if self._active_bonds is not None:
            return self._get_hysteresis_bonds(step_atoms, cell)

        # Check if ASE mode is enabled
        if self.use_ase:
            return self._getbondfromase(step_atoms, cell)

        atomnumber = len(step_atoms)
        if (
            self.pbc
            and atomnumber >= _OPENBABEL_PERIODIC_NEIGHBOR_MIN_ATOMS
            and self._orthorhombic_cell_lengths(cell) is not None
        ):
            accelerated = self._getbondfromperiodicneighborlist(step_atoms, cell)
            if accelerated is not None:
                return accelerated
        return self._getbondfromopenbabel(step_atoms, cell)


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
        iscompleted = False
        N = None
        linecontent = None
        atomtype = None
        seen_atom_ids = None
        index = -1
        stepaindex = None
        for index, line in enumerate(f):
            if line.startswith("ITEM:"):
                linecontent = self.LineType.linecontent(line)
                if linecontent == self.LineType.ATOMS:
                    keys = line.split()
                    self._n_atom_cols = len(keys) - 2
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
                    atomtype = np.empty(N, dtype=int)
                    seen_atom_ids = np.zeros(N, dtype=bool)
                elif linecontent == self.LineType.ATOMS:
                    s = line.split()
                    assert atomtype is not None and seen_atom_ids is not None
                    atom_index, atom_type = self._atom_indices(s, seen_atom_ids)
                    atomtype[atom_index] = atom_type
        else:
            steplinenum = index + 1
        assert N is not None and atomtype is not None
        assert seen_atom_ids is not None
        if not np.all(seen_atom_ids):
            raise ValueError("LAMMPS dump frame is missing one or more atom IDs")
        self.N = N
        self.atomtype = atomtype
        return steplinenum

    def _atom_indices(self, fields, seen_atom_ids):
        """Validate the existing contiguous, one-based dump ID/type contract.

        Mark each ID once so direct array placement cannot silently overwrite an
        atom or leave uninitialized positions, including during first-frame setup.
        """
        atom_index = int(fields[self.id_idx]) - 1
        atom_type = int(fields[self.tidx]) - 1
        if not 0 <= atom_index < len(seen_atom_ids):
            raise ValueError("LAMMPS dump atom ID is out of range")
        if seen_atom_ids[atom_index]:
            raise ValueError("LAMMPS dump contains a duplicate atom ID")
        if not 0 <= atom_type < len(self.atomname):
            raise ValueError("LAMMPS dump atom type is out of range")
        seen_atom_ids[atom_index] = True
        return atom_index, atom_type

    def _try_fast_parse_frame(
        self, step: int, lines
    ) -> tuple[Atoms, np.ndarray, tuple[int, int]] | None:
        """Vectorized parser for standard numeric LAMMPS dump frames."""
        n_cols = getattr(self, "_n_atom_cols", 0)
        if (
            not isinstance(lines, (tuple, list))
            or len(lines) != self.N + 9
            or n_cols <= 0
            or not lines[0].startswith("ITEM: TIMESTEP")
            or not lines[2].startswith("ITEM: NUMBER OF ATOMS")
            or not lines[4].startswith("ITEM: BOX")
            or not lines[8].startswith("ITEM: ATOMS")
        ):
            return None
        try:
            ts = int(lines[1].split()[0])
            b0 = [float(x) for x in lines[5].split()]
            b1 = [float(x) for x in lines[6].split()]
            b2 = [float(x) for x in lines[7].split()]
        except (IndexError, ValueError):
            return None
        if len(b0) == 2 and len(b1) == 2 and len(b2) == 2:
            boxsize = np.array(
                [
                    [b0[1] - b0[0], 0.0, 0.0],
                    [0.0, b1[1] - b1[0], 0.0],
                    [0.0, 0.0, b2[1] - b2[0]],
                ]
            )
        elif (
            len(b0) >= 3
            and len(b1) >= 3
            and len(b2) >= 3
            and len(b0) == len(b1) == len(b2)
        ):
            xy, xz, yz = b0[2], b1[2], b2[2]
            xlo = b0[0] - min(0.0, xy, xz, xy + xz)
            xhi = b0[1] - max(0.0, xy, xz, xy + xz)
            ylo = b1[0] - min(0.0, yz)
            yhi = b1[1] - max(0.0, yz)
            zlo = b2[0]
            zhi = b2[1]
            boxsize = np.array(
                [
                    [xhi - xlo, 0.0, 0.0],
                    [xy, yhi - ylo, 0.0],
                    [xz, yz, zhi - zlo],
                ]
            )
        else:
            return None
        try:
            raw = np.fromstring("".join(lines[9:]), dtype=np.float64, sep=" ")
        except ValueError:
            return None
        if raw.size != self.N * n_cols:
            return None
        table = raw.reshape((self.N, n_cols))
        raw_id = table[:, self.id_idx]
        raw_type = table[:, self.tidx]
        if not (np.all(np.isfinite(raw_id)) and np.all(np.isfinite(raw_type))):
            return None
        atom_index = raw_id.astype(np.intp) - 1
        atom_type = raw_type.astype(np.intp) - 1
        numbers_by_type = self._numbers_by_type
        if (
            np.any(atom_index + 1 != raw_id)
            or np.any((atom_index < 0) | (atom_index >= self.N))
            or np.any(atom_type + 1 != raw_type)
            or np.any((atom_type < 0) | (atom_type >= len(numbers_by_type)))
        ):
            return None
        if self.N == 0 or (
            atom_index[0] == 0
            and atom_index[-1] == self.N - 1
            and np.all(atom_index[1:] - atom_index[:-1] == 1)
        ):
            step_numbers = numbers_by_type[atom_type]
            if self.yidx == self.xidx + 1 and self.zidx == self.yidx + 1:
                step_positions = table[:, self.xidx : self.zidx + 1].copy()
            else:
                step_positions = table[:, (self.xidx, self.yidx, self.zidx)].copy()
        else:
            seen = np.zeros(self.N, dtype=bool)
            seen[atom_index] = True
            if not np.all(seen):
                return None
            step_numbers = np.empty(self.N, dtype=numbers_by_type.dtype)
            step_positions = np.empty((self.N, 3), dtype=np.float64)
            step_numbers[atom_index] = numbers_by_type[atom_type]
            step_positions[atom_index, 0] = table[:, self.xidx]
            step_positions[atom_index, 1] = table[:, self.yidx]
            step_positions[atom_index, 2] = table[:, self.zidx]
        return (
            Atoms(numbers=step_numbers, positions=step_positions),
            boxsize,
            (step, ts),
        )

    def _readstepfunc(self, item) -> tuple[list[bytes], tuple[int, int]]:
        step, lines = item
        parsed = self._try_fast_parse_frame(step, lines)
        if parsed is not None:
            step_atoms, boxsize, timestep = parsed
        else:
            # Place rows in canonical atom order without per-atom ASE objects or sorting.
            numbers_by_type = self._numbers_by_type
            step_numbers = np.empty(self.N, dtype=numbers_by_type.dtype)
            step_positions = np.empty((self.N, 3), dtype=float)
            seen_atom_ids = np.zeros(self.N, dtype=bool)
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
                            atom_index, atom_type = self._atom_indices(s, seen_atom_ids)
                            step_numbers[atom_index] = numbers_by_type[atom_type]
                            step_positions[atom_index] = (
                                float(s[self.xidx]),
                                float(s[self.yidx]),
                                float(s[self.zidx]),
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
            if not np.all(seen_atom_ids):
                raise ValueError("LAMMPS dump frame is missing one or more atom IDs")
            step_atoms = Atoms(numbers=step_numbers, positions=step_positions)
        bond, level = self._getbondfromcrd(step_atoms, boxsize)
        molecules = self._connectmolecule(
            bond, level, frame=timestep[0], timestep=timestep[1]
        )
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

    def _readstepfunc(self, item) -> tuple[list[bytes], tuple[int, int]]:
        step, lines = item
        timestep = step, step
        boxsize = self.cell
        if self.pbc and boxsize is None:
            raise RuntimeError("No cell information is given.")
        line_seq = lines if isinstance(lines, (tuple, list)) else tuple(lines)
        n_atoms = max(0, len(line_seq) - 2)
        step_numbers = np.empty(n_atoms, dtype=int)
        step_positions = np.empty((n_atoms, 3), dtype=float)
        for index, line in enumerate(line_seq[2:]):
            s = line.split()
            step_numbers[index] = atomic_numbers[s[0]]
            step_positions[index] = (float(s[1]), float(s[2]), float(s[3]))
        step_atoms = Atoms(numbers=step_numbers, positions=step_positions)
        bond, level = self._getbondfromcrd(step_atoms, boxsize)
        molecules = self._connectmolecule(
            bond, level, frame=timestep[0], timestep=timestep[1]
        )
        return molecules, timestep


@_Detect.register_subclass("extxyz")
class _Detectextxyz(_DetectCrd):
    """extxyz file. xyz with extended metadata support like cell, force, etc."""

    @staticmethod
    def _parse_lattice(comment: str) -> np.ndarray:
        """Parse the 3x3 lattice stored in an extxyz comment line."""
        try:
            lattice_str = comment.split("Lattice=", 1)[1].split('"', 2)[1]
            lattice = np.array(list(map(float, lattice_str.split())))
            return lattice.reshape((3, 3))
        except (IndexError, ValueError) as exc:
            raise RuntimeError(
                "Missing or invalid Lattice= information in extxyz."
            ) from exc

    def _readNfunc(self, f):
        atomname_dict = dict(zip(self.atomname.tolist(), range(self.atomname.size)))
        N = None
        atomtype = None
        for index, line in enumerate(f):
            if index == 0:
                N = int(line.strip())
                atomtype = np.zeros(N, dtype=int)
            elif index == 1:
                self.cell = self._parse_lattice(line)
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

    def _readstepfunc(self, item) -> tuple[list[bytes], tuple[int, int]]:
        step, lines = item
        timestep = step, step  # Use step as timestep fallback
        line_seq = lines if isinstance(lines, (tuple, list)) else tuple(lines)
        if len(line_seq) < 2:
            raise RuntimeError("No cell information is given in extxyz.")
        # Extxyz trajectories may change their cell between frames,
        # so bond detection must use this frame's metadata.
        boxsize = self._parse_lattice(line_seq[1])
        n_atoms = len(line_seq) - 2
        step_numbers = np.empty(n_atoms, dtype=int)
        step_positions = np.empty((n_atoms, 3), dtype=float)
        for index, line in enumerate(line_seq[2:]):
            s = line.split()
            step_numbers[index] = atomic_numbers[s[0]]
            step_positions[index] = (float(s[1]), float(s[2]), float(s[3]))
        step_atoms = Atoms(numbers=step_numbers, positions=step_positions)
        bond, level = self._getbondfromcrd(step_atoms, boxsize)
        molecules = self._connectmolecule(
            bond, level, frame=timestep[0], timestep=timestep[1]
        )
        return molecules, timestep
