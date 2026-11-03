# SPDX-License-Identifier: LGPL-3.0-or-later
# cython: language_level=3
# cython: linetrace=True
"""HMM Filter.

In order to filter noise, a two-state HMM was adopted, which can be
described as a transition matrix A, an emission matrix B, and an initial
state vector π. The existence of molecules can be converted into 0-1 signals.
In order to predict the state sequence according to the output sequence,
Viterbi Algorithm is used to acquire the path with the most likely hidden
state sequence V called Viterbi path.

References
----------
.. [1] Wang, L.-P.; McGibbon, R. T.; Pande, V. S.; Martínez, T. J. Automated
   discovery and refinement of reactive molecular dynamics pathways. J. Chem.
   Theory Comput. 2016, 12(2), 638-649.
.. [2] Rabiner, L. R. A trtorial on hidden Markov models and selected
   applications in speech recognition. Proc. IEEE 1989, 77(2), 257-286.
.. [3] Forney, G. D. The viterbi algorithm. Porc. IEEE 1973, 61(3), 268-278.
"""

import itertools
import os
import tempfile
from contextlib import ExitStack, closing, suppress
from functools import partial
from multiprocessing import get_start_method
from typing import Any

import numpy as np
from tqdm.auto import tqdm

try:
    # hmmlearn v0.2.8 renamed MultinomialHMM to CategoricalHMM
    from hmmlearn.hmm import CategoricalHMM as MultinomialHMM
except ImportError:
    from hmmlearn.hmm import MultinomialHMM

from ._moleculeorder import _OrderedMoleculeResults
from ._moleculepresence import _iter_molecule_records, _use_frame_indices
from ._step3state import _unsigned_dtype_for_maximum
from .utils import (
    SharedRNGData,
    WriteBuffer,
    appendIfNotNone,
    bytestolist,
    check_zero_signal,
    idx_to_signal,
    listtobytes,
    read_compressed_block,
    run_mp,
)

_HMM_BOUNDED_SIGNAL_FRAMES = 100_000
_HMM_LEGACY_CHUNK_RECORDS = 100
_HMM_WORKER_FUNCTION = None


def _initialize_hmm_worker(step, model, getoriginfile, printfiltersignal):
    """Retain prediction state without serializing the parent RNG context."""
    global _HMM_WORKER_FUNCTION
    _HMM_WORKER_FUNCTION = partial(
        _convert_hmm_record,
        step=step,
        model=model,
        getoriginfile=getoriginfile,
        printfiltersignal=printfiltersignal,
    )


def _convert_hmm_record(item, *, step, model, getoriginfile, printfiltersignal):
    """Expand and filter one record using only prediction state."""
    value = bytestolist(item[-1])
    origin = idx_to_signal(value, step)
    originbytes = listtobytes(origin) if getoriginfile else None
    hmmbytes = None
    if model is not None:
        hmmsignal = model.predict(origin).astype(bool)
        # Origin output retains metadata even for filtered-out molecules.
        if check_zero_signal(hmmsignal) or printfiltersignal or getoriginfile:
            hmmbytes = listtobytes(hmmsignal)
    return originbytes, hmmbytes, item


def _variant_groups(records, atoms, types, miso):
    """Group records on the same atoms using the naming stage's VF2 policy."""
    if miso == 2 or all(record[1] == records[0][1] for record in records):
        return [range(len(records))]

    from ._species import isomorphic, molecule_graph

    groups = []
    representatives = []
    by_endpoints = {}
    for i, record in enumerate(records):
        endpoints = record[1]
        if endpoints not in by_endpoints:
            graph = molecule_graph(atoms, types, bytestolist(endpoints), miso)
            for group_id, representative in enumerate(representatives):
                if isomorphic(graph, representative):
                    break
            else:
                group_id = len(groups)
                groups.append([])
                representatives.append(graph)
            by_endpoints[endpoints] = group_id
        groups[by_endpoints[endpoints]].append(i)
    return groups


def _filter_atom_set(task, *, miso, step, model, getoriginfile, printfiltersignal):
    """Filter each union once, retaining aligned raw definitions and observations."""
    atoms, types, records = task
    # Slots preserve canonical structure order when VF2 classes interleave.
    output = [None] * len(records)
    for indices in _variant_groups(records, atoms, types, miso):
        group_id = indices[0]
        if len(indices) == 1:
            origin, signal, record = _convert_hmm_record(
                records[group_id],
                step=step,
                model=model,
                getoriginfile=getoriginfile,
                printfiltersignal=printfiltersignal,
            )
            if signal is not None:
                output[group_id] = (group_id, origin, signal, record)
            continue
        # A detected component has exactly one raw structure per observed frame.
        owner = np.full(step, -1, dtype=np.intp)
        for i in indices:
            frames = bytestolist(records[i][-1])
            owner[frames] = i
        observed = owner != -1
        kept = model.predict(observed.astype(np.int8).reshape(-1, 1)).astype(bool)
        if not (kept.any() or getoriginfile or printfiltersignal):
            continue
        # Project the pooled path back to disjoint raw ranges. At an inferred
        # frame use the preceding observation, or the first for a leading gap.
        gaps = np.flatnonzero(kept & ~observed)
        if gaps.size:
            visits = np.flatnonzero(observed)
            owner[gaps] = (
                owner[visits[np.maximum(visits.searchsorted(gaps) - 1, 0)]]
                if visits.size
                else group_id
            )
        for i in indices:
            record = records[i]
            origin = (
                listtobytes(idx_to_signal(bytestolist(record[-1]), step))
                if getoriginfile
                else None
            )
            signal = listtobytes(kept & (owner == i))
            output[i] = (group_id, origin, signal, record)
    return [result for result in output if result is not None]


def _filter_hmm_record(record):
    """Filter one record without resending the RNG context with every task."""
    assert _HMM_WORKER_FUNCTION is not None
    return _HMM_WORKER_FUNCTION(record)


def _expand_nohmm_origin(record, *, frame_count):
    """Expand a dense record without sending trajectory metadata to a worker."""
    frames = record[0]
    if _use_frame_indices(frames, frame_count):
        return None
    return listtobytes(idx_to_signal(bytestolist(frames), frame_count))


class _HMMFilter(SharedRNGData):
    HMM_DESCRIPTION = "HMM filter"
    NOHMM_DESCRIPTION = "Pass through original signals (HMM disabled)"

    miso: int = 0
    runHMM: bool
    getoriginfile: bool
    printfiltersignal: bool
    moleculetempfilename: str
    nproc: int
    temp1it: int
    p: np.ndarray
    a: np.ndarray
    b: np.ndarray
    step: int

    def __init__(self, rng):
        SharedRNGData.__init__(
            self,
            rng,
            [
                "runHMM",
                "miso",
                "getoriginfile",
                "printfiltersignal",
                "moleculetempfilename",
                "nproc",
                "temp1it",
                "p",
                "a",
                "b",
                "step",
            ],
            [
                "moleculetemp2filename",
                "originfilename",
                "hmmfilename",
                "hmmit",
                "hmm_route_ids",
            ],
        )

    @property
    def description(self) -> str:
        """Return the description for the active signal-processing mode."""
        if self.runHMM:
            return self.HMM_DESCRIPTION
        return self.NOHMM_DESCRIPTION

    def filter(self):
        """HMM Filters.

        Timesteps of molecules are converted to a visible output sequence.
        O^m=(o_t^m) is given by o_t^m={1, if m exists; 0, otherwise}.
        Similarly, a hidden state sequence I^m=(i_t^m) is given by
        i_t^m={1, if m exists; 0, otherwise.}
        """
        if self.runHMM:
            self._initHMM()
        self._calhmm()
        self.returnkeys()

    def _initHMM(self):
        self._model = MultinomialHMM(n_components=2, algorithm="viterbi")
        self._model.startprob_ = self.p
        self._model.transmat_ = self.a
        self._model.emissionprob_ = self.b

    def _getoriginandhmm(self, item):
        return _convert_hmm_record(
            item,
            step=self.step,
            model=self._model if self.runHMM else None,
            getoriginfile=self.getoriginfile,
            printfiltersignal=self.printfiltersignal,
        )

    def _cal_nohmm(self):
        """Reuse definitions and store only dense fallback signals, in order."""
        self.moleculetemp2filename = self.moleculetempfilename
        self.hmmfilename = None
        self.hmmit = 0
        with (
            open(self.moleculetempfilename, "rb") as definitions,
            WriteBuffer(
                tempfile.NamedTemporaryFile("wb", delete=False), linenumber=0
            ) as origin,
            tqdm(
                total=self.temp1it,
                desc=self.description,
                unit="molecule",
                disable=None,
            ) as progress,
        ):
            self.originfilename = origin.name
            records = _iter_molecule_records(definitions, (3,))
            for record in records:
                if _use_frame_indices(record[0], self.step):
                    self.hmmit += 1
                    progress.update()
                    continue

                # All-sparse inputs never start workers. Once expansion is
                # needed, keep the remaining results in original record order:
                # the shared molecule file cannot be reordered with workers.
                remaining = itertools.chain((record,), records)
                expand = partial(_expand_nohmm_origin, frame_count=self.step)
                if self.nproc == 1:
                    results = map(expand, remaining)
                    context = ExitStack()
                else:
                    results = run_mp(
                        self.nproc,
                        func=expand,
                        l=remaining,
                        total=self.temp1it - self.hmmit,
                        unordered=False,
                        disk_ordered=True,
                        max_inflight=2 * self.nproc,
                        chunksize=1,
                        bar=False,
                    )
                    context = closing(results)
                with context:
                    for signal in results:
                        appendIfNotNone(origin, signal)
                        self.hmmit += 1
                        progress.update()
                break
        if self.hmmit != self.temp1it:
            raise ValueError("Molecule count differs from declared count")

    def _pooled_results(self, definitions, workers):
        """Order input by physical atom set and dispatch each set as one task."""
        with _OrderedMoleculeResults() as ordered:
            for record in _iter_molecule_records(definitions):
                ordered.add(None, None, record)
            records = (record for _, _, record in ordered)
            tasks = (
                (atoms, self.rng.atomtype[list(atoms)], tuple(group))
                for atoms, group in itertools.groupby(
                    records, key=lambda r: tuple(bytestolist(r[0]))
                )
            )
            count = ordered.atom_set_count()
            filter_set = partial(
                _filter_atom_set,
                miso=self.miso,
                step=self.step,
                model=self._model,
                getoriginfile=self.getoriginfile,
                printfiltersignal=self.printfiltersignal,
            )
            workers = min(workers, count)
            if workers > 1:
                results = run_mp(
                    workers,
                    func=filter_set,
                    l=tasks,
                    total=count,
                    unordered=False,
                    disk_ordered=True,
                    max_inflight=2 * workers,
                    chunksize=1,
                    desc=self.description,
                    unit="atom set",
                )
            else:
                results = tqdm(
                    map(filter_set, tasks),
                    total=count,
                    desc=self.description,
                    unit="atom set",
                    disable=None,
                )
            route_ids = np.zeros(
                self.temp1it + 1, dtype=_unsigned_dtype_for_maximum(self.temp1it)
            )
            n = 0
            with closing(results):
                for batch in results:
                    first = {}
                    for group_id, origin, signal, record in batch:
                        n += 1
                        route_ids[n] = first.setdefault(group_id, n)
                        yield origin, signal, record
            self.hmm_route_ids = route_ids[: n + 1]

    def _calhmm(self):
        """Filter aligned streams with bounded long-signal task admission."""
        if not self.runHMM:
            try:
                self._cal_nohmm()
            except BaseException:
                # Only the newly created signal stream belongs to this stage;
                # the reused definitions must survive failure for their owner.
                if self.originfilename is not None:
                    with suppress(OSError):
                        os.remove(self.originfilename)
                raise
            return
        with (
            WriteBuffer(tempfile.NamedTemporaryFile("wb", delete=False))
            if self.getoriginfile or not self.runHMM
            else ExitStack() as fo,
            WriteBuffer(tempfile.NamedTemporaryFile("wb", delete=False))
            if self.runHMM
            else ExitStack() as fh,
            open(self.moleculetempfilename, "rb") as ft,
            WriteBuffer(tempfile.NamedTemporaryFile("wb", delete=False)) as ft2,
        ):
            self.moleculetemp2filename = ft2.name
            if self.getoriginfile or not self.runHMM:
                assert not isinstance(fo, ExitStack)
                self.originfilename = fo.name
            else:
                self.originfilename = None
            if self.runHMM:
                assert not isinstance(fh, ExitStack)
                self.hmmfilename = fh.name
            else:
                self.hmmfilename = None
            workers = self.nproc
            bounded = self.temp1it and self.step >= _HMM_BOUNDED_SIGNAL_FRAMES
            if bounded:
                # The legacy runner executes 100 records per pool task. Avoid
                # increasing simultaneous Viterbi working arrays when replacing
                # those chunks with individual tasks; small inputs run locally.
                workers = min(
                    workers,
                    (self.temp1it + _HMM_LEGACY_CHUNK_RECORDS - 1)
                    // _HMM_LEGACY_CHUNK_RECORDS,
                )
            results = (
                self._pooled_results(ft, workers)
                if self.miso
                else self._independent_results(ft, workers, bounded)
            )
            self.hmmit = 0
            with closing(results):
                for originbytes, hmmbytes, line_c in results:
                    appendIfNotNone(fo, originbytes)
                    appendIfNotNone(fh, hmmbytes)
                    self.hmmit += 1
                    ft2.extend(line_c)

    def _independent_results(self, ft, workers, bounded):
        """Filter individual structures and emit results in canonical order."""
        limits: dict[str, Any] = {}
        if workers == 1:
            # Reuse the initialized model and consume each record together
            # with its origin and filtered signal, without a worker pool.
            records = itertools.zip_longest(*[read_compressed_block(ft)] * 4)
            results = tqdm(
                map(self._getoriginandhmm, records),
                total=self.temp1it,
                desc=self.description,
                unit="molecule",
                disable=None,
            )
        else:
            # A default 100-record pool chunk retains many signal results.
            # Long signals amortize per-record scheduling; bound submitted
            # but unconsumed records. One decoded signal still scales with
            # the number of frames; this is not a total memory bound.
            # Forkserver per-record transport has higher observed overhead.
            # Retain legacy chunks there, but exclude the RNG context from
            # the callable serialized with each chunk as well.
            chunked_forkserver = bounded and get_start_method() == "forkserver"
            function = (
                partial(
                    _convert_hmm_record,
                    step=self.step,
                    model=self._model,
                    getoriginfile=self.getoriginfile,
                    printfiltersignal=self.printfiltersignal,
                )
                if chunked_forkserver
                else _filter_hmm_record
                if bounded
                else self._getoriginandhmm
            )
            limits = (
                {
                    "chunksize": 1,
                    "max_inflight": 2 * workers,
                    "initializer": _initialize_hmm_worker,
                    "initargs": (
                        self.step,
                        self._model,
                        self.getoriginfile,
                        self.printfiltersignal,
                    ),
                }
                if bounded and not chunked_forkserver
                else {}
            )
            results = run_mp(
                workers,
                func=function,
                l=read_compressed_block(ft),
                nlines=4,
                total=self.temp1it,
                desc=self.description,
                unit="molecule",
                **limits,
            )
        with closing(results), _OrderedMoleculeResults() as ordered:
            for originbytes, hmmbytes, line_c in results:
                if originbytes is not None or hmmbytes is not None:
                    ordered.add(originbytes, hmmbytes, line_c)
            yield from ordered
