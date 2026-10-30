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

import numpy as np
from tqdm.auto import tqdm

try:
    # hmmlearn v0.2.8 renamed MultinomialHMM to CategoricalHMM
    from hmmlearn.hmm import CategoricalHMM as MultinomialHMM
except ImportError:
    from hmmlearn.hmm import MultinomialHMM

from ._moleculepresence import _iter_molecule_records, _use_frame_indices
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


def _expand_nohmm_origin(record, *, frame_count):
    """Expand a dense record without sending trajectory metadata to a worker."""
    frames = record[0]
    if _use_frame_indices(frames, frame_count):
        return None
    return listtobytes(idx_to_signal(bytestolist(frames), frame_count))


class _HMMFilter(SharedRNGData):
    HMM_DESCRIPTION = "HMM filter"
    NOHMM_DESCRIPTION = "Pass through original signals (HMM disabled)"

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
            ["moleculetemp2filename", "originfilename", "hmmfilename", "hmmit"],
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
        line_c = item
        value = bytestolist(line_c[-1])
        origin = idx_to_signal(value, self.step)
        originbytes = listtobytes(origin) if self.getoriginfile else None
        hmmbytes = None
        if self.runHMM:
            hmmsignal = self._model.predict(origin).astype(bool)
            # Origin output retains metadata even for filtered-out molecules.
            if (
                check_zero_signal(hmmsignal)
                or self.printfiltersignal
                or self.getoriginfile
            ):
                hmmbytes = listtobytes(hmmsignal)
        return originbytes, hmmbytes, line_c

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

    def _calhmm(self):
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
            if self.nproc == 1:
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
                results = run_mp(
                    self.nproc,
                    func=self._getoriginandhmm,
                    l=read_compressed_block(ft),
                    nlines=4,
                    total=self.temp1it,
                    desc=self.description,
                    unit="molecule",
                )
            hmmit = 0
            for originbytes, hmmbytes, line_c in results:
                if originbytes is not None or hmmbytes is not None:
                    appendIfNotNone(fo, originbytes)
                    appendIfNotNone(fh, hmmbytes)
                    hmmit += 1
                    ft2.extend(line_c)
        self.hmmit = hmmit
