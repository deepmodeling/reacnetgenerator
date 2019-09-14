# cython: language_level=3
# cython: linetrace=True
"""HMM Filter.

In order to filter noise, a two-state HMM was adopted, which can be
described as a transition matrix A, an emission matrix B, and an initial
state vector π. The existence of molecules can be converted into 0-1 signals.
In order to predict the state sequence according to the output sequence,
Viterbi Algorithm is used to acquire the path with the most likely hidden
state sequence V called Viterbi path.

Reference:
[1] Wang, L.-P.; McGibbon, R. T.; Pande, V. S.; Martínez, T. J. Automated
discovery and refinement of reactive molecular dynamics pathways. J. Chem.
Theory Comput. 2016, 12(2), 638-649.
[2] Rabiner, L. R. A trtorial on hidden Markov models and selected
applications in speech recognition. Proc. IEEE 1989, 77(2), 257-286.
[3] Forney, G. D. The viterbi algorithm. Porc. IEEE 1973, 61(3), 268-278.
"""

import tempfile
from contextlib import ExitStack
from multiprocessing import Pool, Semaphore

import numpy as np
from hmmlearn import hmm

from .utils import WriteBuffer, appendIfNotNone, bytestolist, listtobytes, multiopen, SharedRNGData


class _HMMFilter(SharedRNGData):
    def __init__(self, rng):
        SharedRNGData.__init__(self, rng, ['runHMM', 'getoriginfile', 'printfiltersignal',
                                           'moleculetempfilename', 'nproc', 'temp1it', 'p', 'a', 'b', 'step'],
                               ['moleculetemp2filename', 'originfilename', 'hmmfilename', 'hmmit'])

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
        self._model = hmm.MultinomialHMM(n_components=2)
        self._model.startprob_ = self.p
        self._model.transmat_ = self.a
        self._model.emissionprob_ = self.b

    def _getoriginandhmm(self, item):
        line_c, _ = item
        value = bytestolist(line_c[-1])
        origin = np.zeros((self.step, 1), dtype=np.int8)
        origin[value] = 1
        originbytes = listtobytes(
            origin) if self.getoriginfile else None
        hmmbytes = None
        if self.runHMM:
            _, hmm = self._model.decode(origin, algorithm="viterbi")
            if 1 in hmm or self.printfiltersignal:
                hmmbytes = listtobytes(hmm)
        return originbytes, hmmbytes, line_c

    def _calhmm(self):
        with WriteBuffer(tempfile.NamedTemporaryFile('wb', delete=False)) if self.getoriginfile or not self.runHMM else ExitStack() as fo, WriteBuffer(tempfile.NamedTemporaryFile('wb', delete=False)) if self.runHMM else ExitStack() as fh, open(self.moleculetempfilename, 'rb') as ft, WriteBuffer(tempfile.NamedTemporaryFile('wb', delete=False)) as ft2:
            pool = Pool(self.nproc, maxtasksperchild=1000)
            try:
                self.moleculetemp2filename = ft2.name
                self.originfilename = fo.name if self.getoriginfile or not self.runHMM else None
                self.hmmfilename = fh.name if self.runHMM else None
                semaphore = Semaphore(self.nproc*150)
                results = multiopen(pool, self._getoriginandhmm, ft, semaphore,
                                    nlines=3, total=self.temp1it, desc="HMM filter", unit="molecule")
                hmmit = 0
                for originbytes, hmmbytes, line_c in results:
                    if originbytes is not None or hmmbytes is not None:
                        appendIfNotNone(fo, originbytes)
                        appendIfNotNone(fh, hmmbytes)
                        hmmit += 1
                        ft2.extend(line_c)
                    semaphore.release()
            finally:
                pool.close()
                pool.join()
        self.hmmit = hmmit
