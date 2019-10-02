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

from hmmlearn import hmm

from .utils import WriteBuffer, appendIfNotNone, bytestolist, listtobytes, run_mp, SharedRNGData
from .utils_np import idx_to_signal, check_zero_signal


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
        self._model = hmm.MultinomialHMM(n_components=2, algorithm="viterbi")
        self._model.startprob_ = self.p
        self._model.transmat_ = self.a
        self._model.emissionprob_ = self.b

    def _getoriginandhmm(self, item):
        line_c = item
        value = bytestolist(line_c[-1])
        origin = idx_to_signal(value, self.step)
        originbytes = listtobytes(
            origin) if self.getoriginfile else None
        hmmbytes = None
        if self.runHMM:
            hmmsignal = self._model.predict(origin)
            if check_zero_signal(hmmsignal) or self.printfiltersignal:
                hmmbytes = listtobytes(hmmsignal)
        return originbytes, hmmbytes, line_c

    def _calhmm(self):
        with WriteBuffer(tempfile.NamedTemporaryFile('wb', delete=False)) if self.getoriginfile or not self.runHMM else ExitStack() as fo, WriteBuffer(tempfile.NamedTemporaryFile('wb', delete=False)) if self.runHMM else ExitStack() as fh, open(self.moleculetempfilename, 'rb') as ft, WriteBuffer(tempfile.NamedTemporaryFile('wb', delete=False)) as ft2:
            self.moleculetemp2filename = ft2.name
            self.originfilename = fo.name if self.getoriginfile or not self.runHMM else None
            self.hmmfilename = fh.name if self.runHMM else None
            results = run_mp(self.nproc, func=self._getoriginandhmm, l=ft,
                                nlines=3, total=self.temp1it, desc="HMM filter", unit="molecule")
            hmmit = 0
            for originbytes, hmmbytes, line_c in results:
                if originbytes is not None or hmmbytes is not None:
                    appendIfNotNone(fo, originbytes)
                    appendIfNotNone(fh, hmmbytes)
                    hmmit += 1
                    ft2.extend(line_c)
        self.hmmit = hmmit
