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
from tqdm import tqdm


class _HMMFilter:
    def __init__(self, rng):
        self.rng = rng
        self.runHMM = rng.runHMM
        self.originfilename = rng.originfilename
        self.hmmfilename = rng.hmmfilename
        self.getoriginfile = rng.getoriginfile
        self.printfiltersignal = rng.printfiltersignal
        self.moleculetempfilename = rng.moleculetempfilename
        self.nproc = rng.nproc
        self._temp1it = rng.temp1it
        self.p = rng.p
        self.a = rng.a
        self.b = rng.b
        self._step = rng.step
        self._compress = rng.compress
        self._decompress = rng.decompress
        self._produce = rng.produce
        self._model = None
        self.moleculetemp2filename = None
        self._hmmit = None

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
        self.rng.moleculetemp2filename = self.moleculetemp2filename

    def _initHMM(self):
        self._model = hmm.MultinomialHMM(n_components=2)
        self._model.startprob_ = self.p
        self._model.transmat_ = self.a
        self._model.emissionprob_ = self.b

    def _getoriginandhmm(self, item):
        line_c, _ = item
        line = self._decompress(line_c)
        s = line.split()
        value = np.array([int(x) for x in s[-1].split(",")])
        origin = np.zeros(self._step, dtype=np.int8)
        origin[value-1] = 1
        if self.runHMM:
            _, hmm = self._model.decode(
                np.array([origin]).T, algorithm="viterbi")
        return origin, (np.array(hmm) if self.runHMM else np.zeros(0)), line

    def _calhmm(self):
        with open(self.originfilename, 'wb') if self.getoriginfile or not self.runHMM else ExitStack() as fo, open(self.hmmfilename, 'wb') if self.runHMM else ExitStack() as fh, open(self.moleculetempfilename, 'rb') as ft, tempfile.NamedTemporaryFile('wb', delete=False) as ft2, Pool(self.nproc, maxtasksperchild=1000) as pool:
            self.moleculetemp2filename = ft2.name
            semaphore = Semaphore(self.nproc*150)
            results = pool.imap_unordered(
                self._getoriginandhmm, self._produce(semaphore, ft, ()), 100)
            hmmit = 0
            for originsignal, hmmsignal, mlist in tqdm(
                    results, total=self._temp1it, desc="HMM filter",
                    unit="molecule"):
                if 1 in hmmsignal or self.printfiltersignal or not self.runHMM:
                    if self.getoriginfile:
                        fo.write(self._compress(
                            "".join([str(i) for i in originsignal.tolist()])))
                    if self.runHMM:
                        fh.write(self._compress(
                            "".join([str(i) for i in hmmsignal.tolist()])))
                    hmmit += 1
                    ft2.write(self._compress(mlist.strip()))
                semaphore.release()
        self._hmmit = hmmit
        pool.join()
