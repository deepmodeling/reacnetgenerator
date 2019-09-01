# cython: language_level=3
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
import itertools
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
        self._bytestolist = rng.bytestolist
        self._listtobytes = rng.listtobytes
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
        self.rng.originfilename = self.originfilename
        self.rng.hmmfilename = self.hmmfilename
        self.rng.hmmit = self._hmmit

    def _initHMM(self):
        self._model = hmm.MultinomialHMM(n_components=2)
        self._model.startprob_ = self.p
        self._model.transmat_ = self.a
        self._model.emissionprob_ = self.b

    def _getoriginandhmm(self, item):
        line_c, _ = item
        value = self._bytestolist(line_c[-1])
        origin = np.zeros((self._step, 1), dtype=np.int8)
        origin[value] = 1
        originbytes = self._listtobytes(
            origin) if self.getoriginfile else None
        hmmbytes = None
        if self.runHMM:
            _, hmm = self._model.decode(origin, algorithm="viterbi")
            if 1 in hmm or self.printfiltersignal:
                hmmbytes = self._listtobytes(hmm)
        return originbytes, hmmbytes, line_c

    def _calhmm(self):
        with tempfile.NamedTemporaryFile('wb', delete=False) if self.getoriginfile or not self.runHMM else ExitStack() as fo, tempfile.NamedTemporaryFile('wb', delete=False) if self.runHMM else ExitStack() as fh, open(self.moleculetempfilename, 'rb') as ft, tempfile.NamedTemporaryFile('wb', delete=False) as ft2, Pool(self.nproc, maxtasksperchild=1000) as pool:
            self.moleculetemp2filename = ft2.name
            self.originfilename = fo.name if self.getoriginfile or not self.runHMM else None
            self.hmmfilename = fh.name if self.runHMM else None
            semaphore = Semaphore(self.nproc*150)
            results = pool.imap_unordered(
                self._getoriginandhmm, self._produce(semaphore, itertools.zip_longest(*[ft] * 3), ()), 100)
            hmmit = 0
            buffo = []
            buffh = []
            bufft = []
            for originbytes, hmmbytes, line_c in tqdm(
                    results, total=self._temp1it, desc="HMM filter",
                    unit="molecule"):
                if originbytes is not None or hmmbytes is not None:
                    for signalbytes, buff, f in [(originbytes, buffo, fo), (hmmbytes, buffh, fh)]:
                        if signalbytes is not None:
                            buff.append(signalbytes)
                            if len(buff) > 30*self.nproc:
                                f.write(b''.join(buff))
                                buff[:] = []
                    hmmit += 1
                    bufft.extend(line_c)
                    if len(buffh) > 30*self.nproc:
                        ft2.write(b''.join(bufft))
                        bufft[:] = []
                semaphore.release()
            for buff, f in [(buffo, fo), (buffh, fh), (bufft, ft2)]:
                if buff:
                    f.write(b''.join(buff))
        pool.close()
        self._hmmit = hmmit
        pool.join()
