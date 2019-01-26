''' HMM Filter '''

from multiprocessing import Pool, Semaphore
from contextlib import ExitStack
import tempfile

import numpy as np
from tqdm import tqdm
from hmmlearn import hmm


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

    def filter(self):
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
            for originsignal, hmmsignal, mlist in tqdm(results, total=self._temp1it, desc="HMM filter", unit="molecule"):
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
