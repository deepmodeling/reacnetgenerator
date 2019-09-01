# cython: language_level=3
'''A beta version for reactions finder.


'''
import numpy as np
from collections import Counter
from multiprocessing import Pool, Semaphore
from tqdm import tqdm


class ReactionsFinder:
    CONFLICT = -1
    EMPTY = 0

    def __init__(self, rng):
        self.rng = rng
        self._step = rng.step
        self._mname = rng.mname
        self._reactionabcdfilename = rng.reactionabcdfilename
        self.nproc = rng.nproc
        self.produce = rng.produce

    def findreactions(self, atomeach, conflict):
        allreactions = []
        with Pool(self.nproc, maxtasksperchild=1000) as pool:
            semaphore = Semaphore(self.nproc*150)
            # atomeach j, atomeach j+1, conflict j, conflict j+1
            givenarray = [(atomeach[:, j], atomeach[:, j+1], conflict[:, j], conflict[:, j+1]) for j in range(self._step-1)]
            results = pool.imap_unordered(self._getstepreaction, self.produce(
                semaphore, givenarray, ()), 100)
            for networks in tqdm(
                    results, total=self._step-1, desc="Analyze reactions (A+B->C+D)",
                    unit="timestep"):
                allreactions.extend(networks)
                semaphore.release()
            pool.close()
            pool.join()
        # reaction with SMILES
        allreactionswithname = Counter(allreactions).most_common()
        buff = '\n'.join([f"{number} {reaction}" for reaction,
                          number in allreactionswithname if reaction is not None])
        with open(self._reactionabcdfilename, 'w') as f:
            f.write(buff)

    def _getstepreaction(self, item):
        (atomeachj, atomeachjp1, conflictj, conflictjp1), _ = item
        networks = []
        modifiedatoms = np.where(np.not_equal(atomeachj, atomeachjp1))[0]
        for i in modifiedatoms:
            networks.append([({atomeachj[i], self.CONFLICT} if conflictj[i] else {atomeachj[i]}), ({
                            atomeachjp1[i], self.CONFLICT} if conflictjp1[i] else {atomeachjp1[i]})])
        a1 = 0
        while a1 < len(networks):
            n1 = networks[a1]
            for a2, n2 in enumerate(networks[a1+1:], a1+1):
                if n1[0] & n2[0] or n1[1] & n2[1]:
                    networks.pop(a2)
                    networks.pop(a1)
                    networks.append([n1[0] | n2[0], n1[1] | n2[1]])
                    break
            else:
                a1 += 1

        # remove empty AND conflict
        a = 0
        while a < len(networks):
            n = networks[a]
            if self.EMPTY in n[0] or self.EMPTY in n[1] or self.CONFLICT in n[0] or self.CONFLICT in n[1]:
                networks.pop(a)
            else:
                a += 1
        # reaction with SMILES name like A+B->C+D
        networks = [self._filterspec(reaction) for reaction in networks]
        return networks

    def _filterspec(self, reaction):
        leftname, rightname = [[self._mname[spec-1]
                                for spec in side] for side in reaction]
        # remove duplicate species
        i = 0
        while i < len(leftname):
            for j, _ in enumerate(rightname):
                if leftname[i] == rightname[j]:
                    leftname.pop(i)
                    rightname.pop(j)
                    break
            else:
                i += 1
        if leftname and rightname:
            return '->'.join(('+'.join(sorted(side)) for side in (leftname, rightname)))
        return None
