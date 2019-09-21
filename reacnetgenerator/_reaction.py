# cython: language_level=3
# cython: linetrace=True
'''A beta version for reactions finder.


'''
import numpy as np
from collections import Counter

from .utils import WriteBuffer, run_mp, SharedRNGData


class ReactionsFinder(SharedRNGData):
    CONFLICT = -1
    EMPTY = 0

    def __init__(self, rng):
        SharedRNGData.__init__(self, rng, ["step", "mname",
                                           "reactionabcdfilename", "nproc"], [])

    def findreactions(self, atomeach, conflict):
        allreactions = []
        # atomeach j, atomeach j+1, conflict j, conflict j+1
        givenarray = [(atomeach[:, j], atomeach[:, j+1], conflict[:, j],
                        conflict[:, j+1]) for j in range(self.step-1)]
        results = run_mp(self.nproc, func=self._getstepreaction, l=givenarray,
                            total=self.step-1, desc="Analyze reactions (A+B->C+D)", unit="timestep")
        for networks in results:
            allreactions.extend(networks)
        # reaction with SMILES
        allreactionswithname = Counter(allreactions).most_common()
        with WriteBuffer(open(self.reactionabcdfilename, 'w'), sep='\n') as f:
            for reaction, number in allreactionswithname:
                if reaction is not None:
                    f.append(f"{number} {reaction}")

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
        leftname, rightname = [[self.mname[spec-1]
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
