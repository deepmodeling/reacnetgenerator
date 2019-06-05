'''A beta version for reactions finder.


'''
import numpy as np
from collections import Counter


class ReactionsFinder:
    CONFLICT = -1
    EMPTY = 0

    def __init__(self, rng):
        self.rng = rng
        self._step = rng.step
        self._mname = rng.mname
        self._reactionabcdfilename = rng.reactionabcdfilename

    def findreactions(self, atomeach, conflict):
        allreactions = Counter()

        for j in range(self._step-1):
            networks = []
            modifiedatoms = np.where(atomeach[:, j] != atomeach[:, j+1])[0]
            for i in modifiedatoms:
                networks.append([({atomeach[i, j], self.CONFLICT} if conflict[i, j] else {atomeach[i, j]}), ({
                                atomeach[i, j], self.CONFLICT} if conflict[i, j+1] else {atomeach[i, j+1]})])
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

            networks = [(tuple(sorted(r)) for r in n) for n in networks]
            allreactions |= Counter(networks)
        # reaction with SMILES
        allreactionswithname = Counter(
            [self.filterspec(reaction) for reaction in allreactions])
        buff = '\n'.join([f"{'->'.join(('+'.join((spec for spec in side)) for side in reaction))} {number}" for reaction,
                          number in allreactionswithname.most_common() if reaction is not None])
        with open(self._reactionabcdfilename, 'w') as f:
            f.write(buff)

    def filterspec(self, reaction):
        leftname, rightname = [[self._mname[spec-1]
                                for spec in side] for side in reaction]
        i = 0
        while i < len(leftname):
            for j, _ in enumerate(rightname):
                if leftname[i] == rightname[j]:
                    leftname.pop(i)
                    rightname.pop(j)
                    break
            else:
                i += 1
        if len(leftname) and len(rightname):
            return tuple(sorted(leftname)), tuple(sorted(rightname))
        return None
