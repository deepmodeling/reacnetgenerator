# cython: language_level=3
# cython: linetrace=True
'''A beta version for reactions finder.


'''
import numpy as np
from collections import Counter, defaultdict

from .utils import WriteBuffer, run_mp, SharedRNGData, listtobytes, bytestolist
from .dps import dps_reaction


class ReactionsFinder(SharedRNGData):
    CONFLICT = -1
    EMPTY = 0

    def __init__(self, rng):
        SharedRNGData.__init__(self, rng, ["step", "mname",
                                           "reactionabcdfilename", "nproc"], [])

    def findreactions(self, atomeach, conflict):
        allreactions = []
        # atomeach j, atomeach j+1, conflict j, conflict j+1
        givenarray = zip(atomeach[:-1], atomeach[1:],
                         conflict[:-1], conflict[1:])
        results = run_mp(self.nproc, func=self._getstepreaction, l=[listtobytes(x) for x in givenarray],
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
        # atomeachj, atomeachjp1, conflictj, conflictjp1
        item = bytestolist(item)
        modifiedatoms = np.not_equal(item[0], item[1])
        # covert to dict
        reactdict = [defaultdict(list), defaultdict(list)]
        for mol in np.array(item)[:, modifiedatoms].T:
            reactdict[0][mol[0]].append(mol[1])
            reactdict[1][mol[1]].append(mol[0])
            if mol[2]:
                reactdict[0][mol[0]].append(self.CONFLICT)
            if mol[3]:
                reactdict[1][mol[1]].append(self.CONFLICT)
        networks = dps_reaction(reactdict)
        # remove empty AND conflict
        new_networks = []
        for nn in networks:
            if not (self.EMPTY in nn[0] or self.EMPTY in nn[1] or self.CONFLICT in nn[0] or self.CONFLICT in nn[1]):
                new_networks.append(nn)
        # reaction with SMILES name like A+B->C+D
        reactions = [self._filterspec(reaction) for reaction in new_networks]
        return reactions

    def _filterspec(self, reaction):
        leftname, rightname = (
            Counter(self.mname[np.array(side)-1]) for side in reaction)
        # remove duplicate species
        new_leftname = leftname - rightname
        new_rightname = rightname - leftname
        if new_leftname and new_rightname:
            return '->'.join(('+'.join(sorted(side.elements())) for side in (new_leftname, new_rightname)))
        return None
