# SPDX-License-Identifier: LGPL-3.0-or-later
# cython: language_level=3
# cython: linetrace=True
"""Reactions finder."""

import csv
from collections import Counter, defaultdict

import numpy as np

from .dps import dps_reaction  # type:ignore
from .utils import (
    SharedRNGData,
    WriteBuffer,
    bytestolist,
    listtobytes,
    run_mp,
)


class ReactionsFinder(SharedRNGData):
    CONFLICT = -1
    EMPTY = 0

    step: int
    mname: np.ndarray
    reactionabcdfilename: str
    reactioneventfilename: str
    printreactionevent: bool
    nproc: int

    def __init__(self, rng):
        SharedRNGData.__init__(
            self,
            rng,
            [
                "step",
                "mname",
                "reactionabcdfilename",
                "reactioneventfilename",
                "printreactionevent",
                "nproc",
            ],
            [],
        )

    def findreactions(self, atomeach, conflict):
        allreactions = []
        # atomeach j, atomeach j+1, conflict j, conflict j+1
        if self.printreactionevent:
            givenarray = (
                listtobytes((i, *x))
                for i, x in enumerate(
                    zip(atomeach[:-1], atomeach[1:], conflict[:-1], conflict[1:])
                )
            )
        else:
            givenarray = (
                listtobytes(x)
                for x in zip(atomeach[:-1], atomeach[1:], conflict[:-1], conflict[1:])
            )
        results = run_mp(
            self.nproc,
            func=self._getstepreaction,
            l=givenarray,
            unordered=not self.printreactionevent,
            total=self.step - 1,
            desc="Analyze reactions (A+B->C+D)",
            unit="timestep",
        )
        if self.printreactionevent:
            with open(self.reactioneventfilename, "w", newline="") as f_event:
                event_writer = csv.writer(f_event)
                event_writer.writerow(["Timestep_Index", "Reactant", "Product"])
                for events in results:
                    for event in events:
                        reaction = "->".join((event["Reactant"], event["Product"]))
                        allreactions.append(reaction)
                        event_writer.writerow(
                            [
                                event["Timestep_Index"],
                                event["Reactant"],
                                event["Product"],
                            ]
                        )
        else:
            for reactions in results:
                allreactions.extend(reactions)
        # reaction with SMILES
        allreactionswithname = Counter(allreactions).most_common()
        with WriteBuffer(open(self.reactionabcdfilename, "w"), sep="\n") as f:
            for reaction, number in allreactionswithname:
                if reaction is not None:
                    f.append(f"{number} {reaction}")

    def _getstepreaction(self, item):
        # atomeachj, atomeachjp1, conflictj, conflictjp1
        # or stepidx, atomeachj, atomeachjp1, conflictj, conflictjp1
        item = bytestolist(item)
        if self.printreactionevent:
            stepidx = item[0]
            atomeachj, atomeachjp1, conflictj, conflictjp1 = item[1:]
        else:
            stepidx = None
            atomeachj, atomeachjp1, conflictj, conflictjp1 = item
        modifiedatoms = np.not_equal(atomeachj, atomeachjp1)
        # covert to dict
        reactdict = [defaultdict(list), defaultdict(list)]
        for mol in np.array((atomeachj, atomeachjp1, conflictj, conflictjp1))[
            :, modifiedatoms
        ].T:
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
            if not (
                self.EMPTY in nn[0]
                or self.EMPTY in nn[1]
                or self.CONFLICT in nn[0]
                or self.CONFLICT in nn[1]
            ):
                new_networks.append(nn)
        if not self.printreactionevent:
            # reaction with SMILES name like A+B->C+D
            return [self._filterspec(reaction) for reaction in new_networks]
        events = []
        assert stepidx is not None
        for reaction in new_networks:
            reactionpair = self._filterreactionpair(reaction)
            if reactionpair is None:
                continue
            reactant, product = reactionpair
            events.append(
                {
                    "Timestep_Index": int(stepidx),
                    "Reactant": reactant,
                    "Product": product,
                }
            )
        return events

    def _filterreactionpair(self, reaction):
        leftname, rightname = (
            Counter(self.mname[np.array(side) - 1]) for side in reaction
        )
        # remove duplicate species
        new_leftname = leftname - rightname
        new_rightname = rightname - leftname
        if new_leftname and new_rightname:
            return tuple(
                "+".join(sorted(side.elements()))
                for side in (new_leftname, new_rightname)
            )
        return None

    def _filterspec(self, reaction):
        reactionpair = self._filterreactionpair(reaction)
        if reactionpair is None:
            return None
        return "->".join(reactionpair)
