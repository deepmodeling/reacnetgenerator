# SPDX-License-Identifier: LGPL-3.0-or-later
# cython: language_level=3
# cython: linetrace=True
"""Reactions finder."""

import csv
from collections import Counter, defaultdict
from multiprocessing.util import Finalize
from typing import Any

import numpy as np

from ._step3state import _STEP3_SCAN_ROWS, _AtomFrameReader
from .dps import dps_reaction  # type:ignore
from .utils import (
    SharedRNGData,
    WriteBuffer,
    bytestolist,
    listtobytes,
    run_mp,
)

_REACTION_WORKER_STATE = None


def _initialize_reaction_worker(reader_args, names, printreactionevent):
    """Attach shared matrices and compact metadata once per worker."""
    global _REACTION_WORKER_STATE
    reader = _AtomFrameReader(*reader_args)
    Finalize(None, reader.close, exitpriority=10)
    # Workers need only species lookup and reaction formatting, not the RNG
    # object or the output paths retained by the parent-side finder.
    finder = object.__new__(ReactionsFinder)
    finder.mname = names
    finder.printreactionevent = printreactionevent
    _REACTION_WORKER_STATE = reader, finder


def _get_step_reaction_by_index(transition_index):
    """Read one transition in bounded atom blocks from the shared files."""
    assert _REACTION_WORKER_STATE is not None
    reader, finder = _REACTION_WORKER_STATE
    blocks = (
        (
            reader.atomeach[start : start + _STEP3_SCAN_ROWS, transition_index],
            reader.atomeach[start : start + _STEP3_SCAN_ROWS, transition_index + 1],
            reader.conflict[start : start + _STEP3_SCAN_ROWS, transition_index],
            reader.conflict[start : start + _STEP3_SCAN_ROWS, transition_index + 1],
        )
        for start in range(0, reader.atomeach.shape[0], _STEP3_SCAN_ROWS)
    )
    return finder._reactions_from_blocks(
        blocks, transition_index if finder.printreactionevent else None
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

    def findreactions(self, atomeach, conflict, *, matrix_store=None):
        """Analyze indexed shared state, or accept the legacy array inputs."""
        if matrix_store is not None:
            results = run_mp(
                self.nproc,
                func=_get_step_reaction_by_index,
                l=matrix_store.iter_active_transitions(),
                initializer=_initialize_reaction_worker,
                initargs=(
                    matrix_store.reader_args,
                    self.mname,
                    self.printreactionevent,
                ),
                # Preserve the legacy reduction policy: event rows are ordered;
                # count-only output consumes completed workers without ordering.
                unordered=not self.printreactionevent,
                chunksize=1,
                max_inflight=max(2, 2 * self.nproc),
                disk_ordered=self.printreactionevent,
                total=matrix_store.active_transition_count,
                desc="Analyze reactions (A+B->C+D)",
                unit="timestep",
            )
            try:
                self._write_reactions(results)
            finally:
                # Close workers before collect() removes their shared files.
                results.close()
            return
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
        ordered_kwargs: dict[str, Any] = (
            {
                "chunksize": 1,
                "max_inflight": max(2, 2 * self.nproc),
                "disk_ordered": True,
            }
            if self.printreactionevent
            else {}
        )
        results = run_mp(
            self.nproc,
            func=self._getstepreaction,
            l=givenarray,
            unordered=not self.printreactionevent,
            total=self.step - 1,
            desc="Analyze reactions (A+B->C+D)",
            unit="timestep",
            **ordered_kwargs,
        )
        self._write_reactions(results)

    def _write_reactions(self, results):
        """Consume results incrementally in the existing text/CSV formats."""
        reaction_counts = Counter()
        if self.printreactionevent:
            with open(self.reactioneventfilename, "w", newline="") as f_event:
                event_writer = csv.writer(f_event)
                event_writer.writerow(["Timestep_Index", "Reactant", "Product"])
                for events in results:
                    for event in events:
                        reaction = "->".join((event["Reactant"], event["Product"]))
                        reaction_counts[reaction] += 1
                        event_writer.writerow(
                            [
                                event["Timestep_Index"],
                                event["Reactant"],
                                event["Product"],
                            ]
                        )
        else:
            for reactions in results:
                reaction_counts.update(
                    reaction for reaction in reactions if reaction is not None
                )
        # reaction with SMILES
        allreactionswithname = reaction_counts.most_common()
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
        blocks = (
            tuple(
                values[start : start + _STEP3_SCAN_ROWS]
                for values in (atomeachj, atomeachjp1, conflictj, conflictjp1)
            )
            for start in range(0, len(atomeachj), _STEP3_SCAN_ROWS)
        )
        return self._reactions_from_blocks(blocks, stepidx)

    def _reactions_from_blocks(self, blocks, stepidx):
        """Build the same reaction graph without full-frame temporary arrays.

        Blocks visit atoms in their original order. The graph for one transition
        can still grow with its participating atoms; only the scan temporaries
        have a fixed bound.
        """
        reactdict = [defaultdict(list), defaultdict(list)]
        for before, after, before_conflict, after_conflict in blocks:
            for atom in np.flatnonzero(before != after):
                left, right = int(before[atom]), int(after[atom])
                reactdict[0][left].append(right)
                reactdict[1][right].append(left)
                if before_conflict[atom]:
                    reactdict[0][left].append(self.CONFLICT)
                if after_conflict[atom]:
                    reactdict[1][right].append(self.CONFLICT)
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
