# SPDX-License-Identifier: LGPL-3.0-or-later
"""Deduplicate graph edges while preserving conflicts and reaction instances."""

from collections import defaultdict

import numpy as np
import pytest
from reacnetgenerator.dps import dps_reaction

from reacnetgenerator._reaction import (
    ReactionsFinder,
    _add_reaction_neighbor,
    _reaction_pair_first_indices,
)


def _finder(events=True):
    finder = object.__new__(ReactionsFinder)
    finder.mname = np.array(["A", "B", "C", "A", "B", "C"], dtype=object)
    finder.printreactionevent = events
    finder._include_transition_evidence = events
    return finder


def _legacy_events(finder, blocks):
    """Build the original per-atom multigraph as an independent oracle."""
    graph = [defaultdict(list), defaultdict(list)]
    for before, after, left_flags, right_flags in blocks:
        for atom in np.flatnonzero(before != after):
            left, right = int(before[atom]), int(after[atom])
            graph[0][left].append(right)
            graph[1][right].append(left)
            if left_flags[atom]:
                graph[0][left].append(-1)
            if right_flags[atom]:
                graph[1][right].append(-1)
    result = []
    for left, right in dps_reaction(graph):
        if {-1, 0}.intersection(left) or {-1, 0}.intersection(right):
            continue
        names = finder._filterreactionpair((left, right))
        if not finder.printreactionevent:
            result.append("->".join(names) if names is not None else None)
        elif names is not None:
            result.append(
                {
                    "Timestep_Index": 7,
                    "Reactant": names[0],
                    "Product": names[1],
                    "ReactantMoleculeIDs": tuple(sorted(left)),
                    "ProductMoleculeIDs": tuple(sorted(right)),
                }
            )
    return result


@pytest.mark.parametrize("events", [False, True])
@pytest.mark.parametrize("conflict_side", [None, 0, 1])
def test_graph_pairs_are_unique_across_blocks_without_losing_late_conflicts(
    monkeypatch, events, conflict_side
):
    """A non-representative duplicate may carry the only conflict evidence."""
    finder = _finder(events)
    monkeypatch.setattr(
        "reacnetgenerator._reaction._REACTION_PAIR_COMPACTION_MIN_ROWS", 4
    )
    blocks = []
    for _ in range(3):
        before = np.tile([1, 2, 4, 5], 300)
        after = np.tile([3, 3, 6, 6], 300)
        blocks.append(
            (
                before,
                after,
                np.zeros(len(before), dtype=int),
                np.zeros(len(before), dtype=int),
            )
        )
    if conflict_side is not None:
        blocks[-1][2 + conflict_side][-4] = 1
    captured = []

    def inspect(graph):
        captured.append(
            [{key: list(value) for key, value in side.items()} for side in graph]
        )
        for side in graph:
            assert all(len(value) == len(set(value)) for value in side.values())
        return dps_reaction(graph)

    monkeypatch.setattr("reacnetgenerator._reaction.dps_reaction", inspect)
    actual = finder._reactions_from_blocks(blocks, 7 if events else None)
    assert actual == _legacy_events(finder, blocks)
    assert len(actual) == (2 if conflict_side is None else 1)
    assert sum(len(value) for side in captured[0] for value in side.values()) == 8 + (
        conflict_side is not None
    )
    if events and conflict_side is None:
        assert [event["ReactantMoleculeIDs"] for event in actual] == [(1, 2), (4, 5)]


@pytest.mark.parametrize("events", [False, True])
def test_random_graphs_match_original_connected_components(monkeypatch, events):
    """Empty nodes, conflicts, cycles and blocks retain the baseline result order."""
    monkeypatch.setattr(
        "reacnetgenerator._reaction._REACTION_PAIR_COMPACTION_MIN_ROWS", 4
    )
    random = np.random.default_rng(765)
    finder = _finder(events)
    for iteration in range(40):
        values = np.stack((random.choice([1, 2, 4, 5], 45), random.choice([3, 6], 45)))
        if iteration % 3 == 0:
            values[0, ::11] = 0
        if iteration % 4 == 0:
            values[:, ::7] = values[::-1, ::7]
        flags = random.random((2, 45)) < 0.005
        blocks = [
            (
                values[0, start : start + 9],
                values[1, start : start + 9],
                flags[0, start : start + 9],
                flags[1, start : start + 9],
            )
            for start in range(0, 45, 9)
        ]
        assert finder._reactions_from_blocks(
            blocks, 7 if events else None
        ) == _legacy_events(finder, blocks)


@pytest.mark.parametrize("wide", [False, True])
def test_pair_first_indices_preserve_order_and_integer_precision(wide):
    """Packed small IDs and separate wide fields must choose the same first rows."""
    before = np.array([3, 1, 3, 2, 1, 3], dtype=np.int64)
    after = np.array([4, 5, 4, 4, 5, 6], dtype=np.uint64)
    if wide:
        before += 2**53
        after += np.uint64(2**63)
    np.testing.assert_array_equal(
        _reaction_pair_first_indices(before, after), [0, 1, 3, 5]
    )


def test_high_degree_neighbors_promote_without_changing_order():
    """Avoid duplicate storage and quadratic scans for large neighborhoods."""
    graph = defaultdict(list)
    neighbors = [10, 3, 1, 5, 6, 8, 9, 11, 12, -1, 4, 10, -1]
    for neighbor in neighbors:
        _add_reaction_neighbor(graph, 1, neighbor)
    assert isinstance(graph[1], dict)
    assert list(graph[1]) == list(dict.fromkeys(neighbors))
