# SPDX-License-Identifier: LGPL-3.0-or-later
"""Candidate buckets must preserve exact VF2 grouping and name frequencies."""

from types import SimpleNamespace

import networkx as nx
import numpy as np
import pytest

from reacnetgenerator import _path
from reacnetgenerator._path import Molecule, _isomer_candidate_key
from tests.test_isomer_representative import _collect


@pytest.mark.parametrize("seed", range(5))
def test_key_is_invariant_to_atom_ids_and_insertion_order(seed):
    """Relabel shuffled labelled graphs, including loops and isolated atoms."""
    random = np.random.default_rng(seed)
    for size in (0, 1, 7, 20):
        graph = nx.Graph()
        for node in range(size):
            graph.add_node(node, atom=int(random.integers(4)))
        for a in range(size):
            for b in range(a, size):
                if random.random() < 0.18:
                    graph.add_edge(a, b, level=int(random.integers(1, 4)))
        order = random.permutation(size)
        labels = {int(node): int(label) + 101 for node, label in enumerate(order)}
        shuffled = nx.Graph()
        for node in order:
            shuffled.add_node(labels[node], **graph.nodes[node])
        edges = list(graph.edges(data=True))
        for index in random.permutation(len(edges)):
            a, b, data = edges[index]
            shuffled.add_edge(labels[b], labels[a], **data)
        assert _isomer_candidate_key(graph) == _isomer_candidate_key(shuffled)


@pytest.mark.parametrize("miso", [0, 1, 2])
def test_candidate_key_uses_the_selected_isomer_policy(miso):
    """Bond orders are ignored in mode 1 and connectivity in mode 2."""
    context = SimpleNamespace(
        atomtype=np.zeros(4, dtype=int),
        atomnames=np.full(4, "C"),
        atomname=["C"],
        miso=miso,
        convertSMILES=lambda *_: "unused",
    )
    graphs = [
        Molecule(context, np.arange(4), np.asarray(bonds)).graph
        for bonds in (
            [(0, 1, 1), (1, 2, 1), (2, 3, 1)],
            [(0, 1, 2), (1, 2, 1), (2, 3, 1)],
            [(0, 1, 1), (0, 2, 1), (0, 3, 1)],
        )
    ]
    keys = list(map(_isomer_candidate_key, graphs))
    assert (keys[0] == keys[1]) == (miso != 0)
    assert (keys[0] == keys[2]) == (miso == 2)


def test_regular_graph_collision_still_requires_vf2(tmp_path, monkeypatch):
    """A prism and K3,3 share local neighborhoods but are not isomorphic."""
    graphs = [nx.circular_ladder_graph(3), nx.complete_bipartite_graph(3, 3)]
    graphs.append(nx.relabel_nodes(graphs[0], {node: node + 6 for node in graphs[0]}))
    for graph in graphs:
        nx.set_node_attributes(graph, 0, "atom")
        nx.set_edge_attributes(graph, 1, "level")
    assert _isomer_candidate_key(graphs[0]) == _isomer_candidate_key(graphs[1])
    records = [
        (name, list(graph), [(a, b, 1) for a, b in graph.edges], range(count))
        for name, graph, count in zip(
            ("prism", "bipartite", "prism-copy"), graphs, (2, 9, 7)
        )
    ]
    collector = _collect(tmp_path, monkeypatch, records, miso=1)
    assert list(collector.mname) == ["prism-copy", "bipartite", "prism-copy"]


@pytest.mark.parametrize("miso", [0, 1, 2])
def test_forced_hash_collision_preserves_names_and_timeline(
    tmp_path, monkeypatch, miso
):
    """Even a constant key retains structure distinctions and frequency ties."""
    records = [
        ("chain-a", [0, 1, 2, 3], [(0, 1, 1), (1, 2, 1), (2, 3, 1)], [0, 1]),
        ("chain-b", [0, 1, 2, 3], [(0, 1, 2), (1, 2, 1), (2, 3, 1)], [2, 3]),
        ("star", [0, 1, 2, 3], [(0, 1, 1), (0, 2, 1), (0, 3, 1)], range(4, 9)),
        ("star-winner", [4, 5, 6, 7], [(7, 6, 2), (7, 5, 1), (7, 4, 1)], range(9, 16)),
    ]
    normal = tmp_path / "normal"
    normal.mkdir()
    collision = tmp_path / "collision"
    collision.mkdir()
    expected = _collect(normal, monkeypatch, records, miso=miso, print_timeline=True)
    calls = []

    def constant_key(graph):
        """Force collisions after the lazy index has actually been activated."""
        calls.append(graph)
        return 0

    monkeypatch.setattr(_path, "_isomer_candidate_key", constant_key)
    actual = _collect(collision, monkeypatch, records, miso=miso, print_timeline=True)
    assert bool(calls) == (miso == 1)
    assert list(expected.mname) == list(actual.mname)
    for name in ("molecule-names.txt", "timeline.csv"):
        assert (normal / name).read_bytes() == (collision / name).read_bytes()


def test_distinct_neighborhoods_skip_only_impossible_candidates(tmp_path, monkeypatch):
    """Skip the first chain group when a permuted star matches the second group."""
    records = [
        ("chain", [0, 1, 2, 3], [(0, 1, 1), (1, 2, 1), (2, 3, 1)], [0]),
        ("star", [0, 1, 2, 3], [(0, 1, 1), (0, 2, 1), (0, 3, 1)], [1]),
        ("winner", [4, 5, 6, 7], [(7, 6, 2), (7, 5, 1), (7, 4, 1)], [2, 3]),
    ]
    calls = []
    original = Molecule.isomorphic

    def checked(left, right):
        calls.append((left.atoms.tolist(), right.atoms.tolist()))
        return original(left, right)

    monkeypatch.setattr(Molecule, "isomorphic", checked)
    collector = _collect(tmp_path, monkeypatch, records, miso=1)
    assert list(collector.mname) == ["chain", "winner", "winner"]
    # One rejected comparison establishes the second group before promotion.
    assert len(calls) == 2


@pytest.mark.parametrize("miso", [1, 2])
def test_singleton_formulas_do_not_build_keys(tmp_path, monkeypatch, miso):
    """Unique formulas and repeated names retain their original cheap admission."""
    records = [
        ("single", [0], [], [0]),
        ("pair", [0, 1], [(0, 1, 1)], [1]),
        ("pair", [2, 3], [(2, 3, 1)], [2]),
    ]

    def unexpected(graph):
        raise AssertionError("A singleton formula must not compute a key")

    monkeypatch.setattr(_path, "_isomer_candidate_key", unexpected)
    collector = _collect(tmp_path, monkeypatch, records, miso=miso)
    assert list(collector.mname) == ["single", "pair", "pair"]


def test_composition_only_mode_keeps_unindexed_candidates(tmp_path, monkeypatch):
    """Mode 2 has no connectivity to prefilter, including differing canonical names."""
    records = [
        ("chain", [0, 1, 2], [(0, 1, 1), (1, 2, 1)], [0]),
        ("ring", [0, 1, 2], [(0, 1, 1), (1, 2, 1), (2, 0, 1)], [1, 2]),
    ]

    def unexpected(graph):
        raise AssertionError("Composition-only matching must not build graph keys")

    monkeypatch.setattr(_path, "_isomer_candidate_key", unexpected)
    collector = _collect(tmp_path, monkeypatch, records, miso=2)
    assert list(collector.mname) == ["ring", "ring"]


def test_repeated_topology_does_not_build_keys(tmp_path, monkeypatch):
    """Different canonical labels in one structural group do not need an index."""
    records = [
        (f"variant-{order}", [0, 1], [(0, 1, order)], range(order))
        for order in (1, 2, 3)
    ]

    def unexpected(graph):
        raise AssertionError("A single structural group must stay unindexed")

    monkeypatch.setattr(_path, "_isomer_candidate_key", unexpected)
    collector = _collect(tmp_path, monkeypatch, records, miso=1)
    assert list(collector.mname) == ["variant-3"] * 3
