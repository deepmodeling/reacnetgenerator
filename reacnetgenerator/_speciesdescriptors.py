# SPDX-License-Identifier: LGPL-3.0-or-later
"""Searchable properties of the explicit atoms and undirected molecule edges."""

from collections import Counter
from itertools import pairwise

import networkx as nx


def hill_formula(elements):
    """Count recorded atom labels in Hill order, without adding hydrogens."""
    counts = Counter(elements)
    order = sorted(counts)
    if "C" in counts:
        first = ["C", "H"] if "H" in counts else ["C"]
        order = first + sorted(counts.keys() - set(first))
    return "".join(
        symbol if counts[symbol] == 1 else f"{symbol}{counts[symbol]}"
        for symbol in order
    )


def _fused_chain(core):
    """Read consecutive face sizes of a biconnected outerplanar ring block."""
    # Adding a universal vertex tests outerplanarity. In this embedding the
    # faces avoiding the extra vertex are the unique interior rings. A bridged
    # or cage topology must not acquire a fused label from an arbitrary basis.
    apex = -1  # Molecule atom indices are nonnegative.
    graph = core.copy()
    graph.add_edges_from((apex, atom) for atom in core)
    planar, embedding = nx.check_planarity(graph)
    if not planar:
        return "polycyclic"

    visited = set()
    owner = {}
    sizes = []
    dual = nx.Graph()
    for u, v in embedding.edges():
        if (u, v) in visited:
            continue
        face = embedding.traverse_face(u, v, visited)
        if apex in face:
            continue
        ring = len(sizes)
        sizes.append(len(face))
        dual.add_node(ring)
        for a, b in pairwise([*face, face[0]]):
            edge = (min(a, b), max(a, b))
            if edge in owner:
                dual.add_edge(owner[edge], ring)
            else:
                owner[edge] = ring
    ends = [ring for ring, degree in dual.degree if degree == 1]
    if len(ends) != 2:
        return "polycyclic"
    path = [sizes[i] for i in nx.shortest_path(dual, *ends)]
    return "fused-" + "-".join(map(str, max(path, path[::-1])))


def ring_system(atom_count, bonds):
    """Classify rings in a connected RNG molecule using bond endpoints only."""
    # Every detected molecule is one connected component; trees include
    # isolated atoms and need no graph allocation.
    if len(bonds) < atom_count:
        return "acyclic"
    graph = nx.from_edgelist(bonds)
    blocks = [nodes for nodes in nx.biconnected_components(graph) if len(nodes) >= 3]
    if len(blocks) == 1:
        core = graph.subgraph(blocks[0])
        if core.number_of_edges() == len(core):
            return f"monocyclic-{len(core)}"
        return _fused_chain(core)

    # A bridge separates linked rings; an articulation atom joins spiro rings.
    if not blocks or any(
        graph.subgraph(nodes).number_of_edges() != len(nodes) for nodes in blocks
    ):
        return "polycyclic"
    sizes = sorted(map(len, blocks), reverse=True)
    if len(set.union(*blocks)) == sum(sizes):
        kind = "linked"
    elif len(blocks) == 2:
        kind = "spiro"
    else:
        return "polycyclic"
    return kind + "-" + "-".join(map(str, sizes))
