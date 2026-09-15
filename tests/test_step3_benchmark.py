# SPDX-License-Identifier: LGPL-3.0-or-later
"""Small deterministic Step 3 kernels using the existing CodSpeed harness."""

import numpy as np
import pytest

from reacnetgenerator._path import _CollectSMILESPaths
from reacnetgenerator._reaction import ReactionsFinder
from reacnetgenerator.utils import listtobytes


@pytest.mark.benchmark
def test_benchmark_atom_route(benchmark):
    """Measure route compression at the same seam available on the base commit."""
    collector = object.__new__(_CollectSMILESPaths)
    collector.atomname = np.array(["H"])
    collector.selectatoms = ["H"]
    collector.mname = np.array(["A", "B"])
    timeline = np.repeat(np.array([1, 2, 0, 2, 1], dtype=np.uint8), 20000)
    item = (1, (timeline, 0))
    assert collector._getatomroute(item)[1] == ("Atom 1 H: 0 A -> 20000 B -> 60000 A")

    @benchmark
    def bench():
        collector._getatomroute(item)


@pytest.mark.benchmark
def test_benchmark_transition_graph(benchmark):
    """Measure bounded graph scanning with a deterministic connected reaction."""
    finder = object.__new__(ReactionsFinder)
    finder.mname = np.array(["A", "B", "C"])
    finder.printreactionevent = False
    before = np.repeat(np.array([1, 2], dtype=np.uint16), 2048)
    after = np.full(4096, 3, dtype=np.uint16)
    payload = listtobytes(
        (before, after, np.zeros(4096, dtype=bool), np.zeros(4096, dtype=bool))
    )
    assert finder._getstepreaction(payload) == ["A+B->C"]

    @benchmark
    def bench():
        finder._getstepreaction(payload)
