# SPDX-License-Identifier: LGPL-3.0-or-later
"""LAMMPS dump array placement and invalid-input regressions."""

from io import StringIO

import numpy as np
import pytest
from ase import Atom, Atoms

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._detect import _DetectLAMMPSdump


def dump_frame(rows, box=("0 10", "0 10", "0 10"), count=3):
    """Use reordered columns and atom rows, as real LAMMPS dumps permit."""
    return [
        "ITEM: TIMESTEP",
        "40",
        "ITEM: NUMBER OF ATOMS",
        str(count),
        "ITEM: BOX BOUNDS pp pp pp",
        *box,
        "ITEM: ATOMS x y z type q id",
        *rows,
    ]


ROWS = ["1 2 3 2 0 3", "2 2 3 1 0 1", "1 3 3 1 0 2"]


def detector(lines):
    """Initialize the parser through its real first-frame setup."""
    result = _DetectLAMMPSdump(
        ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename="unused",
            atomname=["H", "O"],
            pbc=False,
        )
    )
    result._readNfunc(StringIO("\n".join(lines) + "\n"))
    return result


@pytest.mark.parametrize(
    "box, expected_cell",
    [
        (("-2 8", "1 11", "0 10"), [[10, 0, 0], [0, 10, 0], [0, 0, 10]]),
        (("-2 11 2", "1 12 -1", "0 10 1"), [[10, 0, 0], [2, 10, 0], [-1, 1, 10]]),
    ],
)
@pytest.mark.parametrize("reverse", [False, True])
def test_atom_order_and_cell(box, expected_cell, reverse, monkeypatch):
    """Preserve canonical atom order, unshifted coordinates and tilted cells."""
    lines = dump_frame(ROWS[::-1] if reverse else ROWS, box)
    parser = detector(lines)
    np.testing.assert_array_equal(parser.atomtype, [0, 0, 1])
    original = parser._getbondfromcrd
    captured = []

    def capture(atoms, cell):
        captured.append((atoms.copy(), cell.copy()))
        return original(atoms, cell)

    monkeypatch.setattr(parser, "_getbondfromcrd", capture)
    actual, timestep = parser._readstepfunc((7, lines))
    atoms, cell = captured[0]
    np.testing.assert_array_equal(atoms.numbers, [1, 1, 8])
    np.testing.assert_array_equal(atoms.positions, [[2, 2, 3], [1, 3, 3], [1, 2, 3]])
    np.testing.assert_array_equal(cell, expected_cell)
    # Independent legacy construction verifies records through real bond detection.
    legacy = Atoms([Atom("H", (2, 2, 3)), Atom("H", (1, 3, 3)), Atom("O", (1, 2, 3))])
    bond, level = original(legacy, np.array(expected_cell, dtype=float))
    assert actual == parser._connectmolecule(bond, level, frame=7, timestep=40)
    assert timestep == (7, 40)


@pytest.mark.parametrize(
    "rows, message",
    [
        ([ROWS[0], ROWS[1], ROWS[1]], "duplicate atom ID"),
        (ROWS[:2], "missing one or more atom IDs"),
        ([ROWS[0], ROWS[1], "1 3 3 1 0 0"], "atom ID is out of range"),
        ([ROWS[0], ROWS[1], "1 3 3 1 0 4"], "atom ID is out of range"),
        ([ROWS[0], ROWS[1], "1 3 3 0 0 2"], "atom type is out of range"),
        ([ROWS[0], ROWS[1], "1 3 3 3 0 2"], "atom type is out of range"),
    ],
)
@pytest.mark.parametrize("first_frame", [False, True])
def test_invalid_atom_mapping(rows, message, first_frame):
    """Reject bad metadata before array placement in both parsing passes."""
    bad = dump_frame(rows)
    with pytest.raises(ValueError, match=message):
        if first_frame:
            detector(bad)
        else:
            detector(dump_frame(ROWS))._readstepfunc((1, bad))


def test_first_frame_of_multiframe_dump():
    """The setup pass stops after the next atom-count header, preserving stride."""
    lines = dump_frame(ROWS)
    parser = detector(lines)
    stride = parser._readNfunc(StringIO("\n".join(lines + lines) + "\n"))
    assert stride == len(lines)
    np.testing.assert_array_equal(parser.atomtype, [0, 0, 1])


def test_benchmark_lammps_dump_parser(benchmark):
    """Exercise parsing plus real detection on a deterministic shuffled frame."""
    rng = np.random.default_rng(2026)
    positions = rng.uniform(0, 40, (2048, 3))
    rows = [f"{x} {y} {z} 1 0 {i + 1}" for i, (x, y, z) in enumerate(positions)]
    rng.shuffle(rows)
    lines = dump_frame(rows, ("0 40",) * 3, count=len(rows))
    parser = detector(lines)
    expected = parser._readstepfunc((0, lines))

    @benchmark
    def run():
        assert parser._readstepfunc((0, lines)) == expected
