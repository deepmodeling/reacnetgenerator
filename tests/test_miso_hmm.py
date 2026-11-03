# SPDX-License-Identifier: LGPL-3.0-or-later
"""Regressions for pooling isomer-equivalent observations before Viterbi decoding."""

import csv
from pathlib import Path

import numpy as np
import pytest

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._hmmfilter import _filter_atom_set
from reacnetgenerator.timedoutput import (
    iter_molecule_ranges,
    iter_molecules,
    validate_timed_output,
)
from reacnetgenerator.utils import bytestolist, listtobytes


@pytest.mark.parametrize("nproc", [1, 2])
def test_alternating_oxygen_survives_hmm(tmp_path, monkeypatch, nproc):
    """Keep both observed assignments and one continuous route per physical pair."""
    if nproc == 2:
        # Exercise route-ID propagation to workers even on this 800 atom-frame input.
        monkeypatch.setattr("reacnetgenerator._path._STEP3_SERIAL_ATOM_FRAMES", 0)
    source = tmp_path / "oxygen.dump"
    distances = ([1.23529] * 10 + [1.22136] * 10) * 10
    with source.open("w") as stream:
        for t, d in enumerate(distances):
            stream.write(
                f"ITEM: TIMESTEP\n{t}\nITEM: NUMBER OF ATOMS\n4\n"
                "ITEM: BOX BOUNDS pp pp pp\n0 20\n0 20\n0 20\n"
                f"ITEM: ATOMS id type x y z\n1 1 5 5 5\n2 1 {5 + d} 5 5\n"
                "3 1 12 5 5\n4 1 13.23529 5 5\n"
            )
    path = tmp_path / "timeline.h5"
    rng = ReacNetGenerator(
        inputfilename=source,
        inputfiletype="dump",
        atomname=["O"],
        nproc=nproc,
        miso=1,
        output_dir=tmp_path,
        printmoleculetime=True,
        timed_output=path,
    )
    rng.run_items(["species", "reactions"])
    populations = Path(rng.speciesfilename).read_text().splitlines()
    assert len(populations) == 200
    assert all(line.split(":")[1].split() == ["[O][O]", "2"] for line in populations)
    assert Path(rng.atomroutefilename).read_text().count(" -> ") == 0
    assert Path(rng.reactionabcdfilename).read_text() == ""
    assert validate_timed_output(path).reaction_events == 0
    molecules = {m.molecule_id: m for m in iter_molecules(path)}
    ranges = [
        (molecules[r.molecule_id].bonds, r.start_frame, r.end_frame)
        for r in iter_molecule_ranges(path)
    ]
    assert sorted(ranges) == sorted(
        [(((0, 1, 1 + t // 10 % 2),), t, t + 9) for t in range(0, 200, 10)]
        + [(((2, 3, 1),), 0, 199)]
    )
    with open(rng.moleculetimelinefilename, newline="") as stream:
        observed = list(csv.DictReader(stream))
    assert len(observed) == 400
    assert len({row["BondIDs"] for row in observed}) == 3


class _Model:
    """Capture the pre-HMM signal and return a known path for projection."""

    def __init__(self, path):
        self.inputs = []
        self.path = np.array(path, dtype=bool)

    def predict(self, observations):
        self.inputs.append(observations.ravel().tolist())
        return self.path


@pytest.mark.parametrize("miso", [1, 2])
def test_vf2_classes_and_gap_projection(miso):
    """Pool relabeled graphs; keep a distinct labeled topology separate in mode 1."""
    atoms = (0, 1, 2)
    types = [0, 0, 1]
    # Two C-centered C-C-O chains, then a C-O-C chain on the same atoms.
    pairs = [[(0, 1), (0, 2)], [(0, 1), (1, 2)], [(0, 2), (1, 2)]]
    visits = [[1], [3], [5]]
    records = [
        tuple(listtobytes(v) for v in (atoms, edges, [1, 1], frames))
        for edges, frames in zip(pairs, visits, strict=True)
    ]
    model = _Model([1, 1, 1, 1, 0, 1, 0])
    results = _filter_atom_set(
        (atoms, types, records),
        miso=miso,
        step=7,
        model=model,
        getoriginfile=True,
        printfiltersignal=False,
    )
    if miso == 1:
        assert model.inputs == [[0, 1, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0]]
        assert [row[0] for row in results] == [0, 0, 2]
        expected = [[0, 1, 2], [3, 5], [0, 1, 2, 3, 5]]
    else:
        assert model.inputs == [[0, 1, 0, 1, 0, 1, 0]]
        assert [row[0] for row in results] == [0, 0, 0]
        expected = [[0, 1, 2], [3], [5]]
    for i, (_, origin, signal, record) in enumerate(results):
        assert record == records[i]
        assert np.flatnonzero(bytestolist(origin)).tolist() == visits[i]
        assert np.flatnonzero(bytestolist(signal)).tolist() == expected[i]


@pytest.mark.parametrize("keep_empty", [False, True])
def test_filter_zero_pooled_signal(keep_empty):
    """Optional filtered-signal output retains aligned records even for an empty path."""
    atoms = (0, 1)
    records = [
        tuple(listtobytes(v) for v in (atoms, [(0, 1)], [order], [order - 1]))
        for order in (1, 2)
    ]
    model = _Model([0, 0])
    results = _filter_atom_set(
        (atoms, [0, 0], records),
        miso=1,
        step=2,
        model=model,
        getoriginfile=False,
        printfiltersignal=keep_empty,
    )
    assert model.inputs == [[1, 1]]
    assert len(results) == (2 if keep_empty else 0)
    for _, origin, signal, _ in results:
        assert origin is None
        assert not bytestolist(signal).any()
