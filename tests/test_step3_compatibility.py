# SPDX-License-Identifier: LGPL-3.0-or-later
"""Regression tests for Step 3 signal alignment and scratch storage."""

import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pytest

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._path import _CollectSMILESPaths
from reacnetgenerator.utils import listtobytes


def test_hmm_origin_output_preserves_reactions(tmp_path):
    """Keeping origin signals must not misalign filtered molecule records."""
    results = []
    molecule_counts = []
    for getoriginfile in (False, True):
        directory = tmp_path / str(getoriginfile)
        directory.mkdir()
        input_path = directory / "reaction.bond"
        shutil.copyfile(Path(__file__).parent / "inputs/reaction.bond", input_path)
        rng = ReacNetGenerator(
            inputfilename=str(input_path),
            inputfiletype="lammpsbondfile",
            atomname=["H", "O"],
            nproc=1,
            runHMM=True,
            getoriginfile=getoriginfile,
            printfiltersignal=False,
        )
        rng.run()
        molecule_counts.append(rng.hmmit)
        results.append(
            (
                Path(rng.atomroutefilename).read_bytes(),
                Path(rng.reactionfilename).read_bytes(),
                Path(rng.reactionabcdfilename).read_bytes(),
            )
        )
        assert not list(directory.glob("reacnetgenerator-*.mmap"))

    # The fixture must include filtered molecules to exercise missing records.
    assert molecule_counts[1] > molecule_counts[0]
    assert results[1] == results[0]


@pytest.mark.parametrize("failed_allocation", [1, 2])
def test_discarded_atom_routes_use_scratch_fallback(
    tmp_path, monkeypatch, failed_allocation
):
    """Unwritable output parents must allow fallback without leaking files."""
    molecule_file = tmp_path / "molecules.bin"
    molecule_file.write_bytes(
        b"".join(listtobytes(value) for value in ([0], [], [], []))
    )
    hmm_file = tmp_path / "hmm.bin"
    hmm_file.write_bytes(listtobytes(np.array([True, False])))

    collector = object.__new__(_CollectSMILESPaths)
    collector.N = 1
    collector.step = 2
    collector.hmmit = 1
    collector.runHMM = True
    collector.hmmfilename = str(hmm_file)
    collector.moleculetemp2filename = str(molecule_file)
    collector.atomroutefilename = os.devnull

    preferred = tmp_path / "preferred"
    fallback = tmp_path / "fallback"
    preferred.mkdir()
    fallback.mkdir()
    original_mkstemp = tempfile.mkstemp
    preferred_attempts = 0

    def allocate(*args, **kwargs):
        nonlocal preferred_attempts
        if kwargs.get("dir") is not None:
            assert kwargs["dir"] == os.path.dirname(os.path.abspath(os.devnull))
            preferred_attempts += 1
            if preferred_attempts == failed_allocation:
                raise PermissionError("simulated unwritable output parent")
            kwargs["dir"] = preferred
        else:
            kwargs["dir"] = fallback
        return original_mkstemp(*args, **kwargs)

    monkeypatch.setattr("reacnetgenerator._step3state.tempfile.mkstemp", allocate)

    with collector._getatomeach() as store:
        np.testing.assert_array_equal(store.atomeach, [[1, 0]])
        np.testing.assert_array_equal(store.conflict, [[False, False]])
        assert Path(store.atomeach_path).parent == fallback
        assert Path(store.conflict.path).parent == fallback
        assert not list(preferred.iterdir())

    assert preferred_attempts == failed_allocation
    assert not list(fallback.iterdir())
