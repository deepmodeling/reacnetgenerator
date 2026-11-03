# SPDX-License-Identifier: LGPL-3.0-or-later
"""Order-independent overlap resolution and reproducible HMM artifacts."""

import itertools
import multiprocessing
from contextlib import nullcontext
from pathlib import Path

import numpy as np
import pytest

from reacnetgenerator import ReacNetGenerator, _detect, _moleculeorder
from reacnetgenerator._moleculeorder import _OrderedMoleculeResults, _structure_key
from reacnetgenerator._moleculepresence import _MoleculePresence
from reacnetgenerator._path import _route_changes, _write_molecule_presence
from reacnetgenerator._step3state import _AtomFrameStore
from reacnetgenerator.timedoutput import semantic_manifest
from reacnetgenerator.utils import listtobytes
from tests.test_hmm_scheduling import _use_hmm_context


def _record(atoms, pairs, orders, frames=(0, 1)):
    """Encode a detection record without relying on the sorting implementation."""
    return tuple(listtobytes(value) for value in (atoms, pairs, orders, frames))


def test_structural_order_is_numeric_and_includes_bond_orders(tmp_path, monkeypatch):
    """Compressed representation, tuple prefixes and byte boundaries are distinct."""
    monkeypatch.setattr(_moleculeorder.tempfile, "tempdir", str(tmp_path))
    records = [
        _record([1], [], []),
        _record([1, 2], [(1, 2)], [1]),
        _record([1, 2], [(1, 2)], [2]),
        _record([1, 2, 3], [(1, 2), (2, 3)], [1, 1]),
        _record([255], [], []),
        _record([256], [], []),
        _record([65536], [], []),
    ]
    alternative = tuple(
        listtobytes(np.asarray(value)) for value in ([1, 2], [(1, 2)], [1])
    )
    assert alternative != records[1][:3]
    assert _structure_key(alternative) == _structure_key(records[1])
    for permutation in (list(reversed(range(7))), [3, 0, 5, 2, 6, 4, 1]):
        with _OrderedMoleculeResults() as ordered:
            for index in permutation:
                ordered.add(
                    None if index % 2 else b"origin", bytes([index]), records[index]
                )
            result = list(ordered)
        assert [item[2] for item in result] == records
        assert [item[1] for item in result] == [bytes([i]) for i in range(7)]
        assert [item[0] for item in result] == [
            None if i % 2 else b"origin" for i in range(7)
        ]
        assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("failure", [False, True])
def test_order_index_spills_and_cleans_up(tmp_path, monkeypatch, failure):
    """Large result streams spill before consumption and release scratch on errors."""
    monkeypatch.setattr(_moleculeorder.tempfile, "tempdir", str(tmp_path))
    monkeypatch.setattr(_moleculeorder, "_ORDER_CACHE_KIB", 8)
    with pytest.raises(RuntimeError, match="consumer") if failure else nullcontext():
        with _OrderedMoleculeResults() as ordered:
            for i in reversed(range(4096)):
                ordered.add(None, b"x" * 4096, _record([i], [], []))
            directory = next(tmp_path.iterdir())
            assert (directory / "results.bin").stat().st_size >= 16 * 1024 * 1024
            assert (directory / "index.sqlite").stat().st_size > 8192
            count = 0
            for count, result in enumerate(ordered, 1):
                assert result[2] == _record([count - 1], [], [])
                if failure:
                    raise RuntimeError("consumer failed")
            assert count == 4096
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("failure", ["open", "append"])
def test_order_spool_cleans_up_producer_failures(tmp_path, monkeypatch, failure):
    """Opening the index or serializing a result can fail without leaking files."""
    monkeypatch.setattr(_moleculeorder.tempfile, "tempdir", str(tmp_path))
    record = _record([0], [], [])

    def fail(*args, **kwargs):
        raise OSError("spool failure")

    if failure == "open":
        monkeypatch.setattr(_moleculeorder.sqlite3, "connect", fail)
    else:
        monkeypatch.setattr(_moleculeorder.pickle, "dump", fail)
    with (
        pytest.raises(OSError, match="spool failure"),
        _OrderedMoleculeResults() as ordered,
    ):
        ordered.add(None, b"signal", record)
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("indexed", [False, True])
def test_overlap_winner_and_conflicts_ignore_write_order(
    tmp_path, monkeypatch, indexed
):
    """All permutations agree for identical and partially overlapping atom sets."""
    monkeypatch.setattr("reacnetgenerator._path._MATRIX_WRITE_CELLS", 4)
    monkeypatch.setattr("reacnetgenerator._moleculepresence._PRESENCE_BLOCK_ROWS", 3)
    records = [([0, 1], range(6)), ([0, 1], range(3, 10)), ([1, 2], range(5, 8))]
    candidates = np.zeros((3, 3, 10), dtype=bool)
    for i, (atoms, frames) in enumerate(records):
        candidates[i][np.ix_(atoms, frames)] = True
    expected = np.max(candidates * np.arange(1, 4)[:, None, None], axis=0)
    conflicts = candidates.sum(axis=0) > 1
    for permutation in itertools.permutations(range(3)):
        with _AtomFrameStore((3, 10), 3, directory=tmp_path) as store:
            for i in permutation:
                atoms, frames = records[i]
                values = (
                    np.asarray(frames) if indexed else np.isin(np.arange(10), frames)
                )
                _write_molecule_presence(
                    store, atoms, _MoleculePresence(values, 10, indexed=indexed), i + 1
                )
            np.testing.assert_array_equal(store.atomeach, expected)
            np.testing.assert_array_equal(store.conflict, conflicts)
            positions, molecules = _route_changes(store.atomeach[0])
            assert positions.tolist() == [0, 3]
            assert molecules.tolist() == [1, 2]


def test_repeated_molecule_write_is_not_a_conflict(tmp_path):
    """Writing the same presence twice must not create an ambiguous cell."""
    with _AtomFrameStore((1, 10), 1, directory=tmp_path) as store:
        presence = _MoleculePresence(np.arange(10), 10, indexed=True)
        for _ in range(2):
            _write_molecule_presence(store, [0], presence, 1)
        assert not np.any(store.conflict)


def _write_overlapping_input(path):
    """Two C pairs provide both a genuine HMM overlap and a clean bond event."""
    # Three mutually exclusive bond states on the first pair. Default HMM
    # yields single-bond presence [0, 234] and double-bond presence [225, 299].
    runs = [
        (0, 17),
        (2, 6),
        (0, 27),
        (1, 3),
        (0, 29),
        (1, 6),
        (2, 11),
        (1, 6),
        (2, 6),
        (0, 10),
        (2, 4),
        (0, 24),
        (1, 8),
        (0, 5),
        (2, 4),
        (0, 5),
        (1, 2),
        (0, 17),
        (1, 3),
        (0, 11),
        (2, 6),
        (0, 15),
        (1, 5),
        (0, 5),
        (1, 14),
        (2, 6),
        (1, 37),
        (2, 8),
    ]
    states = [state for state, length in runs for _ in range(length)]
    with path.open("w") as file:
        for frame, state in enumerate(states):
            file.write(f"# Timestep {frame}\n#\n# Number of particles 4\n#\n")
            for first, order in ((1, state + 1), (3, 1 if frame < 150 else 2)):
                for atom, neighbor in ((first, first + 1), (first + 1, first)):
                    file.write(f"{atom} 1 1 {neighbor} 1 {order} {order} 0 0\n")
            file.write("#\n")


def _artifacts(source, output, nproc):
    """Run the real pipeline and compare every semantic HDF5 dataset."""
    rng = ReacNetGenerator(
        inputfilename=source,
        inputfiletype="bond",
        atomname=["C"],
        nproc=nproc,
        miso=0,
        runHMM=True,
        needprintspecies=True,
        output_dir=output,
        timed_output=output / "timeline.h5",
    )
    rng.run()
    result = {
        field: Path(getattr(rng, field)).read_bytes()
        for field in (
            "atomroutefilename",
            "speciesfilename",
            "reactionfilename",
            "tablefilename",
        )
    }
    manifest = semantic_manifest(rng.timed_output)
    assert manifest["datasets"]["reaction_events/transition"]["rows"] > 0
    assert manifest["datasets"]["transition_evidence/transition"]["rows"] > 0
    manifest["metadata"]["configuration"]["parameters"].pop("nproc")
    manifest.pop("semantic_sha256")
    result["manifest"] = manifest
    assert result["atomroutefilename"].decode().splitlines()[:2] == [
        "Atom 1 C: 0 [C][C] -> 225 [C]=[C]",
        "Atom 2 C: 0 [C][C] -> 225 [C]=[C]",
    ]
    return result


@pytest.mark.parametrize("method", ["fork", "spawn"])
@pytest.mark.parametrize("nproc", [1, 2, 8])
def test_pipeline_artifacts_ignore_workers_and_detection_order(
    tmp_path, monkeypatch, method, nproc
):
    """Repeated real pools preserve routes, species and normalized manifests."""
    if method not in multiprocessing.get_all_start_methods():
        pytest.skip(f"{method} unavailable")
    _use_hmm_context(monkeypatch, method)
    monkeypatch.setattr(_moleculeorder.tempfile, "tempdir", str(tmp_path))
    source = tmp_path / "overlap.bond"
    _write_overlapping_input(source)
    expected = _artifacts(source, tmp_path / "serial", 1)
    for repeat in range(2):
        assert _artifacts(source, tmp_path / f"workers-{repeat}", nproc) == expected

    original = _detect.run_mp

    def reversed_detection(workers, **kwargs):
        """Reverse completed frame results, preserving their original indices."""
        yield from reversed(list(original(workers, **kwargs)))

    monkeypatch.setattr(_detect, "run_mp", reversed_detection)
    assert _artifacts(source, tmp_path / "reversed", nproc) == expected
    assert not list(tmp_path.glob("reacnetgenerator-order-*"))
