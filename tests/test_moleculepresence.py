# SPDX-License-Identifier: LGPL-3.0-or-later
"""Signal reuse, shared-reader alignment, and bounded atom-frame writes."""

import pickle
import time
from functools import partial
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._hmmfilter import _HMMFilter
from reacnetgenerator._moleculepresence import (
    _index_ranges,
    _iter_molecule_presence,
    _iter_molecule_records,
    _MoleculePresence,
    _signal_ranges,
    _use_frame_indices,
)
from reacnetgenerator._path import _CollectSMILESPaths, _write_molecule_presence
from reacnetgenerator._timedoutput import _TimedOutputWriter
from reacnetgenerator.timedoutput import iter_molecule_ranges, validate_timed_output
from reacnetgenerator.utils import (
    bytestolist,
    idx_to_signal,
    listtobytes,
    read_compressed_block,
    run_mp,
)


def _input(tmp_path, frames, step=2048, nproc=1):
    """Build actual internal records, retaining independent frame-list oracles."""
    source = tmp_path / "input.bond"
    source.touch()
    rng = ReacNetGenerator(
        inputfilename=str(source),
        inputfiletype="bond",
        atomname=["H"],
        runHMM=False,
        nproc=nproc,
        output_dir=tmp_path,
        timed_output=tmp_path / "timeline.h5",
    )
    rng.N = max(1, len(frames))
    rng.atomtype = np.zeros(rng.N, dtype=int)
    rng.step = step
    rng.timestep = dict(enumerate(range(step)))
    rng.source_frame_counts = [step]
    rng.temp1it = len(frames)
    definitions = tmp_path / "molecules.bin"
    with definitions.open("wb") as file:
        for atom, indices in enumerate(frames):
            for value in ([atom], [], [], np.asarray(indices, dtype=np.int64)):
                file.write(listtobytes(value))
    rng.moleculetempfilename = str(definitions)
    return rng


def _remove_signals(rng):
    for name in (rng.originfilename, rng.hmmfilename):
        if name:
            Path(name).unlink(missing_ok=True)


def test_sparse_nohmm_does_not_expand_or_start_workers(tmp_path, monkeypatch):
    """No-HMM reuses definitions even when callers do not request origin output."""
    frames = [[], [0, 1000, 2047], [7, 2, 7, 3]]
    rng = _input(tmp_path, frames, nproc=2)
    rng.getoriginfile = False
    original = Path(rng.moleculetempfilename).read_bytes()

    def unexpected(*args, **kwargs):
        raise AssertionError(
            "Sparse pass-through must not expand signals or use a pool"
        )

    monkeypatch.setattr("reacnetgenerator._hmmfilter.run_mp", unexpected)
    monkeypatch.setattr("reacnetgenerator._hmmfilter.idx_to_signal", unexpected)
    _HMMFilter(rng).filter()
    try:
        assert rng.moleculetemp2filename == rng.moleculetempfilename
        assert Path(rng.moleculetempfilename).read_bytes() == original
        assert Path(rng.originfilename).stat().st_size == 0
        assert rng.hmmfilename is None
        assert rng.hmmit == len(frames)
        actual = [
            list(presence.ranges()) for _, presence in _iter_molecule_presence(rng)
        ]
        expected = [
            list(_signal_ranges(idx_to_signal(values, rng.step))) for values in frames
        ]
        assert actual == expected
    finally:
        _remove_signals(rng)


def _delayed_origin(record, *, transform, step, completions):
    """Force real workers to finish dense records in the opposite order."""
    count = len(bytestolist(record[0]))
    if count == step:
        time.sleep(1)
    result = transform(record)
    if result is not None:
        with open(completions, "a") as file:
            file.write(f"{count}\n")
    return result


@pytest.mark.parametrize("nproc", [1, 2])
def test_mixed_presence_keeps_matrix_and_hdf5_alignment(tmp_path, monkeypatch, nproc):
    """Sparse, dense, empty, and delayed records share one consistent ID order."""
    step = 2048
    frames = [[1, 5], np.arange(step), [], np.arange(1, step, 2), [step - 1]]
    rng = _input(tmp_path, frames, step, nproc)
    completions = tmp_path / "completed"

    def delayed_workers(nproc, **kwargs):
        # Worker setup must stay small even though rng retains frame metadata.
        assert len(pickle.dumps(kwargs["func"])) < 1024
        kwargs["func"] = partial(
            _delayed_origin,
            transform=kwargs["func"],
            step=step,
            completions=completions,
        )
        return run_mp(nproc, **kwargs)

    monkeypatch.setattr("reacnetgenerator._hmmfilter.run_mp", delayed_workers)
    _HMMFilter(rng).filter()
    try:
        with open(rng.originfilename, "rb") as signals:
            assert len(list(read_compressed_block(signals))) == 2
        if nproc == 2:
            assert completions.read_text().splitlines() == [
                "1024",
                "2048",
            ]
        expected = np.zeros((rng.N, step), dtype=np.uint8)
        for atom, indices in enumerate(frames):
            expected[atom, indices] = atom + 1
        collector = _CollectSMILESPaths(rng)
        with collector._getatomeach() as store:
            np.testing.assert_array_equal(store.atomeach, expected)
            assert not np.any(store.conflict)
        with _TimedOutputWriter(rng) as writer:
            writer.write_molecules(["[H]"] * rng.hmmit)
        actual_ranges = [
            (r.molecule_id, r.start_frame, r.end_frame)
            for r in iter_molecule_ranges(rng.timed_output)
        ]
        expected_ranges = [
            (atom + 1, start, stop)
            for atom, indices in enumerate(frames)
            for start, stop in _signal_ranges(idx_to_signal(indices, step))
        ]
        assert actual_ranges == expected_ranges
        assert validate_timed_output(rng.timed_output).frames == step
        assert Path(rng.moleculetempfilename).exists()
        assert not list(tmp_path.glob("reacnetgenerator-*.mmap"))
    finally:
        _remove_signals(rng)


@pytest.mark.parametrize("frames", [[], [0], [0, 1, 2, 3, 7, 8, 9, 12]])
@pytest.mark.parametrize("block_rows", [1, 2, 3])
def test_index_ranges_join_across_blocks(frames, block_rows):
    """Sparse maximal intervals agree with an independently expanded signal."""
    expected = list(_signal_ranges(idx_to_signal(frames, 15), block_rows=4))
    assert list(_index_ranges(np.asarray(frames), block_rows)) == expected


@pytest.mark.parametrize("indexed", [False, True])
def test_matrix_writes_are_bounded_and_preserve_conflicts(monkeypatch, indexed):
    """Splitting either dimension must preserve last-writer and conflict cells."""
    monkeypatch.setattr("reacnetgenerator._moleculepresence._PRESENCE_BLOCK_ROWS", 7)
    monkeypatch.setattr("reacnetgenerator._path._MATRIX_WRITE_CELLS", 12)
    values = np.zeros((19, 31), dtype=np.uint16)
    conflicts = np.zeros_like(values, dtype=bool)
    reads = []

    class Matrix:
        def __getitem__(self, index):
            selected = values[index]
            reads.append(selected.size)
            assert selected.size <= 12
            return selected

        def __setitem__(self, index, value):
            values[index] = value

    class Conflicts:
        def mark(self, atoms, frames, overlap):
            assert overlap.size <= 12
            conflicts[np.ix_(atoms, frames)] |= overlap

    store = SimpleNamespace(atomeach=Matrix(), conflict=Conflicts())
    expected = np.zeros_like(values)
    expected_conflicts = np.zeros_like(conflicts)
    for molecule_id, atoms, frames in (
        (1, list(range(19)), list(range(31))),
        (2, [18, 0, 7, 7, 0], [0, 2, 2, 11, 12, 30]),
        (3, [1, 4, 7], []),
    ):
        signal = idx_to_signal(frames, 31)
        presence = _MoleculePresence(frames if indexed else signal, 31, indexed=indexed)
        _write_molecule_presence(store, atoms, presence, molecule_id)
        index = np.ix_(np.unique(atoms), np.flatnonzero(signal))
        expected_conflicts[index] |= expected[index] != 0
        expected[index] = molecule_id
    np.testing.assert_array_equal(values, expected)
    np.testing.assert_array_equal(conflicts, expected_conflicts)
    assert len(reads) > 1


@pytest.mark.parametrize("problem", ["missing", "extra", "short", "few", "many"])
def test_presence_reader_rejects_misaligned_streams(tmp_path, problem):
    """Both consumers must reject missing, extra, or incorrectly sized records."""
    rng = _input(tmp_path, [np.arange(2048)])
    _HMMFilter(rng).filter()
    try:
        if problem == "missing":
            Path(rng.originfilename).write_bytes(b"")
        elif problem == "extra":
            with open(rng.originfilename, "ab") as file:
                file.write(listtobytes(np.ones(rng.step, dtype=bool)))
        elif problem == "short":
            Path(rng.originfilename).write_bytes(listtobytes(np.ones(2, dtype=bool)))
        else:
            rng.hmmit = 0 if problem == "many" else 2
        with pytest.raises(ValueError):
            list(_iter_molecule_presence(rng))
        with pytest.raises(ValueError):
            _CollectSMILESPaths(rng)._getatomeach()
        assert not list(tmp_path.glob("reacnetgenerator-*.mmap"))
    finally:
        _remove_signals(rng)


@pytest.mark.parametrize("cut", [1, 65, 180])
def test_truncated_molecule_records_are_not_silently_skipped(tmp_path, cut):
    """Validate even unselected payloads before accepting a molecule record."""
    rng = _input(tmp_path, [[1, 2]])
    path = Path(rng.moleculetempfilename)
    path.write_bytes(path.read_bytes()[:-cut])
    with (
        path.open("rb") as file,
        pytest.raises(ValueError, match="Incomplete molecule"),
    ):
        list(_iter_molecule_records(file, fields=(0,)))


def test_empty_nohmm_input(tmp_path):
    """An empty definition stream retains valid count and shared-file ownership."""
    rng = _input(tmp_path, [])
    _HMMFilter(rng).filter()
    try:
        assert rng.hmmit == 0
        assert rng.moleculetemp2filename == rng.moleculetempfilename
        assert list(_iter_molecule_presence(rng)) == []
    finally:
        _remove_signals(rng)


@pytest.mark.parametrize("nproc", [1, 2])
def test_failed_dense_expansion_preserves_shared_definitions(tmp_path, nproc):
    """Worker or local failure closes the stream without deleting its input."""
    rng = _input(tmp_path, [np.full(2048, 2048)], nproc=nproc)
    hmm_filter = _HMMFilter(rng)
    with pytest.raises(IndexError):
        hmm_filter.filter()
    assert Path(rng.moleculetempfilename).exists()
    assert not Path(hmm_filter.originfilename).exists()


def test_direct_selection_uses_decoded_size():
    """Highly compressible dense indices must still take the bounded signal path."""
    assert _use_frame_indices(listtobytes(np.asarray([1, 9999])), 10000)
    assert not _use_frame_indices(listtobytes(np.arange(10000)), 10000)


@pytest.mark.parametrize("frames", [[-1, 2, 2, 0], [2, 0, 1], []])
def test_index_normalization_matches_original_signal(frames):
    """Duplicates and NumPy negative indices do not create duplicate intervals."""
    actual = list(_MoleculePresence(frames, 4, indexed=True).ranges())
    assert actual == list(_signal_ranges(idx_to_signal(frames, 4)))
