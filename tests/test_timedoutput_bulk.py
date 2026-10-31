# SPDX-License-Identifier: LGPL-3.0-or-later
"""Numeric batching preserves scalar column and atomic timeline contracts."""

from pathlib import Path

import h5py
import numpy as np
import pytest

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator import _timedoutput as timed
from reacnetgenerator._hmmfilter import _HMMFilter
from reacnetgenerator.timedoutput import semantic_manifest, validate_timed_output
from reacnetgenerator.utils import listtobytes


@pytest.mark.parametrize("row_limit,byte_limit", [(7, 1024), (20, 40), (7, 12), (7, 4)])
def test_numeric_batches_match_scalar_flushes(
    tmp_path, monkeypatch, row_limit, byte_limit
):
    """Empty, partial and crossing appends retain order, pending rows and budgets."""
    monkeypatch.setattr(timed, "_BLOCK_ROWS", row_limit)
    monkeypatch.setattr(timed, "_BLOCK_BYTES", byte_limit)
    with h5py.File(tmp_path / "columns.h5", "w") as file:
        scalar = timed._Column(file, "scalar", "<i8")
        bulk = timed._Column(file, "bulk", "<i8")
        scalar.append(-5)
        bulk.append(-5)
        for values in (
            np.asarray([], dtype=np.int64),
            np.asarray([2**63 - 1], dtype=np.uint64),
            np.arange(83, dtype=">i8")[::-2],
            np.asarray([-(2**63), 0, 7], dtype=np.int64),
        ):
            for value in values:
                scalar.append(int(value))
            bulk.extend_numeric(values)
            assert scalar.rows == bulk.rows
            assert scalar.bytes == bulk.bytes
            np.testing.assert_array_equal(scalar.dataset[:], bulk.dataset[:])
        bulk.append(99)
        scalar.append(99)
        scalar.flush()
        bulk.flush()
        np.testing.assert_array_equal(scalar.dataset[:], bulk.dataset[:])
        assert bulk.dataset.dtype.str == "<i8"


def test_large_array_converts_only_buffer_slices(tmp_path, monkeypatch):
    """A large caller-owned array must never become one trajectory-sized list."""
    lengths = []

    class GuardedArray(np.ndarray):
        def tolist(self):
            assert len(self) <= timed._BLOCK_ROWS
            lengths.append(len(self))
            return super().tolist()

    values = np.arange(2**20 + 3, dtype=np.int64).view(GuardedArray)
    peaks = []
    flush = timed._Column.flush

    def checked(column):
        peaks.append((len(column.rows), column.bytes))
        assert len(column.rows) <= timed._BLOCK_ROWS
        assert column.bytes <= timed._BLOCK_BYTES
        return flush(column)

    monkeypatch.setattr(timed._Column, "flush", checked)
    with h5py.File(tmp_path / "large.h5", "w") as file:
        column = timed._Column(file, "values", "<i8")
        column.extend_numeric(values)
        assert len(column.rows) == 3
        column.flush()
        np.testing.assert_array_equal(column.dataset[:], values)
    assert len(lengths) > 100 and max(lengths) == timed._BLOCK_ROWS
    assert max(rows for rows, _ in peaks) == timed._BLOCK_ROWS


@pytest.mark.parametrize(
    "values,error",
    [
        ([1, 2], TypeError),
        (np.asarray([1.0]), TypeError),
        (np.asarray([np.nan]), TypeError),
        (np.asarray([np.inf]), TypeError),
        (np.asarray([True]), TypeError),
        (np.asarray([1 + 0j]), TypeError),
        (np.asarray(["1"]), TypeError),
        (np.asarray([1], dtype=object), TypeError),
        (np.asarray(1), ValueError),
        (np.ones((2, 2), dtype=np.int64), ValueError),
        (np.asarray([2**63], dtype=np.uint64), OverflowError),
        (np.asarray([2**64 - 1], dtype=np.uint64), OverflowError),
    ],
)
def test_numeric_batches_reject_invalid_values(tmp_path, values, error):
    """Do not silently truncate floats, wrap unsigned integers or flatten shapes."""
    with h5py.File(tmp_path / "invalid.h5", "w") as file:
        column = timed._Column(file, "values", "<i8")
        with pytest.raises(error):
            column.extend_numeric(values)
        assert not column.rows and len(column.dataset) == 0


@pytest.mark.parametrize("dtype", ["<f8", ">i8", h5py.string_dtype("utf-8")])
def test_numeric_batches_reject_other_columns(tmp_path, dtype):
    """The internal numeric API cannot alter the schema or string encoding path."""
    with h5py.File(tmp_path / "dtype.h5", "w") as file:
        column = timed._Column(file, "values", dtype)
        with pytest.raises(TypeError, match="signed int64"):
            column.extend_numeric(np.arange(3))


def _fixture(folder, hmm=False, bad_bonds=None):
    """Create valid ordered molecule records with long structures and many ranges."""
    folder.mkdir()
    source = folder / "source.bond"
    source.touch()
    rng = ReacNetGenerator(
        inputfilename=source,
        inputfiletype="bond",
        atomname=["C"],
        nproc=1,
        runHMM=False,
        output_dir=folder,
        timed_output=folder / "timeline.h5",
    )
    rng.N, rng.step, rng.temp1it = 43, 129, 2
    rng.atomtype = np.zeros(rng.N, dtype=np.int64)
    rng.timestep = dict(enumerate(range(rng.step)))
    rng.source_frame_counts = [rng.step]
    rng.moleculetempfilename = str(folder / "records.bin")
    with open(rng.moleculetempfilename, "wb") as file:
        for atoms in (np.arange(2), np.arange(2, rng.N)):
            pairs = np.column_stack((atoms[:-1], atoms[1:]))
            levels = np.ones(len(pairs), dtype=np.int64)
            if len(pairs) > 32:
                if bad_bonds == "length":
                    levels = levels[:-1]
                elif bad_bonds == "shape":
                    pairs = np.column_stack((pairs, pairs[:, 0]))
            for values in (atoms, pairs, levels, np.arange(0, rng.step, 2)):
                file.write(listtobytes(values))
    _HMMFilter(rng).filter()
    if hmm:
        rng.runHMM = True
        rng.parameters["runHMM"] = True
        rng.hmmfilename = str(folder / "hmm.bin")
        with open(rng.hmmfilename, "wb") as file:
            for _ in range(rng.hmmit):
                file.write(listtobytes((np.arange(rng.step) % 4 != 2).astype(np.int8)))
    return rng


@pytest.mark.parametrize("hmm", [False, True])
def test_batched_writer_matches_scalar_manifest(tmp_path, monkeypatch, hmm):
    """All offsets, bonds and maximal effective ranges match across flush boundaries."""
    monkeypatch.setattr(timed, "_BLOCK_ROWS", 17)
    monkeypatch.setattr(timed, "_BLOCK_BYTES", 1024)
    manifests = []
    for variant, threshold in (("scalar", 2**63), ("bulk", 4)):
        monkeypatch.setattr(timed, "_NUMERIC_MIN_ROWS", threshold)
        rng = _fixture(tmp_path / variant, hmm)
        with timed._TimedOutputWriter(rng) as writer:
            writer.write_molecules(["[C][C]", "C-chain"])
        assert validate_timed_output(rng.timed_output).molecules == 2
        manifests.append(semantic_manifest(rng.timed_output))
        with h5py.File(rng.timed_output) as file:
            assert list(file["molecules/atom_offsets"]) == [0, 2, 43]
            assert list(file["molecules/bond_offsets"]) == [0, 1, 41]
            assert np.all(file["molecules/bond_order"][:] == 1)
            if not hmm:
                assert np.all(file["molecule_ranges/start_frame"][:] % 2 == 0)
                np.testing.assert_array_equal(
                    file["molecule_ranges/start_frame"][:],
                    file["molecule_ranges/end_frame"][:],
                )
    assert manifests[0] == manifests[1]


@pytest.mark.parametrize("failure", ["overflow", "write", "length", "shape"])
def test_bulk_failure_preserves_previous_file(tmp_path, monkeypatch, failure):
    """Invalid inputs and failed HDF5 writes leave only a closed incomplete sibling."""
    rng = _fixture(tmp_path / "work", bad_bonds=failure)
    destination = Path(rng.timed_output)
    destination.write_bytes(b"previous result")
    monkeypatch.setattr(timed, "_BLOCK_ROWS", 17)
    monkeypatch.setattr(timed, "_NUMERIC_MIN_ROWS", 4)
    if failure == "write":
        original = h5py.Dataset.__setitem__

        def fail(dataset, key, values):
            if dataset.name == "/molecules/atom_index" and len(dataset) > 17:
                raise OSError("intentional HDF5 write failure")
            return original(dataset, key, values)

        monkeypatch.setattr(h5py.Dataset, "__setitem__", fail)
    with pytest.raises((OverflowError, OSError, ValueError)):
        with timed._TimedOutputWriter(rng) as writer:
            writer.write_molecules(["[C][C]", "C-chain"])
            if failure == "overflow":
                writer.columns["molecules/atom_index"].extend_numeric(
                    np.asarray([0, 2**63], dtype=np.uint64)
                )
    assert destination.read_bytes() == b"previous result"
    (partial,) = destination.parent.glob("*.incomplete")
    with h5py.File(partial, "r+") as file:
        assert file.attrs["status"] == "incomplete"
    with pytest.raises(ValueError, match="status"):
        validate_timed_output(partial)


def test_ranges_consume_only_one_bounded_chunk(tmp_path):
    """A long interval generator must not be materialized before column writes."""
    consumed = produced = peak = 0
    with h5py.File(tmp_path / "ranges.h5", "w") as file:
        writer = object.__new__(timed._TimedOutputWriter)
        writer.columns = {
            f"molecule_ranges/{name}": timed._Column(file, name, "<i8")
            for name in ("molecule_id", "start_frame", "end_frame")
        }
        end = writer.columns["molecule_ranges/end_frame"]
        append, extend = end.append, end.extend_numeric

        def append_end(value):
            nonlocal consumed
            append(value)
            consumed += 1

        def extend_end(values):
            nonlocal consumed
            extend(values)
            consumed += len(values)

        end.append, end.extend_numeric = append_end, extend_end

        def ranges():
            nonlocal produced, peak
            for value in range(65539):
                produced += 1
                peak = max(peak, produced - consumed)
                assert produced - consumed <= 4096
                yield 2 * value, 2 * value

        writer._write_ranges(7, ranges())
        assert consumed == produced == 65539
        assert peak == 4096
        for column in writer.columns.values():
            column.flush()
        np.testing.assert_array_equal(file["start_frame"][:], np.arange(65539) * 2)
        np.testing.assert_array_equal(file["end_frame"][:], file["start_frame"][:])
        assert np.all(file["molecule_id"][:] == 7)
