# SPDX-License-Identifier: LGPL-3.0-or-later
"""Bounded HDF5 presence lookups preserve the scalar search contract."""

from operator import getitem
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from reacnetgenerator import _timedoutputrange as ranges


class _CountedColumn:
    """Count reads and reject unbounded or scalar reads from the cached path."""

    def __init__(self, data):
        self.data = data
        self.reads = 0

    def __len__(self):
        return len(self.data)

    def __getitem__(self, key):
        assert isinstance(key, slice)
        assert 0 < key.stop - key.start <= ranges._RANGE_BLOCK_ROWS
        self.reads += 1
        return self.data[key]


@pytest.mark.parametrize("block_rows", [1, 3, 8192])
def test_cached_ranges_match_scalar_and_membership_oracle(
    tmp_path, monkeypatch, block_rows
):
    """Cover gaps, absent IDs, inclusive endpoints, empty tables and block edges."""
    monkeypatch.setattr(ranges, "_RANGE_BLOCK_ROWS", block_rows)
    random = np.random.default_rng(2521)
    records = []
    for mid in (1, 3, 7, 9):
        active = random.random(71) > 0.5
        starts = np.flatnonzero(active & ~np.r_[False, active[:-1]])
        ends = np.flatnonzero(active & ~np.r_[active[1:], False])
        records.extend((mid, int(start), int(end)) for start, end in zip(starts, ends))
    for rows in (records, []):
        data = np.asarray(rows, dtype=np.int64).reshape((-1, 3))
        with h5py.File(tmp_path / "ranges.h5", "w") as file:
            columns = [file.create_dataset(str(i), data=data[:, i]) for i in range(3)]
            lookup = ranges._MoleculeRangeLookup(*columns)
            for mid in range(11):
                for frame in range(-1, 73):
                    expected = any(m == mid and a <= frame <= b for m, a, b in rows)
                    assert ranges.molecule_present(mid, frame, *columns) == expected
                    assert lookup.contains(mid, frame) == expected


def test_repeated_queries_reuse_blocks():
    """Once a small table is loaded, thousands of checks perform no more reads."""
    columns = [
        _CountedColumn(np.asarray(values, dtype=np.int64))
        for values in ([1, 1, 2], [0, 4, 2], [2, 8, 5])
    ]
    lookup = ranges._MoleculeRangeLookup(*columns)
    for _ in range(1000):
        assert lookup.contains(1, 2)
        assert not lookup.contains(1, 3)
        assert lookup.contains(1, 4)
        assert lookup.contains(2, 5)
        assert not lookup.contains(3, 5)
    assert [column.reads for column in columns] == [1, 1, 1]


def test_shared_cache_evicts_before_reading_and_stays_bounded():
    """A table much larger than the budget stays below 1 MiB across all columns."""
    size = ranges._RANGE_BLOCK_ROWS * (ranges._RANGE_CACHE_BLOCKS + 5) + 1
    data = np.arange(size, dtype=np.int64)
    columns = [_CountedColumn(data) for _ in range(3)]
    lookup = ranges._MoleculeRangeLookup(*columns)
    cache = lookup.columns[0].cache
    assert all(column.cache is cache for column in lookup.columns)
    for offset in range(0, size, ranges._RANGE_BLOCK_ROWS):
        for column in lookup.columns:
            assert column[offset] == offset
            assert len(cache) <= ranges._RANGE_CACHE_BLOCKS
            assert sum(block.nbytes for block in cache.values()) <= 1024 * 1024
    reads = columns[0].reads
    assert lookup.columns[0][0] == 0
    assert columns[0].reads == reads + 1


def test_cache_uses_lru_and_is_scoped_to_one_table(monkeypatch):
    """Keep recently used pages and never reuse values from a different file."""
    monkeypatch.setattr(ranges, "_RANGE_BLOCK_ROWS", 2)
    monkeypatch.setattr(ranges, "_RANGE_CACHE_BLOCKS", 2)
    data = _CountedColumn(np.arange(8, dtype=np.int64))
    lookup = ranges._MoleculeRangeLookup(data, data, data)
    column = lookup.columns[0]
    assert [column[i] for i in (0, 2, 0, 4, 0)] == [0, 2, 0, 4, 0]
    assert data.reads == 3
    assert column[2] == 2
    assert data.reads == 4
    other = ranges._MoleculeRangeLookup(np.ones(8, dtype=np.int64), data, data)
    assert other.columns[0][0] == 1
    with pytest.raises(IndexError):
        getitem(column, -1)
    with pytest.raises(IndexError):
        getitem(column, 8)


def test_read_failure_does_not_poison_cache():
    """A failed bounded read propagates and the same block can be retried."""
    data = _CountedColumn(np.asarray([1], dtype=np.int64))
    lookup = ranges._MoleculeRangeLookup(data, data, data)
    # Inject an I/O failure without changing the column's indexing contract.
    with patch.object(
        _CountedColumn, "__getitem__", side_effect=OSError("read failed")
    ):
        with pytest.raises(OSError, match="read failed"):
            lookup.contains(1, 1)
    assert lookup.contains(1, 1)


def test_cached_molecule_bounds_do_not_read_id_column_again(monkeypatch):
    """Changing frames reuses only the immutable ID bounds, never a presence result."""
    lookup = ranges._MoleculeRangeLookup(
        np.array([1, 1, 2]), np.array([0, 4, 2]), np.array([2, 8, 5])
    )
    assert lookup.contains(1, 1)
    original = ranges._RangeColumn.__getitem__

    def without_id_reads(column, index):
        if column.column == 0:
            raise AssertionError("warm bounds queried molecule IDs again")
        return original(column, index)

    monkeypatch.setattr(ranges._RangeColumn, "__getitem__", without_id_reads)
    assert not lookup.contains(1, 3)
    assert lookup.contains(1, 4)
    assert lookup.contains(1, 8)
    assert not lookup.contains(1, 9)


def test_molecule_bounds_are_lru_bounded_and_table_local(monkeypatch):
    """Evict ID bounds by recency without sharing them between range tables."""
    monkeypatch.setattr(ranges, "_MOLECULE_BOUND_CACHE_ENTRIES", 2)
    lookup = ranges._MoleculeRangeLookup(
        np.array([1, 2, 3]), np.array([0, 0, 0]), np.array([0, 0, 0])
    )
    for mid in (1, 2, 1, 3, 1):
        assert lookup.contains(mid, 0)
        assert len(lookup._bounds) <= 2
    assert list(lookup._bounds) == [3, 1]
    other = ranges._MoleculeRangeLookup(np.array([2]), np.array([0]), np.array([0]))
    assert not other.contains(1, 0)
    assert lookup.contains(1, 0)


def test_failed_presence_read_does_not_cache_a_result(monkeypatch):
    """Read failures remain retryable and cannot turn into cached presence."""
    columns = [_CountedColumn(np.array([v])) for v in (1, 2, 4)]
    lookup = ranges._MoleculeRangeLookup(*columns)
    assert lookup.contains(1, 2)
    lookup.columns[0].cache.clear()
    original = _CountedColumn.__getitem__

    def fail_start(column, index):
        if column is columns[1]:
            raise OSError("start read failed")
        return original(column, index)

    with monkeypatch.context() as context:
        context.setattr(_CountedColumn, "__getitem__", fail_start)
        with pytest.raises(OSError, match="start read failed"):
            lookup.contains(1, 3)
    assert lookup.contains(1, 3)
    assert not lookup.contains(1, 5)


def test_reader_can_request_smaller_blocks():
    """Honor a validator's smaller row budget for every cached read."""

    class TinyColumn(_CountedColumn):
        def __getitem__(self, key):
            assert key.stop - key.start <= 2
            return super().__getitem__(key)

    columns = [TinyColumn(np.arange(16) + delta) for delta in (1, 0, 0)]
    lookup = ranges._MoleculeRangeLookup(*columns, block_rows=2)
    assert lookup.contains(9, 8)
    assert not lookup.contains(9, 9)


@pytest.mark.parametrize("block_rows", [0, -1])
def test_invalid_reader_block_size_is_rejected(block_rows):
    """Reject invalid block sizes before accessing any range column."""
    with pytest.raises(ValueError, match="positive"):
        ranges._MoleculeRangeLookup([], [], [], block_rows=block_rows)
