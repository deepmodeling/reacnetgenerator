# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compare bounded route reduction with the original event-array semantics."""

import multiprocessing
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from reacnetgenerator import utils
from reacnetgenerator._matrix import _GenerateMatrix
from reacnetgenerator._path import (
    _CollectSMILESPaths,
    _compact_route_pairs,
    _SpeciesRouteAccumulator,
)
from reacnetgenerator._step3state import _AtomFrameStore, _MoleculeNameTable
from reacnetgenerator.utils import WriteBuffer


def _legacy_counts(atom_pairs, names, run_hmm):
    """Retain the baseline event arrays as a small, independent test oracle."""
    rows, owners = [], {}
    for atom, pairs in enumerate(atom_pairs):
        for pair in pairs:
            key = tuple(pair)
            if run_hmm or owners.get(key, atom) == atom:
                owners[key] = atom
                rows.append(pair)
    if not rows:
        return Counter()
    rows = np.asarray(rows)
    if run_hmm:
        rows = np.unique(rows, axis=0)
    species = names[rows - 1]
    return Counter(map(tuple, species[species[:, 0] != species[:, 1]]))


@pytest.mark.parametrize("run_hmm", [False, True])
@pytest.mark.parametrize("compact_names", [False, True])
def test_counts_preserve_first_atom_multiplicity_and_species_folding(
    monkeypatch, run_hmm, compact_names
):
    """Repeated occurrences belong to the first atom, not every later atom."""
    monkeypatch.setattr("reacnetgenerator._path._STEP3_SCAN_ROWS", 3)
    monkeypatch.setattr("reacnetgenerator._path._ROUTE_COMPACTION_MIN_ROWS", 4)
    names = np.array(["A", "B", "A", "C"], dtype=object)
    if compact_names:
        names = _MoleculeNameTable.from_names(names)
    inputs = [
        np.array([[1, 2]] * 5 + [[2, 1]] * 3 + [[1, 3]] * 4),
        np.array([[1, 2]] * 7 + [[3, 2]] * 2 + [[2, 4]]),
        np.array([[2, 4]] * 4 + [[4, 1]]),
        np.zeros((0, 2), dtype=int),
    ]
    for compact in (False, True):
        reducer = _SpeciesRouteAccumulator(names, run_hmm)
        for atom, original in enumerate(inputs):
            pairs, counts = (
                _compact_route_pairs(original) if compact else (original, None)
            )
            reducer.add(atom, pairs, counts)
        assert reducer.counts == _legacy_counts(inputs, names, run_hmm)
        assert reducer.counts["A", "B"] == (2 if run_hmm else 7)
        assert reducer.counts["B", "C"] == 1
        assert ("A", "A") not in reducer.counts


def test_name_lookup_is_bounded(monkeypatch):
    """Unique routes must not create a full event-by-species name array."""
    monkeypatch.setattr("reacnetgenerator._path._STEP3_SCAN_ROWS", 7)

    class Names:
        def __len__(self):
            return 100

        def __getitem__(self, index):
            assert index.shape[1] == 2
            assert len(index) <= 7
            return np.asarray(["A", "B"], dtype=object)[index % 2]

    pairs = np.column_stack((np.arange(1, 100), np.arange(2, 101)))
    reducer = _SpeciesRouteAccumulator(Names(), False)
    reducer.add(0, pairs, None)
    assert sum(reducer.counts.values()) == len(pairs)


@pytest.mark.parametrize(
    "dtype,largest",
    [
        (np.uint8, 255),
        (np.uint16, 65535),
        (np.uint32, 2**32 - 1),
        (np.uint64, 2**63 + 1),
    ],
)
def test_compaction_preserves_integer_ids_and_counts_across_blocks(
    monkeypatch, dtype, largest
):
    """Packed narrow IDs and wide fallbacks preserve exact values and counts."""
    monkeypatch.setattr("reacnetgenerator._path._STEP3_SCAN_ROWS", 37)
    pairs = np.tile(np.array([[largest, 1], [1, largest]], dtype=dtype), (1025, 1))
    compact, counts = _compact_route_pairs(pairs)
    assert compact.dtype == pairs.dtype
    assert {tuple(pair): int(count) for pair, count in zip(compact, counts)} == Counter(
        map(tuple, pairs)
    )
    assert compact.nbytes + counts.nbytes < pairs.nbytes // 100


def test_unique_routes_and_compaction_cap_keep_original_array(monkeypatch):
    """A unique sample or an exceeded dictionary budget must use the raw path."""
    pairs = np.column_stack((np.arange(10000), np.arange(1, 10001)))
    result, counts = _compact_route_pairs(pairs)
    assert result is pairs and counts is None
    monkeypatch.setattr("reacnetgenerator._path._ROUTE_COMPACTION_MAX_PAIRS", 8)
    monkeypatch.setattr("reacnetgenerator._path._ROUTE_COMPACTION_SAMPLE_ROWS", 16)
    mixed = pairs.copy()
    mixed[:: max(1, len(mixed) // 16)] = [1, 2]
    result, counts = _compact_route_pairs(mixed)
    assert result is mixed and counts is None


@pytest.mark.parametrize("center", [None, "B"])
def test_matrix_accepts_species_counts_and_legacy_arrays(tmp_path, center):
    """Ranking, species limits, and text/table bytes remain identical."""
    matrix = object.__new__(_GenerateMatrix)
    matrix.mname = _MoleculeNameTable.from_names(["A", "B", "A", "C"])
    matrix.matrix_size = 2
    matrix.n_searchspecies = 1
    matrix.speciescenter = center
    rows = np.array([[1, 2], [1, 2], [3, 2], [2, 4], [4, 1], [1, 3]])
    counts = Counter({("A", "B"): 3, ("B", "C"): 1, ("C", "A"): 1})
    results = []
    for mode, data in enumerate((rows, counts)):
        matrix.reactionfilename = str(tmp_path / f"{mode}.reaction")
        matrix.tablefilename = str(tmp_path / f"{mode}.table")
        assert dict(matrix._getallroute(data)) == counts
        matrix._printtable(matrix._getallroute(data))
        results.append(
            (
                Path(matrix.reactionfilename).read_bytes(),
                Path(matrix.tablefilename).read_bytes(),
            )
        )
    assert results[0] == results[1]


def _expected_routes(matrix, atomtypes, selected, names):
    rows, lines = [], []
    for index, (timeline, atomtype) in enumerate(zip(matrix, atomtypes, strict=True)):
        filtered = timeline[timeline != 0]
        changes = (
            np.r_[0, np.flatnonzero(filtered[1:] != filtered[:-1]) + 1]
            if len(filtered)
            else []
        )
        route = filtered[changes]
        rows.append(
            np.column_stack((route[:-1], route[1:]))
            if atomtype in selected
            else np.zeros((0, 2), dtype=int)
        )
        lines.append(
            f"Atom {index + 1} {atomtype}: "
            + " -> ".join(
                f"{t} {name}" for t, name in zip(changes, names[route - 1], strict=True)
            )
        )
    return rows, "\n".join(lines) + "\n"


@pytest.mark.parametrize("method,nproc", [("spawn", 1), ("spawn", 2), ("fork", 2)])
@pytest.mark.parametrize("run_hmm", [False, True])
def test_real_workers_preserve_counts_route_text_split_and_activity(
    tmp_path, monkeypatch, method, nproc, run_hmm
):
    """Check aggregation through real workers while retaining all atom activity."""
    if method not in multiprocessing.get_all_start_methods():
        pytest.skip(f"{method} is unavailable")
    context = multiprocessing.get_context(method)
    for name in ("Pool", "Event", "Semaphore", "SimpleQueue"):
        monkeypatch.setattr(utils, name, getattr(context, name))
    monkeypatch.setattr("reacnetgenerator._path._STEP3_SERIAL_ATOM_FRAMES", 0)
    matrix = np.stack(
        [np.resize(np.repeat([1, 2], 2), 1025)] * 2
        + [np.r_[np.full(3, 3), np.full(1022, 4)]]
    )
    names = _MoleculeNameTable.from_names(["A", "B", "C", "D"])
    collector = object.__new__(_CollectSMILESPaths)
    collector.N, collector.step = matrix.shape
    collector.nproc, collector.runHMM = nproc, run_hmm
    collector.atomname, collector.atomtype = np.array(["H", "O"]), np.array([0, 0, 1])
    collector.selectatoms, collector.mname = ["H"], names
    collector.atomroutefilename = str(tmp_path / "out.route")
    with _AtomFrameStore(matrix.shape, 4, directory=tmp_path) as store:
        store.atomeach[:] = matrix
        store.prepare_active_transitions()
        for part, (start, stop) in enumerate(((0, 1025), (3, 770))):
            timeaxis = None if part == 0 else 0
            counts = collector._printatomroute(
                store, timeaxis=timeaxis, frame_range=(start, stop)
            )
            rows, text = _expected_routes(
                matrix[:, start:stop], ["H", "H", "O"], ["H"], names
            )
            assert counts == _legacy_counts(rows, names, run_hmm)
            filename = collector.atomroutefilename + ("" if part == 0 else ".0")
            assert Path(filename).read_text() == text
            np.testing.assert_array_equal(
                store.active_transitions,
                np.any(matrix[:, :-1] != matrix[:, 1:], axis=0),
            )
    assert not list(tmp_path.glob("reacnetgenerator-*.mmap"))


@pytest.mark.parametrize("binary", [False, True])
def test_write_buffer_byte_limit_preserves_unicode_and_large_records(binary):
    """Keep output bytes identical while flushing before the next large record."""
    records = ["氢", "氧", "small", "x" * 20, "末"]
    separator = "|"
    if binary:
        records = [value.encode() for value in records]
        separator = separator.encode()

    class Sink:
        name = "memory"

        def __init__(self):
            self.chunks = []

        def writelines(self, lines):
            self.chunks.append((b"" if binary else "").join(lines))

        def __exit__(self, *args):
            return None

    sink = Sink()
    with WriteBuffer(sink, sep=separator, maxbytes=10) as writer:
        writer.extend(iter(records))
        assert writer._buffer_bytes < 10
    empty = b"" if binary else ""
    assert empty.join(sink.chunks) == separator.join(records) + separator
    for chunk in sink.chunks:
        encoded = chunk if binary else chunk.encode()
        assert len(encoded) <= 10 or encoded == b"x" * 20 + b"|"
    with pytest.raises(ValueError, match="positive"):
        WriteBuffer(sink, sep=separator, maxbytes=0)
