# SPDX-License-Identifier: LGPL-3.0-or-later
"""Check indexed Step 3 execution against legacy output and bounded reads."""

import multiprocessing
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from reacnetgenerator import utils
from reacnetgenerator._path import (
    _CollectSMILESPaths,
    _initialize_route_worker,
    _route_changes,
)
from reacnetgenerator._reaction import ReactionsFinder
from reacnetgenerator._step3state import (
    _AtomFrameReader,
    _AtomFrameStore,
    _MoleculeNameTable,
)


class _BoundedTimeline:
    """Reject any attempt to copy a full timeline instead of a small slice."""

    def __init__(self, values, block_rows):
        self.values = values
        self.block_rows = block_rows

    def __len__(self):
        return len(self.values)

    def __getitem__(self, index):
        assert isinstance(index, slice)
        start, stop, step = index.indices(len(self))
        assert step == 1
        assert stop - start <= self.block_rows
        return self.values[index]


@pytest.mark.parametrize("block_rows", [1, 2, 7, 32])
@pytest.mark.parametrize(
    "values",
    [[], [0], [0] * 19, [3] * 19, [0, 1, 1, 0, 0, 2, 2, 0, 2, 1, 0]],
)
def test_route_scan_preserves_filtered_positions_and_raw_changes(values, block_rows):
    """Block boundaries and gaps must not change legacy route coordinates."""
    values = np.asarray(values, dtype=np.uint64)
    filtered = values[values != 0]
    expected_times = (
        np.r_[0, np.flatnonzero(np.diff(filtered)) + 1]
        if len(filtered)
        else np.zeros(0, dtype=int)
    )
    active = np.zeros(max(0, len(values) - 1), dtype=np.uint8)
    times, route = _route_changes(
        _BoundedTimeline(values, block_rows), active, block_rows=block_rows
    )
    np.testing.assert_array_equal(times, expected_times)
    np.testing.assert_array_equal(route, filtered[expected_times])
    np.testing.assert_array_equal(active, values[:-1] != values[1:])


def test_route_scan_handles_decreasing_wide_unsigned_ids():
    """Unsigned subtraction must not change transition detection or IDs."""
    values = np.array([2**63 + 1, 2**63 + 1, 1, 0, 1, 2**63 + 1], dtype=np.uint64)
    times, route = _route_changes(values, block_rows=2)
    np.testing.assert_array_equal(times, [0, 2, 4])
    np.testing.assert_array_equal(route, values[[0, 2, 5]])


def _use_context(monkeypatch, method):
    context = multiprocessing.get_context(method)
    for name in ("Pool", "Event", "Semaphore", "SimpleQueue"):
        monkeypatch.setattr(utils, name, getattr(context, name))


def _collector(tmp_path, names, nproc, shape):
    collector = object.__new__(_CollectSMILESPaths)
    collector.N, collector.step = shape
    collector.nproc = nproc
    collector.atomname = np.array(["H", "O"])
    collector.atomtype = np.array([0, 0, 1, 1])
    collector.selectatoms = ["H"]
    collector.mname = names
    collector.runHMM = True
    collector.atomroutefilename = str(tmp_path / "out.route")
    return collector


@pytest.mark.parametrize("method,nproc", [("spawn", 1), ("spawn", 2), ("fork", 2)])
def test_indexed_workers_preserve_routes_events_and_shared_file_ownership(
    tmp_path, monkeypatch, method, nproc
):
    """Exercise real pools, skipped frames, conflicts, split routes and cleanup."""
    if method not in multiprocessing.get_all_start_methods():
        pytest.skip(f"{method} is unavailable")
    _use_context(monkeypatch, method)
    names = _MoleculeNameTable.from_names(["A", "B", "C", "D", "E", "X", "Y"])
    matrix = np.array(
        [
            [1, 1, 3, 3, 3, 1, 1, 3, 3, 3, 3],
            [2, 2, 3, 3, 3, 2, 2, 3, 3, 3, 3],
            [4, 4, 0, 0, 4, 4, 4, 4, 5, 5, 5],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7],
        ]
    )
    collector = _collector(tmp_path, names, nproc, matrix.shape)
    seen_tasks = {"route": [], "reaction": []}
    original_run_mp = utils.run_mp

    def observe(kind):
        def run(nproc, **kwargs):
            inputs = kwargs.pop("l")

            def indices():
                for value in inputs:
                    assert type(value) is int
                    seen_tasks[kind].append(value)
                    yield value

            yield from original_run_mp(nproc, l=indices(), **kwargs)

        return run

    monkeypatch.setattr("reacnetgenerator._path.run_mp", observe("route"))
    monkeypatch.setattr("reacnetgenerator._reaction.run_mp", observe("reaction"))
    with _AtomFrameStore(matrix.shape, 7, directory=tmp_path) as store:
        store.atomeach[:] = matrix
        store.conflict.mark([3], [9], [[True]])
        store.prepare_active_transitions()
        collector._printatomroute(store)
        # All atoms contribute activity, even when only H routes are selected.
        expected_active = np.any(matrix[:, :-1] != matrix[:, 1:], axis=0)
        np.testing.assert_array_equal(store.active_transitions, expected_active)
        assert list(store.iter_active_transitions(block_rows=2)) == [1, 3, 4, 6, 7, 9]
        assert store.active_transition_count == 6
        assert Path(collector.atomroutefilename).read_text() == (
            "Atom 1 H: 0 A -> 2 C -> 5 A -> 7 C\n"
            "Atom 2 H: 0 B -> 2 C -> 5 B -> 7 C\n"
            "Atom 3 O: 0 D -> 6 E\n"
            "Atom 4 O: 0 X -> 10 Y\n"
        )
        collector._printatomroute(store, timeaxis=0, frame_range=(2, 8))
        assert Path(collector.atomroutefilename + ".0").read_text().splitlines()[0] == (
            "Atom 1 H: 0 C -> 3 A -> 5 C"
        )
        np.testing.assert_array_equal(store.active_transitions, expected_active)

        for events in (False, True):
            finder = ReactionsFinder(
                SimpleNamespace(
                    step=matrix.shape[1],
                    mname=names,
                    nproc=nproc,
                    printreactionevent=events,
                    reactionabcdfilename=str(tmp_path / "out.reactionabcd"),
                    reactioneventfilename=str(tmp_path / "out.events.csv"),
                )
            )
            finder.findreactions(store.atomeach.T, store.conflict.T, matrix_store=store)
            counts = Path(finder.reactionabcdfilename).read_text()
            expected_counts = "2 A+B->C\n1 C->A+B\n1 D->E\n"
            if events:
                assert counts == expected_counts
            else:
                # Legacy count-only mode also orders tied counts by worker
                # completion, so the scientific invariant is the line multiset.
                assert sorted(counts.splitlines()) == sorted(
                    expected_counts.splitlines()
                )
        assert Path(finder.reactioneventfilename).read_text() == (
            "Timestep_Index,Reactant,Product\n1,A+B,C\n4,C,A+B\n6,A+B,C\n7,D,E\n"
        )
        # Worker finalizers must not unlink mappings still owned by this process.
        assert all(
            Path(path).exists()
            for path in (
                store.atomeach_path,
                store.conflict.path,
                store.active_transition_path,
            )
        )
    assert seen_tasks["route"] == list(range(4)) * 2
    assert seen_tasks["reaction"] == [1, 3, 4, 6, 7, 9] * 2
    assert not list(tmp_path.glob("reacnetgenerator-*.mmap"))


def test_reader_is_readonly_and_does_not_remove_parent_files(tmp_path):
    """A worker may close its views independently of the owning store."""
    with _AtomFrameStore((2, 9), 5, directory=tmp_path) as store:
        store.atomeach[:] = 5
        store.conflict.mark([1], [8], [[True]])
        reader = _AtomFrameReader(*store.reader_args)
        try:
            assert reader.atomeach[0, 0] == 5
            assert reader.conflict[1, 8]
            with pytest.raises(ValueError, match="read-only"):
                reader.atomeach[0, 0] = 3
        finally:
            reader.close()
        assert Path(store.atomeach_path).exists()
        assert Path(store.conflict.path).exists()
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("frames", [1, 17])
def test_inactive_trajectory_and_single_frame_have_no_reaction_tasks(tmp_path, frames):
    """An empty active index must retain the event header and empty counts."""
    with _AtomFrameStore((1, frames), 1, directory=tmp_path) as store:
        store.atomeach[:] = 1
        store.prepare_active_transitions()
        assert list(store.iter_active_transitions(block_rows=2)) == []
        assert store.active_transition_count == 0
        finder = ReactionsFinder(
            SimpleNamespace(
                step=frames,
                mname=np.array(["A"]),
                nproc=1,
                printreactionevent=True,
                reactionabcdfilename=str(tmp_path / "counts"),
                reactioneventfilename=str(tmp_path / "events"),
            )
        )
        finder.findreactions(store.atomeach.T, store.conflict.T, matrix_store=store)
        assert Path(finder.reactionabcdfilename).read_bytes() == b""
        assert (
            Path(finder.reactioneventfilename).read_text()
            == "Timestep_Index,Reactant,Product\n"
        )
    assert not list(tmp_path.glob("*.mmap"))


def test_atom_blocks_join_one_reaction_and_preserve_conflict_suppression():
    """A connected reaction spanning blocks must not become multiple events."""
    finder = object.__new__(ReactionsFinder)
    finder.mname = np.array(["A", "B", "C"])
    finder.printreactionevent = False
    blocks = [
        (np.array([1, 1]), np.array([3, 3]), np.array([0, 0]), np.array([0, 0])),
        (np.array([2, 2]), np.array([3, 3]), np.array([0, 0]), np.array([0, 0])),
    ]
    assert finder._reactions_from_blocks(blocks, None) == ["A+B->C"]
    blocks[1][2][0] = 1
    assert finder._reactions_from_blocks(blocks, None) == []


def test_route_write_failure_stops_workers_before_store_cleanup(tmp_path, monkeypatch):
    """A parent-side output failure must not leave workers using deleted files."""
    names = _MoleculeNameTable.from_names(["A", "B"])
    collector = _collector(tmp_path, names, 2, (4, 5))
    children_before = {child.pid for child in multiprocessing.active_children()}

    def fail_write(self, value):
        raise OSError("simulated route write failure")

    with pytest.raises(OSError, match="simulated route write failure"):
        with _AtomFrameStore((4, 5), 2, directory=tmp_path) as store:
            store.atomeach[:] = [1, 1, 2, 2, 1]
            store.prepare_active_transitions()
            monkeypatch.setattr("reacnetgenerator._path.WriteBuffer.append", fail_write)
            collector._printatomroute(store)
    assert {child.pid for child in multiprocessing.active_children()} <= children_before
    assert not list(tmp_path.glob("reacnetgenerator-*.mmap"))


def test_active_allocation_failure_removes_only_partial_file(tmp_path, monkeypatch):
    """A failed new mapping must leave no orphan and keep existing state usable."""
    original_memmap = np.memmap
    with _AtomFrameStore((2, 9), 3, directory=tmp_path) as store:

        def fail_active(path, *args, **kwargs):
            if "reacnetgenerator-active-" in str(path):
                raise OSError("simulated active mapping failure")
            return original_memmap(path, *args, **kwargs)

        monkeypatch.setattr("reacnetgenerator._step3state.np.memmap", fail_active)
        with pytest.raises(OSError, match="simulated active"):
            store.prepare_active_transitions()
        assert not list(tmp_path.glob("reacnetgenerator-active-*"))
        assert store.atomeach.shape == (2, 9)
    assert not list(tmp_path.iterdir())


_INITIALIZED_VALUE = None
_INITIALIZATION_COUNT = 0


def _initialize_value(value):
    global _INITIALIZED_VALUE, _INITIALIZATION_COUNT
    _INITIALIZED_VALUE = value
    _INITIALIZATION_COUNT += 1


def _read_initialized_value(value):
    return _INITIALIZED_VALUE, _INITIALIZATION_COUNT, value


@pytest.mark.parametrize("mode", ["bounded", "ordered"])
def test_run_mp_initializer_runs_once_before_tasks(monkeypatch, mode):
    """Bounded paths should pass initializer arguments once per worker."""
    _use_context(monkeypatch, "spawn")
    options = {}
    if mode != "legacy":
        options = {"max_inflight": 2, "chunksize": 1, "disk_ordered": mode == "ordered"}
    results = list(
        utils.run_mp(
            1,
            func=_read_initialized_value,
            l=[3, 4],
            initializer=_initialize_value,
            initargs=("attached",),
            unordered=mode != "ordered",
            total=2,
            bar=False,
            **options,
        )
    )
    assert results == [("attached", 1, 3), ("attached", 1, 4)]


def test_run_mp_rejects_initializer_without_startup_monitoring():
    """The legacy Pool path cannot safely propagate initialization failures."""
    with pytest.raises(ValueError, match="initializer requires"):
        list(
            utils.run_mp(
                1,
                func=_read_initialized_value,
                l=[0],
                initializer=_initialize_value,
                initargs=("attached",),
            )
        )


def _missing_mapping_case(queue, directory):
    try:
        list(
            utils.run_mp(
                1,
                func=_read_initialized_value,
                l=[0],
                initializer=_initialize_route_worker,
                initargs=(
                    (str(Path(directory) / "missing"), "missing", (1, 2), "u1"),
                    [0],
                    ["H"],
                    ["H"],
                    ["A"],
                    (0, 2),
                    None,
                ),
                max_inflight=1,
                chunksize=1,
                bar=False,
            )
        )
    except RuntimeError as error:
        queue.put(str(error))


def test_mapping_initializer_failure_does_not_hang(tmp_path):
    """Failure before worker readiness must terminate the bounded pool."""
    context = multiprocessing.get_context("spawn")
    queue = context.Queue()
    process = context.Process(target=_missing_mapping_case, args=(queue, str(tmp_path)))
    process.start()
    process.join(timeout=30)
    try:
        assert not process.is_alive(), (
            "run_mp hung after a mapping initialization failure"
        )
        assert process.exitcode == 0
        assert "worker" in queue.get(timeout=2).lower()
    finally:
        if process.is_alive():
            process.terminate()
            process.join(timeout=5)
        queue.close()
