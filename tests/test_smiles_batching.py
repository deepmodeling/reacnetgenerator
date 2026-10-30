# SPDX-License-Identifier: LGPL-3.0-or-later
"""SMILES admission, ordering, resource bounds and real worker failure paths."""

import multiprocessing
import os
import pickle
import time
from contextlib import closing
from pathlib import Path

import numpy as np
import pytest

from reacnetgenerator import _path
from reacnetgenerator import _smilesworker as workers
from reacnetgenerator.utils import listtobytes, run_mp


def _collector(tmp_path, records, miso=0, nproc=2):
    """Build actual compressed definitions and the production SMILES collector."""
    tmp_path.mkdir(exist_ok=True)
    source = tmp_path / "records.bin"
    with source.open("wb") as file:
        for atoms, bonds, frames in records:
            for value in (
                np.asarray(atoms),
                np.asarray([(a, b) for a, b, _ in bonds]),
                np.asarray([level for _, _, level in bonds]),
                np.asarray(frames),
            ):
                file.write(listtobytes(value))
    count = max((max(atoms) for atoms, _, _ in records), default=0) + 1
    collector = object.__new__(_path._CollectSMILESPaths)
    collector.hmmit = len(records)
    collector.nproc, collector.miso = nproc, miso
    collector.atomname = np.asarray(["C"])
    collector.atomtype = np.zeros(count, dtype=int)
    collector.atomnames = np.full(count, "C")
    collector.moleculetemp2filename = str(source)
    collector.moleculefilename = str(tmp_path / "names.txt")
    collector.moleculetimelinefilename = str(tmp_path / "timeline.csv")
    collector.printmoleculetime = True
    collector._moleculeframefilter = None
    collector._moleculetimestepfilter = None
    collector.timestep = {int(f): int(f) * 10 for _, _, fs in records for f in fs}
    return collector


def _records(count=130):
    """Include full structures differing in bond order and equivalent ring forms."""
    result = []
    for index in range(count):
        if index % 3 == 2:
            bonds = [(a, (a + 1) % 6, 1 + (a + index) % 2) for a in range(6)]
            result.append((list(range(6)), bonds, [index, index + 1]))
        else:
            result.append(([0, 1], [(0, 1, 1 + index % 3)], [index]))
    return result


@pytest.fixture
def force_batches(monkeypatch, tmp_path):
    """Force small regression fixtures through the actual initialized pool."""
    monkeypatch.setattr(workers, "_SMILES_SPAWN_RECORDS_PER_WORKER", 1)
    monkeypatch.setattr(workers, "_SMILES_FORK_RECORDS_PER_WORKER", 1)
    monkeypatch.setattr(workers.tempfile, "tempdir", str(tmp_path))
    return tmp_path


def _assert_closed(folder, children):
    assert not list(folder.glob("reacnetgenerator-smiles-*"))
    assert {p.pid for p in multiprocessing.active_children()} == children


@pytest.mark.parametrize("miso", [0, 1, 2])
@pytest.mark.parametrize("spill", [False, True])
def test_real_batched_names_match_local_outputs(
    tmp_path, monkeypatch, force_batches, miso, spill
):
    """Both result transports retain canonical names, representatives and CSV order."""
    records = _records()
    local = _collector(tmp_path / "local", records, miso, 1)
    parallel = _collector(tmp_path / "parallel", records, miso)
    local._printmoleculename()
    if spill:
        monkeypatch.setattr(workers, "_SMILES_RESULT_BYTES", 128)
    calls = []

    def checked(nproc, **kwargs):
        assert kwargs["initializer"] is workers._initialize_smiles_worker
        assert kwargs["disk_ordered"] and not kwargs["unordered"]
        assert kwargs["max_inflight"] == 2 * nproc
        calls.append(kwargs["total"])
        return run_mp(nproc, **kwargs)

    monkeypatch.setattr(workers, "run_mp", checked)
    children = {p.pid for p in multiprocessing.active_children()}
    parallel._printmoleculename()
    assert calls and 1 < calls[0] < len(records)
    assert list(parallel.mname) == list(local.mname)
    assert (
        Path(parallel.moleculefilename).read_bytes()
        == Path(local.moleculefilename).read_bytes()
    )
    assert (
        Path(parallel.moleculetimelinefilename).read_bytes()
        == Path(local.moleculetimelinefilename).read_bytes()
    )
    _assert_closed(tmp_path, children)


def test_small_spawn_retains_original_ordered_runner(tmp_path, monkeypatch):
    """The fallback is exercised at the collector entry, without a batch pre-scan."""
    collector = _collector(tmp_path, _records(4))
    monkeypatch.setattr(workers, "get_start_method", lambda: "spawn")
    calls = []

    def original(nproc, **kwargs):
        assert kwargs["func"].__self__ is collector
        assert kwargs["nlines"] == 4 and kwargs["total"] == 4
        assert "initializer" not in kwargs
        calls.append(nproc)
        return run_mp(nproc, **kwargs)

    def unexpected(*args, **kwargs):
        raise AssertionError("Small spawn inputs must not plan batches")

    monkeypatch.setattr(_path, "run_mp", original)
    monkeypatch.setattr(workers, "_smiles_batches", unexpected)
    collector._printmoleculename()
    assert calls == [2]


def test_admission_is_bounded_by_decoded_and_wire_sizes(tmp_path, monkeypatch):
    """A tiny compressed frame range must not share a batch with other records."""
    collector = _collector(
        tmp_path, [*_records(2), ([0], [], np.arange(100000)), *_records(2)]
    )
    with open(collector.moleculetemp2filename, "rb") as file:
        batches = list(workers._smiles_batches(file))
    assert list(map(len, batches)) == [2, 1, 2]
    # Force the wire budget to split records that otherwise share a decoded budget.
    monkeypatch.setattr(workers, "_SMILES_INPUT_BYTES", 5500)
    with open(collector.moleculetemp2filename, "rb") as file:
        tasks = list(workers._encoded_smiles_batches(file))
    for _, payload in tasks:
        batch = pickle.loads(payload)
        assert len(batch) == 1 or len(payload) <= 5500
        if len(batch) > 1:
            costs = [workers._smiles_record_cost(record) for record in batch]
            assert (
                sum(cost[1] for cost in costs) <= workers._SMILES_DECODED_ESTIMATE_BYTES
            )
    assert sum(len(pickle.loads(payload)) for _, payload in tasks) == 5


def test_result_writer_spills_before_retaining_oversized_write(tmp_path):
    """Disk fallback accepts an indivisible large write without growing the prefix."""
    path = tmp_path / "large.pickle"
    with closing(workers._SMILESResultWriter(path, 128)) as writer:
        writer.write(b"a" * 64)
        writer.write(b"b" * 10000)
        assert writer.memory.closed and writer.result() is None
    assert path.read_bytes() == b"a" * 64 + b"b" * 10000


def test_worker_streams_large_decoded_results(tmp_path, monkeypatch):
    """Large frame arrays are single-record tasks with bounded result envelopes."""
    collector = _collector(tmp_path, [([0], [], np.arange(100000))])
    monkeypatch.setattr(workers, "_SMILES_WORKER_STATE", None)
    workers._initialize_smiles_worker(
        collector.atomname, collector.atomnames, True, False, str(tmp_path), 1024
    )
    with open(collector.moleculetemp2filename, "rb") as file:
        task = next(workers._encoded_smiles_batches(file))
    index, count, payload = workers._convert_smiles_batch(task)
    assert count == 1 and payload is None
    with (tmp_path / f"{index}.pickle").open("rb") as file:
        name, atoms, bonds, frames = pickle.load(file)
        assert file.read() == b""
    assert name == "[C]" and bonds == []
    np.testing.assert_array_equal(atoms, [0])
    np.testing.assert_array_equal(frames, np.arange(100000))


def _delayed_batch(task, folder):
    """Let a later batch finish before the first while recording completions."""
    index, _ = task
    if index == 0:
        time.sleep(0.5)
    result = workers._convert_smiles_batch(task)
    with open(Path(folder) / "completed", "a") as file:
        file.write(f"{index}\n")
    return result


def test_ordered_window_and_early_close(tmp_path, monkeypatch, force_batches):
    """A delayed first result cannot trigger unbounded submission or orphan files."""
    from functools import partial

    monkeypatch.setattr(workers, "_SMILES_BATCH_RECORDS", 1)
    monkeypatch.setattr(workers, "_SMILES_RESULT_BYTES", 128)
    collector = _collector(tmp_path, _records(12))
    submitted = []
    windows = []

    def instrument(nproc, **kwargs):
        source = kwargs["l"]
        windows.append(kwargs["max_inflight"])

        def inputs():
            for item in source:
                submitted.append(item[0])
                yield item

        kwargs["l"] = inputs()
        kwargs["func"] = partial(_delayed_batch, folder=str(tmp_path))
        return run_mp(nproc, **kwargs)

    monkeypatch.setattr(workers, "run_mp", instrument)
    children = {p.pid for p in multiprocessing.active_children()}
    with open(collector.moleculetemp2filename, "rb") as file:
        with closing(
            collector._getSMILESresults(file, collector._calmoleculeSMILESname)
        ) as results:
            first = next(results)
            assert list(first[-1]) == [0]
            assert len(submitted) <= windows[0]
    completions = (tmp_path / "completed").read_text().splitlines()
    assert completions.index("1") < completions.index("0")
    _assert_closed(tmp_path, children)


def _failed_initializer(*args):
    raise RuntimeError("intentional initializer failure")


def _failed_task(task):
    raise RuntimeError("intentional conversion failure")


def _exited_task(task):
    os._exit(9)


def _unpicklable_result(record):
    return lambda: None


def _serialization_initializer(*args):
    workers._initialize_smiles_worker(*args)
    _, folder, limit = workers._SMILES_WORKER_STATE
    workers._SMILES_WORKER_STATE = _unpicklable_result, folder, limit


@pytest.mark.parametrize(
    "failure", ["initializer", "task", "exit", "serialize", "consumer"]
)
def test_batch_failures_close_workers_and_streams(
    tmp_path, monkeypatch, force_batches, failure
):
    """All failures cross the collector/runner boundary without leaking resources."""
    monkeypatch.setattr(workers, "_SMILES_RESULT_BYTES", 128)
    collector = _collector(tmp_path, _records())

    def inject(nproc, **kwargs):
        if failure == "initializer":
            kwargs["initializer"] = _failed_initializer
        elif failure == "serialize":
            kwargs["initializer"] = _serialization_initializer
        elif failure in ("task", "exit"):
            kwargs["func"] = _failed_task if failure == "task" else _exited_task
        return run_mp(nproc, **kwargs)

    monkeypatch.setattr(workers, "run_mp", inject)
    if failure == "consumer":

        def fail_consumer(*args, **kwargs):
            raise RuntimeError("intentional consumer failure")

        monkeypatch.setattr(
            _path._CollectSMILESPaths, "_formatmoleculename", fail_consumer
        )
    children = {p.pid for p in multiprocessing.active_children()}
    with pytest.raises(
        (RuntimeError, AttributeError), match=r"intentional|exited|pickl|local object"
    ):
        collector._printmoleculename()
    assert Path(collector.moleculetemp2filename).exists()
    _assert_closed(tmp_path, children)


@pytest.mark.parametrize("delta", [-1, 1])
def test_batch_preflight_rejects_record_count_mismatch(
    tmp_path, monkeypatch, force_batches, delta
):
    """Pre-counting batches preserves the runner's exact-count input contract."""
    collector = _collector(tmp_path, _records(5))
    collector.hmmit += delta

    def unexpected(*args, **kwargs):
        raise AssertionError("A count mismatch must fail before starting workers")

    monkeypatch.setattr(workers, "run_mp", unexpected)
    with pytest.raises(ValueError, match="record count"):
        collector._printmoleculename()


def _selected_failure(self, atoms, bonds):
    """Fail two distinct structures without hiding the successful candidate."""
    if bonds[0][2] in (1, 2):
        raise ValueError("forced SMILES failure")
    return _path._CollectPaths.convertSMILES(self, atoms, bonds)


def _fallback_initializer(*args):
    from types import MethodType

    workers._initialize_smiles_worker(*args)
    function, _, _ = workers._SMILES_WORKER_STATE
    function.__self__.convertSMILES = MethodType(_selected_failure, function.__self__)


@pytest.mark.parametrize("miso", [0, 1, 2])
def test_failed_candidates_keep_distinct_global_frequencies(
    tmp_path, monkeypatch, force_batches, miso
):
    """Batched failed names preserve VF2 identities before representative selection."""
    from types import MethodType

    records = [
        ([0, 1], [(0, 1, order)], list(range(count)))
        for order, count in ((1, 4), (2, 4), (3, 7))
    ]
    local = _collector(tmp_path / "local", records, miso, 1)
    local.convertSMILES = MethodType(_selected_failure, local)
    local._printmoleculename()
    parallel = _collector(tmp_path / "parallel", records, miso)

    def inject(nproc, **kwargs):
        kwargs["initializer"] = _fallback_initializer
        return run_mp(nproc, **kwargs)

    monkeypatch.setattr(workers, "run_mp", inject)
    parallel._printmoleculename()
    assert list(parallel.mname) == list(local.mname)
    if miso:
        assert list(parallel.mname) == ["[C]#[C]"] * 3
    else:
        assert list(parallel.mname) == [
            "C2_unknownSMILES_0",
            "C2_unknownSMILES_1",
            "[C]#[C]",
        ]


def _state_probe(record, *, function):
    return (
        os.getpid(),
        id(function.__self__),
        tuple(vars(function.__self__)),
        function(record),
    )


def _probe_initializer(*args):
    from functools import partial

    workers._initialize_smiles_worker(*args)
    function, folder, limit = workers._SMILES_WORKER_STATE
    workers._SMILES_WORKER_STATE = (
        partial(_state_probe, function=function),
        folder,
        limit,
    )


def test_converter_reused_without_rng_metadata(tmp_path, monkeypatch, force_batches):
    """Multiple batches in each real worker share only the needed converter state."""
    monkeypatch.setattr(workers, "_SMILES_BATCH_RECORDS", 1)
    collector = _collector(tmp_path, _records(20))

    def instrument(nproc, **kwargs):
        kwargs["initializer"] = _probe_initializer
        return run_mp(nproc, **kwargs)

    monkeypatch.setattr(workers, "run_mp", instrument)
    with open(collector.moleculetemp2filename, "rb") as file:
        results = list(workers._run_smiles_batches(collector, file, False))
    states = {}
    for pid, identity, keys, _ in results:
        states.setdefault(pid, set()).add(identity)
        assert set(keys) == {
            "atomname",
            "atomnames",
            "printmoleculetime",
            "_moleculeframefilter",
            "_moleculetimestepfilter",
        }
    assert len(results) == 20
    assert all(len(identities) == 1 for identities in states.values())
