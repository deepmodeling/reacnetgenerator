# SPDX-License-Identifier: LGPL-3.0-or-later
"""Bounded long-signal HMM scheduling must retain aligned molecule records."""

import multiprocessing
import os
import time
from functools import partial
from multiprocessing.reduction import ForkingPickler
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from reacnetgenerator import _hmmfilter, utils
from reacnetgenerator._moleculepresence import _iter_molecule_records
from tests.test_indexed_step3 import _use_context


def _use_hmm_context(monkeypatch, method):
    """Select both the actual pool and the HMM transport policy explicitly."""
    _use_context(monkeypatch, method)
    monkeypatch.setattr(_hmmfilter, "get_start_method", lambda: method)


def _filter(folder, nproc, origin=True, filtered=True, frame_count=128, count=12):
    """Build sparse, empty and dense signals with distinguishable definitions."""
    folder.mkdir()
    source = folder / "records.bin"
    with source.open("wb") as file:
        for index in range(count):
            frames = (
                np.zeros(0, dtype=int)
                if index % 3 == 0
                else np.arange(frame_count)
                if index % 3 == 1
                else np.arange(16, frame_count - 16)
            )
            for value in (
                np.array([index]),
                np.empty((0, 2), dtype=int),
                np.zeros(0, dtype=int),
                frames,
            ):
                file.write(utils.listtobytes(value))
    result = object.__new__(_hmmfilter._HMMFilter)
    result.nproc, result.step, result.temp1it = nproc, frame_count, count
    result.runHMM, result.getoriginfile, result.printfiltersignal = (
        True,
        origin,
        filtered,
    )
    result.moleculetempfilename = str(source)
    result.p = np.array([0.5, 0.5])
    result.a = np.array([[0.999, 0.001], [0.001, 0.999]])
    result.b = np.array([[0.99, 0.01], [0.01, 0.99]])
    result._initHMM()
    return result


def _read(filter):
    """Join the separately written streams by consumption order, then by atom ID."""
    result = {}
    with open(filter.moleculetemp2filename, "rb") as file:
        records = list(_iter_molecule_records(file))
    signals = []
    for filename in (filter.originfilename, filter.hmmfilename):
        if filename is None:
            signals.append([None] * len(records))
        else:
            with open(filename, "rb") as file:
                signals.append(
                    [
                        utils.bytestolist(block).tobytes()
                        for block in utils.read_compressed_block(file)
                    ]
                )
    for record, origin, filtered in zip(records, *signals, strict=True):
        key = int(utils.bytestolist(record[0])[0])
        assert key not in result
        result[key] = (tuple(record), origin, filtered)
    assert len(result) == filter.hmmit
    return result


@pytest.mark.parametrize("method", ["fork", "spawn", "forkserver"])
@pytest.mark.parametrize(
    "origin,filtered", [(False, False), (True, False), (False, True), (True, True)]
)
def test_long_hmm_keeps_aligned_outputs(
    tmp_path, monkeypatch, method, origin, filtered
):
    """Actual pools retain zero-signal filtering and both optional signal outputs."""
    if method not in multiprocessing.get_all_start_methods():
        pytest.skip(f"{method} unavailable")
    _use_hmm_context(monkeypatch, method)
    monkeypatch.setattr(_hmmfilter, "_HMM_BOUNDED_SIGNAL_FRAMES", 64)
    monkeypatch.setattr(_hmmfilter, "_HMM_LEGACY_CHUNK_RECORDS", 4)
    monkeypatch.setattr(_hmmfilter.tempfile, "tempdir", str(tmp_path))
    serial = _filter(tmp_path / "serial", 1, origin, filtered)
    parallel = _filter(tmp_path / "parallel", 2, origin, filtered)
    before = {p.pid for p in multiprocessing.active_children()}
    serial._calhmm()
    parallel._calhmm()
    assert _read(serial) == _read(parallel)
    assert parallel.hmmit == (8 if not origin and not filtered else 12)
    assert {p.pid for p in multiprocessing.active_children()} == before


@pytest.mark.parametrize("method", ["spawn", "forkserver"])
@pytest.mark.parametrize("frames,bounded", [(99999, False), (100000, True)])
def test_long_hmm_dispatch_boundary(tmp_path, monkeypatch, frames, bounded, method):
    """Select bounded long-signal tasks while retaining forkserver chunks."""
    monkeypatch.setattr(_hmmfilter, "get_start_method", lambda: method)
    filter = _filter(tmp_path / "input", 4, frame_count=frames, count=1)
    filter.temp1it = 201
    calls = []

    def capture(nproc, **kwargs):
        """Inspect dispatch parameters without running long boundary signals."""
        calls.append((nproc, kwargs))
        return (item for item in ())

    monkeypatch.setattr(_hmmfilter, "run_mp", capture)
    monkeypatch.setattr(_hmmfilter.tempfile, "tempdir", str(tmp_path))
    filter._calhmm()
    assert calls[0][0] == (3 if bounded else 4)
    options = calls[0][1]
    per_record = bounded and method != "forkserver"
    assert options.get("max_inflight") == (6 if per_record else None)
    assert options.get("chunksize") == (1 if per_record else None)
    if per_record:
        assert options["func"] is _hmmfilter._filter_hmm_record
        assert options["initializer"] is _hmmfilter._initialize_hmm_worker
        assert options["initargs"] == (
            filter.step,
            filter._model,
            filter.getoriginfile,
            filter.printfiltersignal,
        )
    else:
        assert "initializer" not in options
        if bounded:
            assert options["func"].func is _hmmfilter._convert_hmm_record
            assert options["func"].keywords["model"] is filter._model
            assert len(ForkingPickler.dumps(options["func"])) < 8192


def _delayed_record(record, *, function):
    """Make the first definition finish after its neighbors."""
    if int(utils.bytestolist(record[0])[0]) == 0:
        time.sleep(0.2)
    return function(record)


@pytest.mark.parametrize("bounded", [False, True])
def test_out_of_order_hmm_keeps_streams_together(tmp_path, monkeypatch, bounded):
    """Delayed workers must preserve canonical order and aligned signals."""
    _use_hmm_context(monkeypatch, "spawn")
    monkeypatch.setattr(
        _hmmfilter, "_HMM_BOUNDED_SIGNAL_FRAMES", 1 if bounded else 100000
    )
    monkeypatch.setattr(_hmmfilter, "_HMM_LEGACY_CHUNK_RECORDS", 4)
    monkeypatch.setattr(_hmmfilter.tempfile, "tempdir", str(tmp_path))
    count = 12 if bounded else 201
    serial = _filter(tmp_path / "serial", 1, count=count)
    parallel = _filter(tmp_path / "parallel", 2, count=count)
    serial._calhmm()

    def delayed(nproc, **kwargs):
        """Delay an input inside the real worker pool."""
        kwargs["func"] = partial(_delayed_record, function=kwargs["func"])
        return utils.run_mp(nproc, **kwargs)

    monkeypatch.setattr(_hmmfilter, "run_mp", delayed)
    parallel._calhmm()
    assert _read(parallel) == _read(serial)
    assert list(_read(parallel)) == list(range(count))


def _failed_record(record):
    raise RuntimeError("intentional HMM worker failure")


def _exited_record(record):
    os._exit(9)


@pytest.mark.parametrize("failure", ["worker", "exit", "consumer"])
def test_long_hmm_failure_stops_workers(tmp_path, monkeypatch, failure):
    """Worker failure and failed result writing must release the bounded pool."""
    _use_hmm_context(monkeypatch, "spawn")
    monkeypatch.setattr(_hmmfilter, "_HMM_BOUNDED_SIGNAL_FRAMES", 1)
    monkeypatch.setattr(_hmmfilter, "_HMM_LEGACY_CHUNK_RECORDS", 4)
    monkeypatch.setattr(_hmmfilter.tempfile, "tempdir", str(tmp_path))
    filter = _filter(tmp_path / "input", 2)

    def inject(nproc, **kwargs):
        """Inject only the selected worker-side fault."""
        if failure != "consumer":
            kwargs["func"] = _failed_record if failure == "worker" else _exited_record
        return utils.run_mp(nproc, **kwargs)

    def failed_write(*args, **kwargs):
        """Fail consumption after the pool has started."""
        raise RuntimeError("intentional HMM consumer failure")

    monkeypatch.setattr(_hmmfilter, "run_mp", inject)
    if failure == "consumer":
        monkeypatch.setattr(utils.WriteBuffer, "extend", failed_write)
    before = {p.pid for p in multiprocessing.active_children()}
    with pytest.raises(RuntimeError, match=r"intentional|exited"):
        filter._calhmm()
    assert {p.pid for p in multiprocessing.active_children()} == before
    assert Path(filter.moleculetempfilename).exists()
    assert not list(tmp_path.glob("reacnetgenerator-order-*"))


@pytest.mark.parametrize("workers", [4, 32])
def test_one_legacy_chunk_of_long_signals_runs_locally(tmp_path, monkeypatch, workers):
    """Few records must not activate more HMM working arrays than the old chunk."""
    monkeypatch.setattr(_hmmfilter, "_HMM_BOUNDED_SIGNAL_FRAMES", 64)
    monkeypatch.setattr(_hmmfilter.tempfile, "tempdir", str(tmp_path))
    serial = _filter(tmp_path / "serial", 1)
    selected = _filter(tmp_path / "selected", workers)

    def unexpected(*args, **kwargs):
        """Reject unnecessary pool creation for one legacy chunk."""
        raise AssertionError("One legacy chunk must use the local HMM path")

    monkeypatch.setattr(_hmmfilter, "run_mp", unexpected)
    serial._calhmm()
    selected._calhmm()
    assert _read(selected) == _read(serial)
    assert selected.nproc == workers


def test_long_hmm_tasks_do_not_resend_parent_metadata(tmp_path, monkeypatch):
    """Neither worker initialization nor individual tasks carry the RNG context."""
    _use_hmm_context(monkeypatch, "spawn")
    monkeypatch.setattr(_hmmfilter, "_HMM_BOUNDED_SIGNAL_FRAMES", 64)
    monkeypatch.setattr(_hmmfilter, "_HMM_LEGACY_CHUNK_RECORDS", 4)
    monkeypatch.setattr(_hmmfilter.tempfile, "tempdir", str(tmp_path))
    filter = _filter(tmp_path / "input", 2)
    filter.rng = SimpleNamespace(metadata=b"x" * (1024 * 1024))
    original = utils._submit_pool_item
    sizes = []
    initialization_sizes = []

    def capture(nproc, **kwargs):
        """Measure the actual initializer envelope, including its model."""
        initialization_sizes.append(
            len(ForkingPickler.dumps((kwargs["initializer"], kwargs["initargs"])))
        )
        return utils.run_mp(nproc, **kwargs)

    monkeypatch.setattr(_hmmfilter, "run_mp", capture)

    def submit(pool, completed, pending, function, index, item):
        """Measure the actual worker call envelope before its submission."""
        sizes.append(len(ForkingPickler.dumps((function, item))))
        return original(pool, completed, pending, function, index, item)

    monkeypatch.setattr(utils, "_submit_pool_item", submit)
    filter._calhmm()
    assert filter.hmmit == len(sizes) == 12
    assert max(sizes) < 32768
    assert len(initialization_sizes) == 1
    assert initialization_sizes[0] < 8192


@pytest.mark.parametrize("method", ["spawn", "forkserver"])
def test_hmm_worker_preserves_model_without_parent_context(
    tmp_path, monkeypatch, method
):
    """Keep changed prediction parameters, while excluding unpicklable RNG data."""
    if method not in multiprocessing.get_all_start_methods():
        pytest.skip(f"{method} unavailable")
    _use_hmm_context(monkeypatch, method)
    monkeypatch.setattr(_hmmfilter, "_HMM_BOUNDED_SIGNAL_FRAMES", 64)
    monkeypatch.setattr(_hmmfilter, "_HMM_LEGACY_CHUNK_RECORDS", 4)
    monkeypatch.setattr(_hmmfilter.tempfile, "tempdir", str(tmp_path))
    default = _filter(tmp_path / "default", 1)
    serial = _filter(tmp_path / "serial", 1)
    parallel = _filter(tmp_path / "parallel", 2)
    for filter in (serial, parallel):
        # These values deliberately differ from the original p/a/b settings.
        filter._model.algorithm = "map"
        filter._model.emissionprob_ = filter._model.emissionprob_[::-1].copy()
        filter.rng = SimpleNamespace(unrelated_callback=lambda: None)
    default._calhmm()
    serial._calhmm()
    parallel._calhmm()
    assert _read(default) != _read(serial)
    assert _read(serial) == _read(parallel)


@pytest.mark.parametrize("method", ["spawn", "forkserver"])
@pytest.mark.parametrize("long_signal", [False, True])
def test_consumer_failure_closes_hmm_pool_with_retained_traceback(
    tmp_path, monkeypatch, method, long_signal
):
    """Keeping a write-error traceback alive must not keep pool workers alive."""
    if method not in multiprocessing.get_all_start_methods():
        pytest.skip(f"{method} unavailable")
    _use_hmm_context(monkeypatch, method)
    monkeypatch.setattr(
        _hmmfilter, "_HMM_BOUNDED_SIGNAL_FRAMES", 64 if long_signal else 256
    )
    monkeypatch.setattr(_hmmfilter, "_HMM_LEGACY_CHUNK_RECORDS", 4)
    monkeypatch.setattr(_hmmfilter.tempfile, "tempdir", str(tmp_path))
    filter = _filter(tmp_path / "input", 2)

    def failed_write(*args, **kwargs):
        """Fail definition output after the pool has delivered a result."""
        raise OSError("intentional HMM output failure")

    monkeypatch.setattr(utils.WriteBuffer, "extend", failed_write)
    before = {child.pid for child in multiprocessing.active_children()}
    iterators = []

    def track_results(*args, **kwargs):
        """Keep failed-test cleanup under the iterator's own pool lifecycle."""
        results = utils.run_mp(*args, **kwargs)
        iterators.append(results)
        return results

    monkeypatch.setattr(_hmmfilter, "run_mp", track_results)
    try:
        with pytest.raises(OSError, match="intentional HMM output failure") as error:
            filter._calhmm()
        # Retain the exception and its _calhmm frame through the worker check.
        # Otherwise garbage collection can hide missing iterator cleanup.
        assert error.value.__traceback__ is not None
        assert {child.pid for child in multiprocessing.active_children()} == before
        assert Path(filter.moleculetempfilename).exists()
    finally:
        # Clean up after a red regression test without terminating individual
        # workers that a still-running Pool would otherwise replace.
        for results in iterators:
            results.close()
