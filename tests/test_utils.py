# SPDX-License-Identifier: LGPL-3.0-or-later
"""Test multiprocessing utilities."""

import multiprocessing
import os
import socket
import time

import pytest

from reacnetgenerator.utils import run_mp


class _CountingIterator:
    def __init__(self, stop):
        self.current = 0
        self.stop = stop
        self.produced = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.current == self.stop:
            raise StopIteration
        value = self.current
        self.current += 1
        self.produced += 1
        return value


def _identity(value):
    return value


def _finish_out_of_order(value):
    if value == 0:
        time.sleep(0.2)
    else:
        time.sleep(0.01)
    return value, time.monotonic()


def _fail_out_of_order(value):
    if value == 0:
        time.sleep(0.2)
    if value == 2:
        raise RuntimeError("intentional worker failure")
    return value


def _exit_worker(value):
    if value == 0:
        os._exit(7)
    return value


def _exit_worker_cleanly(value):
    if value == 0:
        os._exit(0)
    return value


def _raise_system_exit(value):
    if value == 0:
        raise SystemExit
    return value


class _ExitWhilePickling:
    def __reduce__(self):
        os._exit(0)


class _RaiseSystemExitWhilePickling:
    def __reduce__(self):
        raise SystemExit


class _ExitWhilePicklingError(Exception):
    def __reduce__(self):
        os._exit(0)


class _RaiseSystemExitWhilePicklingError(Exception):
    def __reduce__(self):
        raise SystemExit


def _return_exit_while_pickling(value):
    return _ExitWhilePickling()


def _return_system_exit_while_pickling(value):
    return _RaiseSystemExitWhilePickling()


def _raise_exit_while_pickling(value):
    raise _ExitWhilePicklingError


def _raise_system_exit_while_pickling(value):
    raise _RaiseSystemExitWhilePicklingError


def _raise_system_exit_when_unpickled():
    raise SystemExit


class _ExitWhenUnpickled:
    def __reduce__(self):
        return os._exit, (0,)


class _RaiseSystemExitWhenUnpickled:
    def __reduce__(self):
        return _raise_system_exit_when_unpickled, ()


def _run_serialization_exit_case(result_queue, func, disk_ordered, spool_dir):
    """Run the potential hang in a watchdog process and report its outcome."""
    try:
        list(
            run_mp(
                2,
                func=func,
                l=range(2),
                unordered=not disk_ordered,
                chunksize=1,
                max_inflight=2,
                disk_ordered=disk_ordered,
                ordered_spool_dir=spool_dir,
                total=2,
                bar=False,
            )
        )
    except BaseException as error:
        result_queue.put((type(error).__name__, str(error)))
    else:
        result_queue.put(("success", ""))


def _run_input_deserialization_exit_case(result_queue, mode, disk_ordered, spool_dir):
    """Exercise worker input deserialization inside an outer watchdog."""
    item = (
        _ExitWhenUnpickled()
        if mode == "os-exit-zero"
        else (_RaiseSystemExitWhenUnpickled())
    )
    try:
        list(
            run_mp(
                2,
                func=_identity,
                l=[item],
                unordered=not disk_ordered,
                chunksize=1,
                max_inflight=1,
                disk_ordered=disk_ordered,
                ordered_spool_dir=spool_dir,
                total=1,
                bar=False,
            )
        )
    except BaseException as error:
        result_queue.put((type(error).__name__, str(error)))
    else:
        result_queue.put(("success", ""))


def test_run_mp_default_ordering_is_unchanged():
    """Keep the existing ordered ``run_mp`` behavior without opt-in arguments."""
    assert list(
        run_mp(
            2,
            func=_identity,
            l=range(6),
            unordered=False,
            total=6,
            bar=False,
        )
    ) == list(range(6))


def test_disk_ordered_run_mp_bounds_inflight_and_restores_order(tmp_path):
    """Bound submitted results while workers finish out of order."""
    source = _CountingIterator(12)
    max_inflight = 4
    results = run_mp(
        2,
        func=_finish_out_of_order,
        l=source,
        unordered=False,
        chunksize=1,
        max_inflight=max_inflight,
        disk_ordered=True,
        ordered_spool_dir=str(tmp_path),
        total=12,
        bar=False,
    )

    records = []
    for record in results:
        records.append(record)
        assert source.produced - len(records) < max_inflight
        if len(records) == 1:
            assert source.produced == max_inflight
            assert list(tmp_path.iterdir())

    assert [record[0] for record in records] == list(range(12))
    assert records[1][1] < records[0][1]
    assert not list(tmp_path.iterdir())


def test_disk_ordered_run_mp_propagates_error_and_cleans_spool(tmp_path):
    """Propagate ordinary worker errors without retaining temporary files."""
    with pytest.raises(RuntimeError, match="intentional worker failure"):
        list(
            run_mp(
                2,
                func=_fail_out_of_order,
                l=range(8),
                unordered=False,
                chunksize=1,
                max_inflight=4,
                disk_ordered=True,
                ordered_spool_dir=str(tmp_path),
                total=8,
                bar=False,
            )
        )

    assert not list(tmp_path.iterdir())


def test_disk_ordered_run_mp_detects_worker_exit_and_cleans_spool(tmp_path):
    """Detect a worker process that exits before reporting its result."""
    with pytest.raises(RuntimeError, match=r"worker .* exited unexpectedly"):
        list(
            run_mp(
                2,
                func=_exit_worker,
                l=range(4),
                unordered=False,
                chunksize=1,
                max_inflight=2,
                disk_ordered=True,
                ordered_spool_dir=str(tmp_path),
                total=4,
                bar=False,
            )
        )

    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize(
    "func",
    [_exit_worker_cleanly, _raise_system_exit],
    ids=["os-exit-zero", "system-exit"],
)
def test_disk_ordered_run_mp_detects_clean_worker_exit(tmp_path, func):
    """Do not hang when a task terminates its worker with a zero exit code."""
    with pytest.raises(RuntimeError, match=r"worker .* exited"):
        list(
            run_mp(
                2,
                func=func,
                l=range(4),
                unordered=False,
                chunksize=1,
                max_inflight=2,
                disk_ordered=True,
                ordered_spool_dir=str(tmp_path),
                total=4,
                bar=False,
            )
        )

    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize(
    "func",
    [
        _return_exit_while_pickling,
        _return_system_exit_while_pickling,
        _raise_exit_while_pickling,
        _raise_system_exit_while_pickling,
    ],
    ids=[
        "result-os-exit-zero",
        "result-system-exit",
        "error-os-exit-zero",
        "error-system-exit",
    ],
)
@pytest.mark.parametrize("disk_ordered", [False, True], ids=["unordered", "ordered"])
def test_run_mp_detects_clean_exit_during_result_serialization(
    tmp_path, func, disk_ordered
):
    """Fail promptly if result serialization terminates a worker cleanly."""
    context = multiprocessing.get_context("spawn")
    result_queue = context.Queue()
    process = context.Process(
        target=_run_serialization_exit_case,
        args=(result_queue, func, disk_ordered, str(tmp_path)),
    )
    process.start()
    process.join(timeout=10)
    if process.is_alive():
        process.terminate()
        process.join()
        pytest.fail("run_mp hung after a clean exit during result serialization")

    assert process.exitcode == 0
    error_type, message = result_queue.get(timeout=2)
    assert error_type == "RuntimeError"
    assert "worker" in message and "exited" in message
    result_queue.close()
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("mode", ["os-exit-zero", "system-exit"])
@pytest.mark.parametrize("disk_ordered", [False, True], ids=["unordered", "ordered"])
def test_run_mp_detects_clean_exit_during_input_deserialization(
    tmp_path, mode, disk_ordered
):
    """Track a task before user-controlled input deserialization begins."""
    context = multiprocessing.get_context("spawn")
    result_queue = context.Queue()
    process = context.Process(
        target=_run_input_deserialization_exit_case,
        args=(result_queue, mode, disk_ordered, str(tmp_path)),
    )
    process.start()
    process.join(timeout=10)
    if process.is_alive():
        process.terminate()
        process.join()
        pytest.fail("run_mp hung after a clean exit during input deserialization")

    assert process.exitcode == 0
    error_type, message = result_queue.get(timeout=2)
    assert error_type == "RuntimeError"
    assert "worker" in message and "exited" in message
    result_queue.close()
    assert not list(tmp_path.iterdir())


def test_bounded_run_mp_preserves_pool_result_reducers():
    """Keep support for types registered only with ForkingPickler."""
    result = next(
        iter(
            run_mp(
                1,
                func=socket.socket,
                l=[socket.AF_INET],
                unordered=True,
                chunksize=1,
                max_inflight=1,
                total=1,
                bar=False,
            )
        )
    )
    try:
        assert isinstance(result, socket.socket)
        assert result.family == socket.AF_INET
    finally:
        result.close()


def test_disk_ordered_run_mp_allows_configured_worker_recycling(tmp_path):
    """Do not mistake ``maxtasksperchild`` recycling for a worker failure."""
    assert list(
        run_mp(
            1,
            func=_identity,
            l=range(1001),
            unordered=False,
            chunksize=1,
            max_inflight=4,
            disk_ordered=True,
            ordered_spool_dir=str(tmp_path),
            total=1001,
            bar=False,
        )
    ) == list(range(1001))
    assert not list(tmp_path.iterdir())


def test_disk_ordered_run_mp_early_close_cleans_spool(tmp_path):
    """Terminate outstanding work and remove the spool when consumption stops."""
    results = run_mp(
        2,
        func=_finish_out_of_order,
        l=range(6),
        unordered=False,
        chunksize=1,
        max_inflight=4,
        disk_ordered=True,
        ordered_spool_dir=str(tmp_path),
        total=6,
        bar=False,
    )

    assert next(results)[0] == 0
    assert list(tmp_path.iterdir())
    results.close()

    assert not list(tmp_path.iterdir())


def test_disk_ordered_run_mp_rejects_inexact_total_and_cleans_spool(tmp_path):
    """Reject inaccurate ordering metadata without retaining temporary files."""
    with pytest.raises(RuntimeError, match="count does not match"):
        list(
            run_mp(
                2,
                func=_identity,
                l=range(3),
                unordered=False,
                chunksize=1,
                max_inflight=2,
                disk_ordered=True,
                ordered_spool_dir=str(tmp_path),
                total=4,
                bar=False,
            )
        )

    assert not list(tmp_path.iterdir())
