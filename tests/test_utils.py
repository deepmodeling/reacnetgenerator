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


def _exit_worker_initializer(_task_events):
    os._exit(9)


def _clean_exit_worker_initializer(_task_events):
    os._exit(0)


def _system_exit_worker_initializer(_task_events):
    raise SystemExit


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


def _sleep_and_return(value):
    time.sleep(0.5)
    return value


class _ExitWhenUnpickled:
    def __reduce__(self):
        return os._exit, (0,)


class _RaiseSystemExitWhenUnpickled:
    def __reduce__(self):
        return _raise_system_exit_when_unpickled, ()


class _SlowWhenUnpickled:
    def __reduce__(self):
        return _sleep_and_return, ("slow",)


class _ExitWhenRepickled:
    """Exit if a parent process serializes a worker result a second time."""

    def __init__(self, armed=False):
        self.armed = armed

    def __reduce__(self):
        if self.armed:
            os._exit(0)
        return type(self), (True,)


def _return_resource_reducer_at_legacy_boundary(index):
    if index == 998:
        return _SlowWhenUnpickled()
    if index == 999:
        return socket.socket()
    return index


def _return_out_of_order_socket(index):
    if index == 0:
        time.sleep(0.5)
        return "first"
    return socket.socket()


def _return_out_of_order_connection(index):
    if index == 0:
        time.sleep(0.5)
        return "first"
    receiver, sender = multiprocessing.Pipe(duplex=False)
    sender.send("from worker")
    sender.close()
    return receiver


def _return_out_of_order_repickle_guard(index):
    if index == 0:
        time.sleep(0.5)
        return "first"
    return _ExitWhenRepickled()


def _worker_pid(value):
    return os.getpid()


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
    with pytest.raises(RuntimeError, match="intentional worker failure") as error_info:
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

    assert type(error_info.value.__cause__).__name__ == "RemoteTraceback"
    assert "_fail_out_of_order" in str(error_info.value.__cause__)
    assert not list(tmp_path.iterdir())


def test_bounded_unordered_run_mp_preserves_remote_traceback():
    """Keep legacy worker traceback diagnostics in bounded unordered mode."""
    with pytest.raises(RuntimeError, match="intentional worker failure") as error_info:
        list(
            run_mp(
                2,
                func=_fail_out_of_order,
                l=range(8),
                unordered=True,
                chunksize=1,
                max_inflight=4,
                total=8,
                bar=False,
            )
        )

    assert type(error_info.value.__cause__).__name__ == "RemoteTraceback"
    assert "_fail_out_of_order" in str(error_info.value.__cause__)


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
    "initializer",
    [
        _exit_worker_initializer,
        _clean_exit_worker_initializer,
        _system_exit_worker_initializer,
    ],
    ids=["nonzero", "clean", "system-exit"],
)
@pytest.mark.parametrize("disk_ordered", [False, True], ids=["unordered", "ordered"])
def test_bounded_run_mp_detects_worker_initializer_exit(
    monkeypatch, tmp_path, initializer, disk_ordered
):
    """Report a worker that exits before it can announce its first task."""
    monkeypatch.setattr(
        "reacnetgenerator.utils._init_bounded_pool_worker",
        initializer,
    )
    with pytest.raises(RuntimeError, match=r"worker .* exited unexpectedly"):
        list(
            run_mp(
                1,
                func=_identity,
                l=[1],
                unordered=not disk_ordered,
                chunksize=1,
                max_inflight=1,
                disk_ordered=disk_ordered,
                ordered_spool_dir=str(tmp_path),
                total=1,
                bar=False,
            )
        )


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


@pytest.mark.parametrize("disk_ordered", [False, True], ids=["unordered", "ordered"])
def test_bounded_run_mp_preserves_reducers_at_legacy_recycle_boundary(
    tmp_path, disk_ordered
):
    """Keep reducer resources alive past the legacy 1000-task boundary."""
    counts = {"integer": 0, "slow": 0, "socket": 0}
    results = run_mp(
        1,
        func=_return_resource_reducer_at_legacy_boundary,
        l=range(1000),
        unordered=not disk_ordered,
        chunksize=1,
        max_inflight=4,
        disk_ordered=disk_ordered,
        ordered_spool_dir=str(tmp_path),
        total=1000,
        bar=False,
    )
    for result in results:
        if isinstance(result, socket.socket):
            counts["socket"] += 1
            result.close()
        elif result == "slow":
            counts["slow"] += 1
        else:
            counts["integer"] += 1

    assert counts == {"integer": 998, "slow": 1, "socket": 1}
    assert not list(tmp_path.iterdir())


def test_disk_ordered_run_mp_spools_socket_reducer_bytes(tmp_path):
    """Spool a genuinely out-of-order socket without re-pickling it."""
    results = iter(
        run_mp(
            2,
            func=_return_out_of_order_socket,
            l=range(2),
            unordered=False,
            chunksize=1,
            max_inflight=2,
            disk_ordered=True,
            ordered_spool_dir=str(tmp_path),
            total=2,
            bar=False,
        )
    )
    assert next(results) == "first"
    result_socket = next(results)
    try:
        assert isinstance(result_socket, socket.socket)
        assert result_socket.family == socket.AF_INET
    finally:
        result_socket.close()
    with pytest.raises(StopIteration):
        next(results)
    assert not list(tmp_path.iterdir())


def test_disk_ordered_run_mp_spools_connection_reducer_bytes(tmp_path):
    """Transfer ownership of an out-of-order Connection exactly once."""
    results = iter(
        run_mp(
            2,
            func=_return_out_of_order_connection,
            l=range(2),
            unordered=False,
            chunksize=1,
            max_inflight=2,
            disk_ordered=True,
            ordered_spool_dir=str(tmp_path),
            total=2,
            bar=False,
        )
    )
    assert next(results) == "first"
    connection = next(results)
    try:
        assert connection.recv() == "from worker"
    finally:
        connection.close()
    with pytest.raises(StopIteration):
        next(results)
    assert not list(tmp_path.iterdir())


def test_disk_ordered_run_mp_does_not_repickle_out_of_order_results(tmp_path):
    """Keep user reduction code out of the parent-side spool write path."""
    context = multiprocessing.get_context("spawn")
    result_queue = context.Queue()
    process = context.Process(
        target=_run_serialization_exit_case,
        args=(result_queue, _return_out_of_order_repickle_guard, True, str(tmp_path)),
    )
    process.start()
    process.join(timeout=10)
    if process.is_alive():
        process.terminate()
        process.join()
        pytest.fail("run_mp hung while spooling an out-of-order result")

    assert process.exitcode == 0
    assert result_queue.get(timeout=2) == ("success", "")
    result_queue.close()
    assert not list(tmp_path.iterdir())


def test_disk_ordered_run_mp_keeps_bounded_worker_alive(tmp_path):
    """Preserve reducer resources by disabling bounded worker recycling."""
    worker_pids = list(
        run_mp(
            1,
            func=_worker_pid,
            l=range(1001),
            unordered=False,
            chunksize=1,
            max_inflight=4,
            disk_ordered=True,
            ordered_spool_dir=str(tmp_path),
            total=1001,
            bar=False,
        )
    )
    assert len(set(worker_pids)) == 1
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
