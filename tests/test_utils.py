# SPDX-License-Identifier: LGPL-3.0-or-later
"""Test multiprocessing utilities."""

import os
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
