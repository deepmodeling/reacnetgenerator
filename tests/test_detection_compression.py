# SPDX-License-Identifier: LGPL-3.0-or-later
"""Bound compression payloads and ordered delivery without losing observations."""

import multiprocessing
import pickle
import time
from array import array
from contextlib import closing
from pathlib import Path

import pytest

from reacnetgenerator import _detect as detection
from reacnetgenerator import utils
from reacnetgenerator.utils import bytestolist, listtobytes, read_compressed_block

from .test_detection_records import _detector


def _values():
    """Mix widths and retain duplicate, unsorted visits in the wire round trip."""
    return [array("H", [7, 0, 7]), array("I", [2**16, 3]), array("Q", [2**40, 1])]


def _slow_first_batch(payload):
    """Make later real-worker results arrive before the first batch."""
    if pickle.loads(payload)[0][0] == 0:
        time.sleep(0.3)
    return detection._compress_frame_batch_original(payload)


@pytest.mark.parametrize("method", ["spawn", "forkserver", "fork"])
@pytest.mark.parametrize("nproc", [1, 2])
def test_parallel_selection_preserves_small_dense_and_spawn_paths(
    monkeypatch, method, nproc
):
    """Only sufficiently many short records under fork start compression workers."""
    monkeypatch.setattr(detection, "get_start_method", lambda: method)
    monkeypatch.setattr(detection, "_DETECT_COMPRESS_MIN_RECORDS", 4)
    short = [array("H", [i]) for i in range(4)]
    selected = detection._compression_batch_size(short, nproc)
    assert bool(selected) == (method == "fork" and nproc > 1)
    assert detection._compression_batch_size(short[:3], nproc) == 0
    long = [*short, array("H", [0]) * 513]
    assert detection._compression_batch_size(long, nproc) == 0
    assert detection._compression_batch_size([], nproc) == 0


def test_batches_bound_raw_input_serialized_input_and_serialized_output(monkeypatch):
    """Budget actual bytes as well as records, including empty and wide arrays."""
    monkeypatch.setattr(detection, "get_start_method", lambda: "fork")
    monkeypatch.setattr(detection, "_DETECT_COMPRESS_MIN_RECORDS", 0)
    values = [array("H") for _ in range(2300)]
    values.extend(
        array("Q", [i * 65537 + j * 104729 for j in range(4)]) for i in range(2300)
    )
    values.extend(_values())
    size = detection._compression_batch_size(values, 2)
    restored = []
    for payload in detection._frame_batches(values, size):
        assert len(payload) <= detection._DETECT_COMPRESS_INPUT_BYTES
        batch = pickle.loads(payload)
        assert len(batch) <= detection._DETECT_COMPRESS_MAX_RECORDS
        assert (
            sum(len(v) * v.itemsize for v in batch)
            <= detection._DETECT_COMPRESS_RAW_BYTES
        )
        result = detection._compress_frame_batch(payload)
        assert len(result) <= detection._DETECT_COMPRESS_RESULT_BYTES
        restored.extend(bytestolist(block).tolist() for block in pickle.loads(result))
    assert restored == [list(v) for v in values]


def test_byte_budget_reduces_batch_count_for_larger_records(monkeypatch):
    """A fixed task count alone must not determine the input-memory bound."""
    monkeypatch.setattr(detection, "get_start_method", lambda: "fork")
    monkeypatch.setattr(detection, "_DETECT_COMPRESS_MIN_RECORDS", 0)
    small = [array("H", [1])] * 2000
    large = [array("Q", range(128))] * 2000
    assert detection._compression_batch_size(
        large, 2
    ) < detection._compression_batch_size(small, 2)
    for payload in detection._frame_batches(
        large, detection._compression_batch_size(large, 2)
    ):
        assert sum(len(v) * v.itemsize for v in pickle.loads(payload)) <= 65536


@pytest.mark.parametrize("budget", ["input", "result"])
def test_excess_serialized_payload_fails_before_queueing(monkeypatch, budget):
    """Unexpected serializer growth must fail instead of silently exceeding the cap."""
    if budget == "input":
        monkeypatch.setattr(detection, "_DETECT_COMPRESS_INPUT_BYTES", 1)
        with pytest.raises(ValueError, match=r"input.*budget"):
            next(detection._frame_batches(_values(), 3))
    else:
        payload = next(detection._frame_batches(_values(), 3))
        monkeypatch.setattr(detection, "_DETECT_COMPRESS_RESULT_BYTES", 1)
        with pytest.raises(ValueError, match=r"result.*budget"):
            detection._compress_frame_batch(payload)


def test_serial_fallback_keeps_single_large_record_lazy(tmp_path, monkeypatch):
    """An indivisible oversized record stays in the serial path and starts no pool."""
    detector = _detector(tmp_path)
    monkeypatch.setattr(detection, "get_start_method", lambda: "fork")
    monkeypatch.setattr(detection, "_DETECT_COMPRESS_MIN_RECORDS", 0)
    monkeypatch.setattr(
        detection, "run_mp", lambda *a, **k: pytest.fail("Unexpected pool")
    )
    values = [array("H", [7, 1, 7]) * 1000]
    with closing(detector._compressed_values(values)) as results:
        assert [bytestolist(block).tolist() for block in results] == [list(values[0])]


@pytest.mark.parametrize("failure", ["compression", "write"])
def test_parallel_failure_closes_results_and_removes_partial_file(
    tmp_path, monkeypatch, failure
):
    """Closing the detection path must cancel suspended compression on every failure."""
    detector = _detector(tmp_path)
    monkeypatch.setattr(detection, "get_start_method", lambda: "fork")
    monkeypatch.setattr(detection, "_DETECT_COMPRESS_MIN_RECORDS", 0)
    monkeypatch.setattr(detection, "_DETECT_COMPRESS_MAX_RECORDS", 1)
    monkeypatch.setattr(detection.tempfile, "tempdir", str(tmp_path))
    definitions = [b"".join(listtobytes(v) for v in ([i], [], [])) for i in (0, 1)]
    closed = []
    retained = []

    def run(nproc, **kwargs):
        if kwargs["desc"] != "Save molecules":
            return iter_results()
        assert kwargs["max_inflight"] == 4 and kwargs["disk_ordered"]

        def results():
            try:
                yield detection._compress_frame_batch(next(iter(kwargs["l"])))
                raise RuntimeError("failed compression")
            finally:
                closed.append("compression")

        result = results()
        retained.append(result)
        return result

    def iter_results():
        yield definitions, (0, 0)

    if failure == "write":

        def failed_write(self, items):
            raise OSError("failed write")

        monkeypatch.setattr(detection.WriteBuffer, "extend", failed_write)
    monkeypatch.setattr(detection, "run_mp", run)
    previous = tmp_path / "previous.bin"
    previous.write_bytes(b"published")
    detector.moleculetempfilename = str(previous)
    with pytest.raises(RuntimeError if failure == "compression" else OSError):
        detector.detect()
    assert closed == ["compression"]
    assert detector.moleculetempfilename == str(previous)
    assert previous.read_bytes() == b"published"
    assert list(tmp_path.iterdir()) == [previous]


@pytest.mark.skipif(
    multiprocessing.get_start_method() != "fork", reason="real fork scheduling test"
)
def test_real_out_of_order_batches_are_bounded_and_written_progressively(
    tmp_path, monkeypatch
):
    """Delay batch zero to fill the window; verify ordered output and bounded spooling."""
    detector = _detector(tmp_path)
    monkeypatch.setattr(detection, "_DETECT_COMPRESS_MIN_RECORDS", 0)
    monkeypatch.setattr(detection, "_DETECT_COMPRESS_MAX_RECORDS", 2)
    monkeypatch.setattr(
        detection,
        "_compress_frame_batch_original",
        detection._compress_frame_batch,
        raising=False,
    )
    monkeypatch.setattr(detection, "_compress_frame_batch", _slow_first_batch)
    values = [array("H", [i, 0, i]) for i in range(30)]
    original_batches = detection._frame_batches
    submitted, written, maximum_pending = 0, 0, 0

    def batches(*args):
        """Count submissions against fully written batches, not completed workers."""
        nonlocal submitted
        for payload in original_batches(*args):
            submitted += 1
            assert submitted - written // 2 <= 4
            yield payload

    class Spool(utils._DiskOrderedResultSpool):
        def put(self, index, result):
            nonlocal maximum_pending
            super().put(index, result)
            maximum_pending = max(maximum_pending, self.pending_count)
            assert self.pending_count <= 3
            assert sum(Path(p).stat().st_size for p in self._paths.values()) <= 3 * (
                detection._DETECT_COMPRESS_RESULT_BYTES + 1024
            )

    class Writer(detection.WriteBuffer):
        def extend(self, items):
            nonlocal written
            assert self.maxbytes == 1024 * 1024
            written += 1
            super().extend(items)

    monkeypatch.setattr(detection, "_frame_batches", batches)
    monkeypatch.setattr(utils, "_DiskOrderedResultSpool", Spool)
    monkeypatch.setattr(detection, "WriteBuffer", Writer)
    with closing(detector._compressed_values(values)) as results:
        detector._writemoleculetempfile(([b""] * len(values), results))
    path = Path(detector.moleculetempfilename)
    with path.open("rb") as f:
        restored = [bytestolist(block).tolist() for block in read_compressed_block(f)]
    assert restored == [list(value) for value in values]
    assert submitted == 15 and written == 30 and maximum_pending > 0
    path.unlink()


def test_default_parallel_threshold_is_conservative(monkeypatch):
    """Keep inputs below the measured 24,576-record control on the serial path."""
    monkeypatch.setattr(detection, "get_start_method", lambda: "fork")
    values = [array("H", [0])] * 24576
    assert detection._compression_batch_size(values[:-1], 2) == 0
    assert detection._compression_batch_size(values, 2) > 0


def test_batch_reuses_only_exact_histories_without_dropping_records(monkeypatch):
    """Reuse bytes within a batch while preserving dtype, order and multiplicity."""
    values = [
        array("H", [2, 0, 2]),
        array("H", [2, 0, 2]),
        array("H", [2, 2, 0]),
        array("H", [2, 0]),
        array("Q", [2, 0, 2]),
    ]
    calls = []
    original = detection.listtobytes

    def compress(value):
        calls.append((value.dtype.str, value.tolist()))
        return original(value)

    monkeypatch.setattr(detection, "listtobytes", compress)
    payload = next(detection._frame_batches(values, len(values)))
    blocks = pickle.loads(detection._compress_frame_batch(payload))
    assert len(calls) == 4
    assert len(blocks) == len(values)
    for block, value in zip(blocks, values, strict=True):
        restored = bytestolist(block)
        assert restored.itemsize == value.itemsize
        assert restored.tolist() == list(value)


@pytest.mark.parametrize("parallel", [False, True])
def test_collection_widens_before_serial_or_parallel_saving(
    tmp_path, monkeypatch, parallel
):
    """Exercise append overflow through collection, retaining arrival order and repeats."""
    detector = _detector(tmp_path)
    monkeypatch.setattr(
        detection, "get_start_method", lambda: "fork" if parallel else "spawn"
    )
    monkeypatch.setattr(detection, "_DETECT_COMPRESS_MIN_RECORDS", 0)
    definition = b"".join(listtobytes(v) for v in ([0], [], []))
    frames = [65535, 65536, 0, 2**32, 65535]

    def run(nproc, **kwargs):
        """Inject width-crossing arrivals while executing the selected save path."""
        if kwargs["desc"] == "Save molecules":
            yield from map(kwargs["func"], kwargs["l"])
        else:
            for frame in frames:
                yield [definition], (frame, frame * 10)

    monkeypatch.setattr(detection, "run_mp", run)
    detector.detect()
    path = Path(detector.rng.moleculetempfilename)
    with path.open("rb") as stream:
        blocks = list(read_compressed_block(stream))
    restored = bytestolist(blocks[3])
    assert restored.dtype.kind == "u" and restored.itemsize == 8
    assert restored.tolist() == frames
    assert detector.rng.step == 4
    path.unlink()
