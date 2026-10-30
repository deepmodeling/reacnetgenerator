# SPDX-License-Identifier: LGPL-3.0-or-later
"""Preserve detection records while bounding compression and write buffering."""

from array import array
from pathlib import Path

import numpy as np
import pytest

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator import _detect as detection
from reacnetgenerator.utils import (
    WriteBuffer,
    bytestolist,
    listtobytes,
    read_compressed_block,
)

INPUTS = Path(__file__).parent / "inputs"


def _detector(tmp_path, **kwargs):
    rng = ReacNetGenerator(
        inputfilename=INPUTS / "water.bond",
        inputfiletype="bond",
        atomname=["H", "O"],
        output_dir=tmp_path,
        timed_output=tmp_path / "timeline.h5",
        nproc=2,
        **kwargs,
    )
    return detection._Detect.gettype(rng)


def test_detection_streams_compact_records_in_first_seen_order(tmp_path, monkeypatch):
    """Writing must begin before later compression, preserving arrival order."""
    monkeypatch.setattr(detection.tempfile, "tempdir", str(tmp_path))
    detector = _detector(tmp_path)
    definitions = [
        b"".join(listtobytes(x) for x in ([atom], [], [])) for atom in (0, 1)
    ]
    calls, writes, frame_arrays = [], [], []

    def results(nproc, **kwargs):
        calls.append((nproc, kwargs))
        if kwargs["desc"] == "Save molecules":
            yield from map(kwargs["func"], kwargs["l"])
            return
        for frame, ids in [(2, [1, 0]), (0, [0, 1]), (1, [1, 1])]:
            yield [definitions[index] for index in ids], (frame, frame * 10)

    original_compress = detector._compressvalue

    def compress(values):
        if frame_arrays:
            assert writes, "All compression results were accumulated before writing"
        assert isinstance(values, array)
        frame_arrays.append(list(values))
        assert values.itemsize * len(values) <= 2 * len(values)
        return original_compress(values)

    class ObservedBuffer(WriteBuffer):
        def extend(self, items):
            writes.append(True)
            super().extend(items)

    monkeypatch.setattr(detection, "run_mp", results)
    monkeypatch.setattr(detector, "_compressvalue", compress)
    monkeypatch.setattr(detection, "WriteBuffer", ObservedBuffer)
    detector.detect()
    assert len(calls) == 1  # Compression does not start a second process pool.
    assert calls[0][0] == 2 and calls[0][1]["max_inflight"] == 4
    assert frame_arrays == [[2, 0, 1, 1], [2, 0]]
    assert detector.rng.timestep == {2: 20, 0: 0, 1: 10}
    assert detector.rng.step == 3 and detector.rng.temp1it == 2
    with open(detector.rng.moleculetempfilename, "rb") as source:
        blocks = list(read_compressed_block(source))
    assert b"".join(blocks[:3]) == definitions[1]
    assert b"".join(blocks[4:7]) == definitions[0]
    assert [bytestolist(blocks[index]).tolist() for index in (3, 7)] == frame_arrays
    Path(detector.rng.moleculetempfilename).unlink()


@pytest.mark.parametrize("frame", [2**16, 2**32, 2**64 - 1])
def test_frame_width_expands_without_reordering_or_loss(frame):
    """Sparse jumps and boundary values survive compression without huge fixtures."""
    before = array("H", [7, 0, 7, 65535])
    widened = detection._widen_frame_indices(before, frame)
    assert list(before) == [7, 0, 7, 65535]
    assert list(widened) == [7, 0, 7, 65535, frame]
    assert widened.itemsize > before.itemsize
    restored = bytestolist(detection._Detect._compressvalue(widened))
    assert restored.dtype.kind == "u"
    assert restored.tolist() == list(widened)


@pytest.mark.parametrize("frame", [-1, 2**64])
def test_unrepresentable_frames_do_not_wrap(frame):
    """An unsupported frame must fail without changing earlier observations."""
    before = array("Q", [3, 1])
    with pytest.raises(OverflowError):
        detection._widen_frame_indices(before, frame)
    assert list(before) == [3, 1]


@pytest.mark.parametrize("failure", ["compression", "write", "short", "long"])
def test_failed_save_closes_and_removes_partial_output(tmp_path, monkeypatch, failure):
    """A partial or misaligned definition stream must never become an artifact."""
    detector = _detector(tmp_path)
    handles = []
    create = detection.tempfile.NamedTemporaryFile

    def temporary(*args, **kwargs):
        handle = create(*args, dir=tmp_path, **kwargs)
        handles.append(handle)
        return handle

    def values():
        yield b"first"
        if failure == "compression":
            raise RuntimeError("failed compression")
        if failure != "short":
            yield b"second"
        if failure == "long":
            yield b"extra"

    def failed_flush(self):
        raise OSError("failed write")

    monkeypatch.setattr(detection.tempfile, "NamedTemporaryFile", temporary)
    if failure == "write":
        monkeypatch.setattr(detection.WriteBuffer, "flush", failed_flush)
    error = {"compression": RuntimeError, "write": OSError}.get(failure, ValueError)
    with pytest.raises(error):
        detector._writemoleculetempfile(([b"A", b"B"], values()))
    assert detector.moleculetempfilename is None
    assert all(handle.closed and not Path(handle.name).exists() for handle in handles)


def test_collection_failure_closes_workers_before_source(tmp_path, monkeypatch):
    """An accumulation error must close a suspended worker-result generator."""
    detector = _detector(tmp_path)
    closed, retained = [], []

    def lines(_):
        try:
            yield "line"
        finally:
            closed.append("source")

    def run(nproc, **kwargs):
        def results():
            try:
                next(kwargs["l"])
                yield [[]], (0, 0)  # Invalid key forces a parent-side failure.
            finally:
                closed.append("workers")

        result = results()
        retained.append(result)  # Prevent reference counting from hiding a leak.
        return result

    monkeypatch.setattr(detector, "_source_lines", lines)
    monkeypatch.setattr(detection, "run_mp", run)
    with pytest.raises(TypeError):
        detector.detect()
    assert closed == ["workers", "source"]


@pytest.mark.parametrize("key_size", [5, 20])
def test_save_buffer_is_bounded_by_payload_bytes(tmp_path, monkeypatch, key_size):
    """Multiple records cannot accumulate past the budget; one payload may exceed it."""
    detector = _detector(tmp_path)
    create = detection.tempfile.NamedTemporaryFile
    chunks = []

    class File:
        def __init__(self, handle):
            self.handle = handle

        def __getattr__(self, name):
            return getattr(self.handle, name)

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return self.handle.__exit__(*args)

        def writelines(self, values):
            data = b"".join(values)
            assert len(data) <= max(16, key_size)
            chunks.append(data)
            self.handle.write(data)

    monkeypatch.setattr(detection, "_DETECT_WRITE_BYTES", 16)
    monkeypatch.setattr(
        detection.tempfile,
        "NamedTemporaryFile",
        lambda *args, **kwargs: File(create(*args, dir=tmp_path, **kwargs)),
    )
    detector._writemoleculetempfile(([b"a" * key_size] * 5, iter([b"b" * 7] * 5)))
    path = Path(detector.moleculetempfilename)
    assert path.read_bytes() == (b"a" * key_size + b"b" * 7) * 5 == b"".join(chunks)
    assert len(chunks) > 1
    path.unlink()


def test_no_detected_molecules_writes_empty_record_file(tmp_path, monkeypatch):
    """Empty molecule collections keep frame metadata and need no save worker."""
    detector = _detector(tmp_path)
    calls = []

    def results(nproc, **kwargs):
        calls.append(kwargs["desc"])
        yield [], (0, 99)

    monkeypatch.setattr(detection, "run_mp", results)
    detector.detect()
    path = Path(detector.rng.moleculetempfilename)
    assert path.read_bytes() == b""
    assert detector.rng.timestep == {0: 99} and detector.rng.temp1it == 0
    assert len(calls) == 1
    path.unlink()


@pytest.mark.parametrize("nproc", [1, 2])
@pytest.mark.parametrize("filetype", ["bond", "dump", "xyz", "extxyz"])
def test_real_detection_formats_keep_definition_and_frame_alignment(
    tmp_path, filetype, nproc
):
    """Real workers produce four-field molecule records in every supported format."""
    rng = ReacNetGenerator(
        inputfilename=INPUTS / f"water.{filetype}",
        inputfiletype=filetype,
        atomname=["H", "O"],
        output_dir=tmp_path,
        nproc=nproc,
        pbc=False,
        timed_output=tmp_path / "timeline.h5",
    )
    detection._Detect.gettype(rng).detect()
    try:
        with open(rng.moleculetempfilename, "rb") as source:
            blocks = list(read_compressed_block(source))
        assert len(blocks) == 4 and rng.temp1it == 1 and rng.N == 3
        np.testing.assert_array_equal(rng.atomtype, [0, 0, 1])
        assert sorted(bytestolist(blocks[0])) == [0, 1, 2]
        frames = bytestolist(blocks[3])
        assert frames.dtype.kind == "u"
        assert sorted(frames.tolist()) == list(range(rng.step))
        assert rng.source_frame_counts == [rng.step]
    finally:
        Path(rng.moleculetempfilename).unlink()
