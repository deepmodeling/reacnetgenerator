# SPDX-License-Identifier: LGPL-3.0-or-later
"""Bound ordered SMILES batches while keeping converter state in each worker."""

import io
import os
import pickle
import tempfile
from contextlib import closing
from multiprocessing import get_start_method

import lz4.frame
from tqdm.auto import tqdm

from ._moleculepresence import _iter_molecule_records
from .utils import run_mp

_SMILES_BATCH_RECORDS = 64
_SMILES_INPUT_BYTES = 256 * 1024
_SMILES_DECODED_ESTIMATE_BYTES = 4 * 1024 * 1024
_SMILES_RESULT_BYTES = 1024 * 1024
_SMILES_SPAWN_RECORDS_PER_WORKER = 8192
_SMILES_FORK_RECORDS_PER_WORKER = 256
_SMILES_WORKER_STATE = None


def _should_batch_smiles(molecule_count, nproc):
    """Retain ordered per-record scheduling for small startup-sensitive inputs."""
    threshold = (
        _SMILES_FORK_RECORDS_PER_WORKER
        if get_start_method() == "fork"
        else _SMILES_SPAWN_RECORDS_PER_WORKER
    )
    return nproc > 1 and molecule_count >= threshold * nproc


def _smiles_record_cost(record):
    """Estimate decoded storage without expanding compressed molecule fields.

    Pickle content lengths are a planning proxy, not an RSS guarantee. A large
    or unknown-size record is isolated; workers and consumers decode only one
    result at a time, and actual output bytes are separately capped/spilled.
    """
    decoded = 0
    for block in record:
        size = lz4.frame.get_frame_info(block[64:])["content_size"]
        if not size:
            return sum(map(len, record)), _SMILES_DECODED_ESTIMATE_BYTES + 1
        decoded += size
    # Allow for graph lists, Python integers, names and temporary NumPy copies.
    return sum(map(len, record)), 64 * decoded + 4096


def _smiles_batches(molecule_file):
    """Plan finite batches, with a singleton exception for indivisible records."""
    batch = []
    # At most 256 byte strings, tuple headers and pickle frames fit this reserve.
    compressed = 4096
    decoded = 0
    for record in _iter_molecule_records(molecule_file):
        record_bytes, record_decoded = _smiles_record_cost(record)
        if batch and (
            len(batch) >= _SMILES_BATCH_RECORDS
            or compressed + record_bytes > _SMILES_INPUT_BYTES
            or decoded + record_decoded > _SMILES_DECODED_ESTIMATE_BYTES
        ):
            yield batch
            batch, compressed, decoded = [], 4096, 0
        batch.append(record)
        compressed += record_bytes
        decoded += record_decoded
        if compressed > _SMILES_INPUT_BYTES or decoded > _SMILES_DECODED_ESTIMATE_BYTES:
            # Do not accumulate another record alongside a large input.
            yield batch
            batch, compressed, decoded = [], 4096, 0
    if batch:
        yield batch


def _encoded_smiles_batches(molecule_file):
    """Verify serialized task bytes before handing them to the shared runner."""
    for index, batch in enumerate(_smiles_batches(molecule_file)):
        payload = pickle.dumps(batch, protocol=pickle.HIGHEST_PROTOCOL)
        if len(payload) > _SMILES_INPUT_BYTES and len(batch) != 1:
            raise ValueError("SMILES batch exceeds its input byte budget")
        yield index, payload


def _initialize_smiles_worker(
    atomname, atomnames, need_frames, frequency, folder, result_bytes
):
    """Install only conversion metadata; no RNG object travels with each task."""
    from ._path import _CollectSMILESPaths

    global _SMILES_WORKER_STATE
    converter = object.__new__(_CollectSMILESPaths)
    converter.atomname = atomname
    converter.atomnames = atomnames
    converter.printmoleculetime = need_frames
    converter._moleculeframefilter = None
    converter._moleculetimestepfilter = None
    function = (
        converter._calmoleculeSMILESfrequency
        if frequency
        else converter._calmoleculeSMILESname
    )
    _SMILES_WORKER_STATE = function, folder, result_bytes


class _SMILESResultWriter:
    """Keep at most a fixed serialized prefix, spilling before a write exceeds it.

    Independent pickle records let the parent decode large results one at a
    time. The parent owns the directory until all workers have stopped, so
    interrupted workers cannot leave orphaned result files.
    """

    def __init__(self, path, limit):
        self.path = path
        self.limit = limit
        self.memory = io.BytesIO()
        self.file = None

    def write(self, data):
        """Stream oversized writes directly to disk rather than growing RAM."""
        if isinstance(data, pickle.PickleBuffer):
            # Protocol 5 passes large NumPy storage as a buffer, not bytes.
            data = data.raw()
        if self.file is None and self.memory.tell() + len(data) > self.limit:
            self.file = open(self.path, "wb")
            self.file.write(self.memory.getvalue())
            self.memory.close()
        if self.file is not None:
            return self.file.write(data)
        return self.memory.write(data)

    def result(self):
        """Return bounded bytes, or a marker for the completed disk stream."""
        return self.memory.getvalue() if self.file is None else None

    def close(self):
        """Close both normal and exceptional writes before returning a result."""
        if self.file is not None:
            self.file.close()
        self.memory.close()


def _convert_smiles_batch(task):
    """Convert without retaining a batch of decoded graphs or large IPC results."""
    assert _SMILES_WORKER_STATE is not None
    function, folder, result_bytes = _SMILES_WORKER_STATE
    index, payload = task
    records = pickle.loads(payload)
    if len(payload) > _SMILES_INPUT_BYTES and len(records) != 1:
        raise ValueError("SMILES batch exceeds its input byte budget")
    path = os.path.join(folder, f"{index}.pickle")
    with closing(_SMILESResultWriter(path, result_bytes)) as output:
        for record in records:
            # A fresh pickle per record avoids retaining prior decoded objects
            # in a Pickler memo and allows one-at-a-time parent consumption.
            pickle.dump(function(record), output, protocol=pickle.HIGHEST_PROTOCOL)
        return index, len(records), output.result()


def _run_smiles_batches(collector, molecule_file, frequency):
    """Use exact ordered task counts and close workers before removing streams.

    Two batches per worker may be submitted but not consumed. Normal tasks
    carry <=256 KiB input and <=1 MiB output, plus runner framing; larger output
    uses a disk reference. Only indivisible input records may exceed the input
    and decoded planning budgets. Converter metadata and individual decoded
    molecules remain proportional to the trajectory, not to these byte caps.
    """
    position = molecule_file.tell()
    batch_count = record_count = 0
    try:
        # The existing ordered runner needs an exact task count. Re-scan the
        # compressed stream instead of retaining a trajectory-sized batch index.
        for batch in _smiles_batches(molecule_file):
            batch_count += 1
            record_count += len(batch)
    finally:
        molecule_file.seek(position)
    if batch_count:
        del batch
    if record_count != collector.hmmit:
        raise ValueError("SMILES record count differs from the declared count")
    with (
        tempfile.TemporaryDirectory(prefix="reacnetgenerator-smiles-") as folder,
        closing(
            run_mp(
                collector.nproc,
                func=_convert_smiles_batch,
                l=_encoded_smiles_batches(molecule_file),
                initializer=_initialize_smiles_worker,
                initargs=(
                    collector.atomname,
                    collector.atomnames,
                    collector._needmoleculetimeline(),
                    frequency,
                    folder,
                    _SMILES_RESULT_BYTES,
                ),
                unordered=False,
                disk_ordered=True,
                max_inflight=2 * collector.nproc,
                chunksize=1,
                total=batch_count,
                bar=False,
            )
        ) as results,
        tqdm(
            total=collector.hmmit,
            desc="Indentify isomers",
            unit="molecule",
            disable=None,
        ) as progress,
    ):
        for expected_index, (index, count, payload) in enumerate(results):
            if index != expected_index:
                raise ValueError("SMILES batches are out of order")
            path = os.path.join(folder, f"{index}.pickle")
            if payload is not None and len(payload) > _SMILES_RESULT_BYTES:
                raise ValueError("SMILES batch exceeds its result byte budget")
            with open(path, "rb") if payload is None else io.BytesIO(payload) as stream:
                for _ in range(count):
                    result = pickle.load(stream)
                    progress.update()
                    yield result
                    del result
                if stream.read(1):
                    raise ValueError("Extra SMILES results in a batch")
            if payload is None:
                os.unlink(path)
