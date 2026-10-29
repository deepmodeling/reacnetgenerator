# SPDX-License-Identifier: LGPL-3.0-or-later
"""Read effective molecule presence consistently for matrices and timelines."""

import os

import lz4.frame
import numpy as np

from .utils import bytestolist, read_compressed_block

_PRESENCE_BLOCK_ROWS = 65536


def _use_frame_indices(frame_block, frame_count):
    """Reuse indices only when their decoded payload fits a full byte signal.

    LZ4 records include the uncompressed pickle size. The small allowance
    covers array metadata, including empty and short trajectories. Unknown
    sizes conservatively take the signal path. Producer and both readers must
    use this same decision: direct records have no entry in the signal file.
    """
    size = lz4.frame.get_frame_info(frame_block[64:])["content_size"]
    return 0 < size <= frame_count + 256


def _iter_molecule_records(handle, fields=(0, 1, 2, 3)):
    """Read selected fields of complete four-block molecule records.

    Skip bond payloads during matrix construction and all structure payloads
    during no-HMM filtering. Check lengths even for skipped fields so a
    truncated last record cannot silently disappear.
    """
    file_size = os.fstat(handle.fileno()).st_size
    while handle.tell() < file_size:
        record = []
        for field in range(4):
            header = handle.read(64)
            if len(header) != 64:
                raise ValueError("Incomplete molecule record header")
            size = int.from_bytes(header, byteorder="little")
            if handle.tell() + size > file_size:
                raise ValueError("Incomplete molecule record payload")
            if field in fields:
                record.append(header + handle.read(size))
            else:
                handle.seek(size, os.SEEK_CUR)
        yield tuple(record)


def _signal_ranges(values, block_rows=8192):
    """Find closed presence intervals using bounded scan temporaries."""
    if block_rows <= 0:
        raise ValueError("Scan block size must be positive")
    previous = False
    begin = None
    for start in range(0, len(values), block_rows):
        block = np.asarray(values[start : start + block_rows], dtype=np.bool_)
        # HMM/origin records use (frames, 1); callers may also supply 1-D
        # signals. Flatten only a single value per frame, within this block.
        if block.size != len(block):
            raise ValueError("Expected one signal value per frame")
        block = block.reshape(-1)
        changes = np.empty(len(block), dtype=np.bool_)
        changes[0] = block[0] != previous
        changes[1:] = block[1:] != block[:-1]
        for local in np.flatnonzero(changes):
            frame = start + int(local)
            if block[local]:
                begin = frame
            else:
                assert begin is not None
                yield begin, frame - 1
                begin = None
        previous = bool(block[-1])
    if begin is not None:
        yield begin, len(values) - 1


def _index_ranges(values, block_rows=8192):
    """Coalesce sorted unique indices, including runs across block boundaries."""
    begin = previous = None
    for start in range(0, len(values), block_rows):
        block = values[start : start + block_rows]
        if previous is None or int(block[0]) != previous + 1:
            if begin is not None:
                yield begin, previous
            begin = int(block[0])
        for end in np.flatnonzero(np.diff(block) != 1):
            yield begin, int(block[end])
            begin = int(block[end + 1])
        previous = int(block[-1])
    if begin is not None:
        yield begin, previous


class _MoleculePresence:
    """One molecule's effective presence, with bounded views for its consumers."""

    def __init__(self, values, frame_count, *, indexed):
        self.indexed = indexed
        values = np.asarray(values)
        if indexed:
            if values.ndim != 1 or (
                values.size and not np.issubdtype(values.dtype, np.integer)
            ):
                raise ValueError("Expected one-dimensional integer frame indices")
            values = values.astype(np.int64, copy=False)
            if values.size and values.min() < 0:
                # Match NumPy indexing used by idx_to_signal, without changing
                # the stored observed-frame record used by legacy CSV output.
                values = np.where(values < 0, values + frame_count, values)
            if values.size and (values.min() < 0 or values.max() >= frame_count):
                raise ValueError("Molecule frame index is outside analyzed frames")
            if np.any(values[1:] <= values[:-1]):
                values = np.unique(values)
        else:
            if values.ndim not in (1, 2) or values.size != len(values):
                raise ValueError("Expected one signal value per frame")
            if len(values) != frame_count:
                raise ValueError(
                    "Molecule signal length differs from analyzed frame count"
                )
            values = values.reshape(-1)
        self.values = values

    def ranges(self):
        """Yield maximal closed intervals, independent of CSV frame filters."""
        if self.indexed:
            yield from _index_ranges(self.values)
        else:
            yield from _signal_ranges(self.values)

    def selections(self):
        """Yield bounded frame indices or slices for completely occupied spans."""
        for start in range(0, len(self.values), _PRESENCE_BLOCK_ROWS):
            block = self.values[start : start + _PRESENCE_BLOCK_ROWS]
            if self.indexed:
                if int(block[-1]) - int(block[0]) == len(block) - 1:
                    yield slice(int(block[0]), int(block[-1]) + 1)
                else:
                    yield block
            elif np.all(block):
                yield slice(start, start + len(block))
            else:
                frames = np.flatnonzero(block)
                if frames.size:
                    yield frames + start


def _iter_molecule_presence(rng, *, structures=False):
    """Pair molecule records with effective presence in molecule-ID order.

    HMM records always consume a signal. Without HMM, sparse records reuse
    observed indices and only dense records consume the ordered origin stream.
    Files remain owned by the pipeline; both consumers can read them in turn.
    """
    fields = (0, 1, 2, 3) if structures else (0, 3)
    filename = rng.hmmfilename if rng.runHMM else rng.originfilename
    count = 0
    with (
        open(rng.moleculetemp2filename, "rb") as definitions,
        open(filename, "rb") as signals_file,
    ):
        signals = read_compressed_block(signals_file)
        for count, record in enumerate(_iter_molecule_records(definitions, fields), 1):
            if count > rng.hmmit:
                raise ValueError("Molecule count differs from declared count")
            indexed = not rng.runHMM and _use_frame_indices(record[-1], rng.step)
            block = record[-1] if indexed else next(signals, None)
            if block is None:
                raise ValueError("Missing molecule signal record")
            presence = _MoleculePresence(bytestolist(block), rng.step, indexed=indexed)
            yield record, presence
        if next(signals, None) is not None:
            raise ValueError("Extra molecule signal record")
    if count != rng.hmmit:
        raise ValueError("Molecule count differs from declared count")
