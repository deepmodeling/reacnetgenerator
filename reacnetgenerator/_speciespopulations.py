# SPDX-License-Identifier: LGPL-3.0-or-later
"""Aggregate observed species populations without a Counter for every frame.

Small/sparse inputs stay in a bounded lazy dictionary. Larger inputs use compact
counts or block-partitioned events. All paths preserve first-molecule order within
each frame, including duplicate observations. Effective HMM presence is separate:
this output deliberately reads the stored observed frame indices.
"""

import os
import tempfile
from collections import Counter
from contextlib import closing, suppress
from itertools import islice

import numpy as np

from ._step3state import _MoleculeNameTable, _unsigned_dtype_for_maximum
from .utils import WriteBuffer, bytestolist

_SMALL_COUNT_BYTES = 8 * 1024 * 1024
_SMALL_OBSERVATION_ROWS = 65536


_SPECIES_COUNT_BLOCK_BYTES = 16 * 1024 * 1024
_SPECIES_EVENT_CHUNK_ROWS = 65536
_SPECIES_OBSERVATION_SPOOL_BYTES = 256 * 1024 * 1024


def _frame_array(frames, frame_count):
    """Preserve duplicate visits and Python-list negative indexing semantics."""
    frames = np.asarray(frames)
    if frames.ndim != 1:
        raise TypeError("Expected one-dimensional molecule frame indices")
    if not frames.size:
        return np.zeros(0, dtype=np.int64)
    if frames.dtype.kind not in "iub":
        raise TypeError("Molecule frame indices must be integers")
    if int(frames.min()) < -frame_count or int(frames.max()) >= frame_count:
        raise IndexError("Molecule frame index is outside the trajectory")
    frames = frames.astype(np.int64, copy=False)
    if np.any(frames < 0):
        frames = np.where(frames < 0, frames + frame_count, frames)
    return frames


def _frame_records(filename, keys, frame_count, *, normalize=True):
    """Read observation fields in molecule order, skipping structure payloads."""
    with open(filename, "rb") as source:
        size = os.fstat(source.fileno()).st_size
        position = 0
        for key in keys:
            if position == size:
                return  # Match the legacy zip when fewer records than names exist.
            for field in range(4):
                header = source.read(64)
                if len(header) != 64:
                    raise ValueError("Incomplete molecule record header")
                length = int.from_bytes(header, "little")
                position += 64 + length
                if position > size:
                    raise ValueError("Incomplete molecule record payload")
                if field < 3:
                    source.seek(length, os.SEEK_CUR)
                else:
                    frames = bytestolist(header + source.read(length))
                    yield (
                        key,
                        _frame_array(frames, frame_count) if normalize else frames,
                    )


def _small_counts(records, frame_count):
    """Try a lazy dictionary with budgets for both live cells and work.

    Empty frames allocate nothing. The conservative per-cell/per-frame estimate
    selects a small-input optimization, not an exact Python allocator RSS bound.
    """
    counts, cells, observations = {}, 0, 0
    for name, frames in records:
        frames = np.asarray(frames)
        if frames.ndim != 1 or (frames.size and frames.dtype.kind not in "iub"):
            raise TypeError("Expected one-dimensional integer molecule frame indices")
        observations += len(frames)
        if observations > _SMALL_OBSERVATION_ROWS:
            return None
        for frame in frames.tolist():
            if frame < -frame_count or frame >= frame_count:
                raise IndexError("Molecule frame index is outside the trajectory")
            if frame < 0:
                frame += frame_count
            row = counts.get(frame)
            if row is None:
                row = counts[frame] = {}
            previous = row.get(name, 0)
            if not previous:
                cells += 1
                if cells * 320 + len(counts) * 192 > _SMALL_COUNT_BYTES:
                    return None
            row[name] = previous + 1
    return counts, observations, cells * 320 + len(counts) * 192


def _write_rows(filename, timesteps, rows, names=None):
    """Stream lines to retain the original per-frame formatting and order."""
    if names is not None:
        names = list(names)
    with WriteBuffer(open(filename, "w")) as output:
        for frame, counts in rows:
            output.append(f"Timestep {timesteps[frame]}:")
            # Batch formatting without materializing a species-wide line or
            # paying a buffer-limit check for every individual population.
            counts = iter(counts)
            while batch := tuple(islice(counts, 1200)):
                if names is None:
                    output.extend(f" {name} {count}" for name, count in batch)
                else:
                    output.extend(
                        f" {names[species]} {count}" for species, count in batch
                    )
            output.append("\n")


def _rescan_counts(filename, names, frame_count):
    """Bound count blocks by rereading the source when scratch storage fails."""
    width = len(names.names)
    block_frames = max(1, _SPECIES_COUNT_BLOCK_BYTES // max(1, width * 16))
    for start in range(0, frame_count, block_frames):
        stop = min(start + block_frames, frame_count)
        counts = np.zeros((width, stop - start), dtype=np.uint64)
        priority = np.zeros_like(counts)
        for molecule, (species, frames) in enumerate(
            _frame_records(filename, names.ids, frame_count)
        ):
            for offset in range(0, len(frames), _SPECIES_EVENT_CHUNK_ROWS):
                block = frames[offset : offset + _SPECIES_EVENT_CHUNK_ROWS]
                selected = block[(block >= start) & (block < stop)] - start
                np.add.at(counts[species], selected, 1)
                np.maximum.at(priority[species], selected, len(names) - molecule)
        for offset in range(stop - start):
            active = np.flatnonzero(counts[:, offset])
            active = active[np.argsort(priority[active, offset])[::-1]]
            yield (
                start + offset,
                ((int(species), int(counts[species, offset])) for species in active),
            )


def _write_species_populations(filename, names, timesteps, output):
    """Write observed populations; return resource diagnostics for benchmarks.

    The source, names, and timesteps remain caller-owned. Only scratch files
    created here are deleted. Count blocks need at least one species-wide row;
    input records and name/timestep metadata can independently exceed the budget.
    """
    frame_count = len(timesteps)
    with closing(
        _frame_records(filename, names, frame_count, normalize=False)
    ) as records:
        small = _small_counts(records, frame_count)
    if small is not None:
        counts, observations, working_bytes = small
        rows = ((frame, counts.get(frame, {}).items()) for frame in range(frame_count))
        _write_rows(output, timesteps, rows)
        return {
            "layout": "small",
            "observations": observations,
            "temporary_bytes": 0,
            "working_bytes": working_bytes,
        }

    names = (
        names
        if isinstance(names, _MoleculeNameTable)
        else _MoleculeNameTable.from_names(names)
    )
    store = None
    try:
        try:
            store = _SpeciesTimelineStore(frame_count, len(names.names), len(names))
            for molecule, (species, frames) in enumerate(
                _frame_records(filename, names.ids, frame_count)
            ):
                store.observe(molecule, species, frames)
            while store.refine_sparse_blocks():
                for molecule, (species, frames) in enumerate(
                    _frame_records(filename, names.ids, frame_count)
                ):
                    store.observe(molecule, species, frames)
            store.allocate()
            if store.requires_second_pass:
                for molecule, (species, frames) in enumerate(
                    _frame_records(filename, names.ids, frame_count)
                ):
                    store.add(molecule, species, frames)
        except OSError:
            if store is not None:
                with suppress(OSError):
                    store.close()
                store = None
            _write_rows(
                output,
                timesteps,
                _rescan_counts(filename, names, frame_count),
                names.names,
            )
            return {"layout": "rescan", "temporary_bytes": 0}
        _write_rows(output, timesteps, store.iter_frame_counts(), names.names)
        return {
            "layout": store.mode,
            "observations": store.total_rows,
            "temporary_bytes": store.peak_temporary_bytes,
            "working_bytes": store.maximum_count_block_bytes,
            "observation_reused": store.observation_reused,
        }
    finally:
        if store is not None:
            store.close()


class _SpeciesTimelineStore:
    """Bound aggregation buffers and choose compact dense or sparse scratch data."""

    def __init__(self, frame_count, species_count, molecule_count):
        """Initialize dimension-dependent dtypes, partition counts, and a spool.

        The first pass measures counts before choosing the final layout. An
        unavailable spool disables reuse without preventing source-file rescans.
        """
        self.frame_count = int(frame_count)
        self.species_count = int(species_count)
        self.molecule_count = int(molecule_count)
        if min(self.frame_count, self.species_count, self.molecule_count) < 0:
            raise ValueError("Species timeline dimensions must be non-negative")
        self.count_dtype = np.dtype(np.uint64)
        self.priority_dtype = _unsigned_dtype_for_maximum(self.molecule_count)
        self.force_dense = False
        self.maximum_count = 0
        self.cell_dtype = np.dtype(
            [
                ("count", self.count_dtype),
                ("priority", self.priority_dtype),
            ]
        )
        count_row_bytes = max(1, self.species_count * self.cell_dtype.itemsize)
        self.block_frame_count = max(
            1,
            _SPECIES_COUNT_BLOCK_BYTES // count_row_bytes,
        )
        self.block_frame_count = min(max(1, self.frame_count), self.block_frame_count)
        self.block_count = (
            self.frame_count + self.block_frame_count - 1
        ) // self.block_frame_count
        self.block_rows = np.zeros(self.block_count, dtype=np.uint64)
        self.total_rows = 0
        self.added_rows = 0
        self.mode = None
        self.path = None
        self.data: np.memmap | None = None
        self.block_offsets = None
        self.block_cursors = None
        self.frame_dtype = _unsigned_dtype_for_maximum(max(0, self.frame_count - 1))
        self.species_dtype = _unsigned_dtype_for_maximum(max(0, self.species_count - 1))
        self.event_dtype = np.dtype(
            [
                ("frame", self.frame_dtype),
                ("species", self.species_dtype),
            ]
        )
        self.observation_dtype = np.dtype(
            [
                ("frame", self.frame_dtype),
                ("species", self.species_dtype),
                ("priority", self.priority_dtype),
            ]
        )
        self.observation_path = None
        self.observation_file = None
        self.observation_buffer = None
        self.observation_buffer_rows = 0
        self.observation_rows = 0
        self.observation_complete = True
        self.observation_reused = False
        self.peak_temporary_bytes = 0
        self.dense_bytes = (
            self.frame_count * self.species_count * self.cell_dtype.itemsize
        )
        try:
            self._create_observation_spool()
        except OSError:
            with suppress(OSError):
                self._close_observation_spool(incomplete=True)

    def __enter__(self):
        """Return the store whose scratch files are owned by this context."""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Release owned scratch files on both successful and failed operations."""
        self.close()

    def _frame_array(self, frames):
        """Normalize observed indices against this store's trajectory length."""
        return _frame_array(frames, self.frame_count)

    @staticmethod
    def _chunks(values):
        """Limit vectorized work to the configured number of observation rows."""
        for start in range(0, len(values), _SPECIES_EVENT_CHUNK_ROWS):
            yield values[start : start + _SPECIES_EVENT_CHUNK_ROWS]

    def _create_observation_spool(self):
        """Own a temporary event file and fixed-size buffer for first-pass reuse."""
        handle, self.observation_path = tempfile.mkstemp(
            prefix="reacnetgenerator-species-observed-",
            suffix=".bin",
        )
        try:
            self.observation_file = os.fdopen(handle, "w+b")
        except BaseException:
            os.close(handle)
            raise
        self.observation_buffer = np.empty(
            _SPECIES_EVENT_CHUNK_ROWS,
            dtype=self.observation_dtype,
        )

    def _flush_observation_buffer(self):
        """Write buffered events before reusing their storage or reading the spool."""
        if self.observation_file is not None and self.observation_buffer_rows:
            assert self.observation_buffer is not None
            self.observation_file.write(
                self.observation_buffer[: self.observation_buffer_rows].tobytes()
            )
            self.observation_buffer_rows = 0

    def _close_observation_spool(self, *, incomplete=False):
        """Discard owned spool resources, optionally requiring a source reread.

        The buffer and path are released even when closing the file fails.
        """
        observation_file = self.observation_file
        self.observation_file = None
        try:
            if observation_file is not None:
                observation_file.close()
        finally:
            self.observation_buffer = None
            self.observation_buffer_rows = 0
            if self.observation_path is not None:
                with suppress(FileNotFoundError):
                    os.unlink(self.observation_path)
                self.observation_path = None
            if incomplete:
                self.observation_complete = False

    def _record_observation(self, molecule_index, species_id, frames):
        """Buffer events in molecule order, abandoning reuse if the spool cap is hit.

        Descending molecule priority preserves first appearance in dense counts.
        The cap includes events still waiting in the write buffer.
        """
        if self.observation_file is None:
            return
        projected_rows = self.observation_rows + len(frames)
        if (
            projected_rows * self.observation_dtype.itemsize
            > _SPECIES_OBSERVATION_SPOOL_BYTES
        ):
            self._close_observation_spool(incomplete=True)
            return
        assert self.observation_buffer is not None
        priority = self.molecule_count - molecule_index
        if len(frames) == 1:
            self.observation_buffer[self.observation_buffer_rows] = (
                int(frames[0]),
                species_id,
                priority,
            )
            self.observation_buffer_rows += 1
            if self.observation_buffer_rows == len(self.observation_buffer):
                self._flush_observation_buffer()
            self.observation_rows = projected_rows
            self._update_peak_temporary_bytes()
            return
        source_start = 0
        while source_start < len(frames):
            available = len(self.observation_buffer) - self.observation_buffer_rows
            row_count = min(available, len(frames) - source_start)
            source_stop = source_start + row_count
            target_stop = self.observation_buffer_rows + row_count
            target = self.observation_buffer[self.observation_buffer_rows : target_stop]
            target["frame"] = frames[source_start:source_stop]
            target["species"] = species_id
            target["priority"] = priority
            self.observation_buffer_rows = target_stop
            source_start = source_stop
            if self.observation_buffer_rows == len(self.observation_buffer):
                self._flush_observation_buffer()
        self.observation_rows = projected_rows
        self._update_peak_temporary_bytes()

    def _iter_observations(self):
        """Replay a complete spool in bounded chunks after flushing buffered events."""
        if not self.observation_complete or self.observation_file is None:
            raise RuntimeError("Species observation spool is incomplete")
        self._flush_observation_buffer()
        self.observation_file.flush()
        self.observation_file.seek(0)
        while True:
            rows = np.fromfile(
                self.observation_file,
                dtype=self.observation_dtype,
                count=_SPECIES_EVENT_CHUNK_ROWS,
            )
            if rows.size == 0:
                return
            yield rows

    def _create_temp_path(self):
        """Reserve an owned path for the final dense or sparse memory mapping."""
        handle, self.path = tempfile.mkstemp(
            prefix="reacnetgenerator-species-",
            suffix=".mmap",
        )
        os.close(handle)

    def _allocate_dense(self):
        """Map zero-filled count/priority cells and clean up if allocation fails."""
        self.mode = "dense"
        self._create_temp_path()
        assert self.path is not None
        try:
            self.data = np.memmap(
                self.path,
                mode="w+",
                dtype=self.cell_dtype,
                shape=(self.frame_count, self.species_count),
            )
            self._update_peak_temporary_bytes()
        except BaseException:
            self.close()
            raise

    def _update_peak_temporary_bytes(self):
        """Account for simultaneous final storage and buffered or spooled events."""
        final_bytes = os.path.getsize(self.path) if self.path is not None else 0
        observation_bytes = (
            self.observation_rows * self.observation_dtype.itemsize
            if self.observation_file is not None
            else 0
        )
        self.peak_temporary_bytes = max(
            self.peak_temporary_bytes,
            final_bytes + observation_bytes,
        )

    def _add_dense(self, molecule_index, species_id, frames):
        """Accumulate duplicate visits and retain the first contributing molecule."""
        assert self.data is not None
        priority = self.molecule_count - molecule_index
        for chunk in self._chunks(frames):
            np.add.at(self.data["count"][:, species_id], chunk, 1)
            np.maximum.at(
                self.data["priority"][:, species_id],
                chunk,
                priority,
            )

    def observe(self, molecule_index, species_id, frames):
        """Count frame-block rows before choosing the smaller disk layout."""
        if self.mode is not None:
            raise RuntimeError("Cannot observe species rows after allocation")
        species_id = int(species_id)
        if species_id < 0 or species_id >= self.species_count:
            raise IndexError("Species ID is outside the name table")
        molecule_index = int(molecule_index)
        if molecule_index < 0 or molecule_index >= self.molecule_count:
            raise IndexError("Molecule index is outside the name table")
        frames = self._frame_array(frames)
        # A sorted unique record contributes at most one count per frame.
        # Otherwise len(frames) is a safe upper bound, including repeats.
        previous = None
        unique = True
        for block in self._chunks(frames):
            if (previous is not None and int(block[0]) <= previous) or np.any(
                block[1:] <= block[:-1]
            ):
                unique = False
                break
            previous = int(block[-1])
        self.maximum_count += (1 if unique else len(frames)) if len(frames) else 0
        try:
            self._record_observation(molecule_index, species_id, frames)
        except OSError:
            with suppress(OSError):
                self._close_observation_spool(incomplete=True)
        if len(frames) == 1:
            block_id = int(frames[0]) // self.block_frame_count
            self.block_rows[block_id] += 1
            self.total_rows += 1
            return
        for chunk in self._chunks(frames):
            block_ids, counts = np.unique(
                chunk // self.block_frame_count,
                return_counts=True,
            )
            np.add.at(
                self.block_rows,
                block_ids.astype(np.intp, copy=False),
                counts.astype(np.uint64, copy=False),
            )
        self.total_rows += len(frames)

    def refine_sparse_blocks(self):
        """Bound sparse sorting with logarithmically many refinement passes.

        Return True when the caller must recount source records for the smaller
        partitions. A complete spool is recounted here instead. Return False
        once partitions fit or a single crowded frame requires dense counts.
        """
        self.count_dtype = _unsigned_dtype_for_maximum(self.maximum_count)
        self.cell_dtype = np.dtype(
            [("count", self.count_dtype), ("priority", self.priority_dtype)]
        )
        self.dense_bytes = (
            self.frame_count * self.species_count * self.cell_dtype.itemsize
        )
        if self.mode == "dense" or self.total_rows == 0:
            return False
        if self.dense_bytes <= self.total_rows * self.event_dtype.itemsize:
            return False
        bytes_per_sorted_row = self._sparse_working_bytes_per_row
        row_limit = max(1, _SPECIES_COUNT_BLOCK_BYTES // bytes_per_sorted_row)
        while True:
            maximum_rows = int(self.block_rows.max(initial=0))
            if maximum_rows <= row_limit:
                return False
            if self.block_frame_count == 1:
                # A single duplicate-heavy frame must not force an unbounded
                # event sort. Dense reduction scales with species, not events.
                self.force_dense = True
                return False
            # A near-limit crowded frame can keep maximum_rows almost
            # constant. Halving guarantees logarithmically many full passes,
            # even when the spool is unavailable and every pass rereads input.
            refined_frame_count = max(
                1,
                min(
                    self.block_frame_count // 2,
                    self.block_frame_count * row_limit // maximum_rows,
                ),
            )
            self.block_frame_count = refined_frame_count
            self.block_count = (
                self.frame_count + self.block_frame_count - 1
            ) // self.block_frame_count
            self.block_rows = np.zeros(self.block_count, dtype=np.uint64)
            if (
                not self.observation_complete
                or self.observation_rows != self.total_rows
            ):
                self.total_rows = 0
                self.maximum_count = 0
                return True
            for rows in self._iter_observations():
                block_ids, counts = np.unique(
                    rows["frame"] // self.block_frame_count,
                    return_counts=True,
                )
                np.add.at(
                    self.block_rows,
                    block_ids.astype(np.intp, copy=False),
                    counts.astype(np.uint64, copy=False),
                )

    def allocate(self):
        """Choose dense counts or sparse events from their exact disk sizes."""
        if self.mode == "dense":
            return
        if self.mode is not None:
            raise RuntimeError("Species timeline store is already allocated")
        if self.total_rows == 0 or self.frame_count == 0 or self.species_count == 0:
            self.mode = "empty"
            self._close_observation_spool()
            return
        if (
            self.dense_bytes <= self.total_rows * self.event_dtype.itemsize
            or self.force_dense
        ):
            self._allocate_dense()
            self._populate_from_observations()
            return
        self.mode = "sparse"
        self._create_temp_path()
        assert self.path is not None
        try:
            self.block_offsets = np.empty(self.block_count + 1, dtype=np.uint64)
            self.block_offsets[0] = 0
            np.cumsum(self.block_rows, out=self.block_offsets[1:])
            self.block_cursors = self.block_offsets[:-1].copy()
            self.data = np.memmap(
                self.path,
                mode="w+",
                dtype=self.event_dtype,
                shape=(self.total_rows,),
            )
            self._update_peak_temporary_bytes()
            self._populate_from_observations()
        except BaseException:
            self.close()
            raise

    def _populate_from_observations(self):
        """Populate the selected layout from a complete spool, then release it.

        An incomplete spool leaves the second source pass to the caller. Reuse
        is reported only after all spooled events have been added successfully.
        """
        if (
            not self.observation_complete
            or self.observation_file is None
            or self.observation_rows != self.total_rows
        ):
            return
        assert self.data is not None
        for rows in self._iter_observations():
            frame_ids = rows["frame"].astype(np.intp)
            species_ids = rows["species"].astype(np.intp)
            if self.mode == "dense":
                np.add.at(self.data["count"], (frame_ids, species_ids), 1)
                np.maximum.at(
                    self.data["priority"],
                    (frame_ids, species_ids),
                    rows["priority"],
                )
            else:
                self._populate_sparse_rows(rows, frame_ids)
            self.added_rows += len(rows)
        self.observation_reused = True
        self._close_observation_spool()

    def _populate_sparse_rows(self, rows, frame_ids):
        """Scatter one observation chunk by block without repeated full scans."""
        assert self.data is not None and self.block_cursors is not None
        row_count = len(rows)
        if row_count == 0:
            return
        block_ids = frame_ids // self.block_frame_count
        if row_count == 1 or np.all(block_ids[1:] >= block_ids[:-1]):
            order = None
            sorted_block_ids = block_ids
        else:
            order = np.argsort(block_ids, kind="stable")
            sorted_block_ids = block_ids[order]
        boundaries = np.flatnonzero(sorted_block_ids[1:] != sorted_block_ids[:-1]) + 1
        group_starts = np.concatenate((np.zeros(1, dtype=np.intp), boundaries))
        group_stops = np.concatenate(
            (boundaries, np.asarray([row_count], dtype=np.intp))
        )
        group_counts = group_stops - group_starts
        unique_blocks = sorted_block_ids[group_starts]
        cursor_starts = self.block_cursors[unique_blocks].astype(
            np.intp,
            copy=False,
        )
        target_positions = np.arange(row_count, dtype=np.intp)
        target_positions += np.repeat(
            cursor_starts - group_starts,
            group_counts,
        )
        if order is None:
            self.data["frame"][target_positions] = rows["frame"]
            self.data["species"][target_positions] = rows["species"]
        else:
            self.data["frame"][target_positions] = rows["frame"][order]
            self.data["species"][target_positions] = rows["species"][order]
        self.block_cursors[unique_blocks] += group_counts.astype(
            np.uint64,
            copy=False,
        )

    def add(self, molecule_index, species_id, frames):
        """Add one molecule timeline during the second sequential file pass."""
        if self.mode is None:
            raise RuntimeError("Species timeline store has not been allocated")
        species_id = int(species_id)
        if species_id < 0 or species_id >= self.species_count:
            raise IndexError("Species ID is outside the name table")
        molecule_index = int(molecule_index)
        if molecule_index < 0 or molecule_index >= self.molecule_count:
            raise IndexError("Molecule index is outside the name table")
        frames = self._frame_array(frames)
        if self.mode == "empty":
            if frames.size:
                raise RuntimeError("Observed and written species row counts differ")
            return
        if self.mode == "dense":
            self._add_dense(molecule_index, species_id, frames)
        else:
            assert self.data is not None and self.block_cursors is not None
            for chunk in self._chunks(frames):
                rows = np.empty(len(chunk), dtype=self.event_dtype)
                rows["frame"] = chunk
                rows["species"] = species_id
                self._populate_sparse_rows(rows, chunk)
        self.added_rows += len(frames)

    @property
    def requires_second_pass(self):
        """Report whether observed rows still need to be populated from the source."""
        return self.added_rows != self.total_rows

    def _check_complete(self):
        """Reject missing rows or incomplete sparse partitions before writing output."""
        if self.added_rows != self.total_rows:
            raise RuntimeError("Observed and written species row counts differ")
        if self.mode == "sparse":
            assert self.block_offsets is not None
            if not np.array_equal(self.block_cursors, self.block_offsets[1:]):
                raise RuntimeError("Species event block offsets are incomplete")

    def iter_frame_counts(self):
        """Yield ordered species counts one frame at a time with bounded memory."""
        self._check_complete()
        for block_index in range(self.block_count):
            frame_start = block_index * self.block_frame_count
            frame_stop = min(
                frame_start + self.block_frame_count,
                self.frame_count,
            )
            if self.mode == "dense":
                assert self.data is not None
                # Row-major storage exposes a view instead of copying a full
                # count block just to transpose it for output.
                count_block = np.asarray(self.data[frame_start:frame_stop])
                for frame_offset, row in enumerate(count_block):
                    active_species = np.flatnonzero(row["count"])
                    active_species = active_species[
                        np.argsort(row["priority"][active_species])[::-1]
                    ]
                    yield (
                        frame_start + frame_offset,
                        zip(
                            active_species.tolist(),
                            row["count"][active_species].tolist(),
                            strict=True,
                        ),
                    )
                continue
            if self.mode == "empty":
                for frame in range(frame_start, frame_stop):
                    yield frame, ()
                continue
            assert self.block_offsets is not None and self.data is not None
            row_start = int(self.block_offsets[block_index])
            row_stop = int(self.block_offsets[block_index + 1])
            rows = np.asarray(self.data[row_start:row_stop])
            order = np.argsort(rows["frame"], kind="stable")
            sorted_frames = rows["frame"][order]
            sorted_species = rows["species"][order]
            # Compute all boundaries once. A scalar searchsorted for every
            # frame can repeatedly promote the entire compact integer array.
            boundaries = np.flatnonzero(sorted_frames[1:] != sorted_frames[:-1]) + 1
            stops = iter(boundaries)
            event_start = 0
            for frame in range(frame_start, frame_stop):
                if (
                    event_start >= len(sorted_frames)
                    or int(sorted_frames[event_start]) != frame
                ):
                    yield frame, ()
                    continue
                event_stop = next(stops, len(sorted_frames))
                yield (
                    frame,
                    self._iter_sparse_counts(sorted_species[event_start:event_stop]),
                )
                event_start = event_stop

    @staticmethod
    def _iter_sparse_counts(species_ids):
        """Aggregate one frame without Python objects per active species."""
        if len(species_ids) <= 128:
            yield from Counter(species_ids.tolist()).items()
            return
        unique_species, first_indices, counts = np.unique(
            species_ids,
            return_index=True,
            return_counts=True,
        )
        for index in np.argsort(first_indices):
            yield int(unique_species[index]), int(counts[index])

    @property
    def _sparse_working_bytes_per_row(self):
        """Conservatively estimate stable-sort and unique scratch per event."""
        return (
            10 * np.dtype(np.intp).itemsize
            + 4 * self.event_dtype.itemsize
            + np.dtype(np.bool_).itemsize
        )

    @property
    def temporary_bytes(self):
        """Return the current final-layout file size, excluding the observation spool."""
        if self.path is None:
            return 0
        return os.path.getsize(self.path)

    @property
    def maximum_count_block_bytes(self):
        """Estimate count or sorting workspace, excluding metadata and mapped pages."""
        if self.mode == "sparse":
            return int(self.block_rows.max(initial=0)) * (
                self._sparse_working_bytes_per_row
            )
        return min(self.frame_count, self.block_frame_count) * (
            self.species_count * self.cell_dtype.itemsize
        )

    def close(self):
        """Release mappings and unlink temporary files even after failed writes."""
        try:
            self._close_observation_spool()
        finally:
            try:
                if self.data is not None:
                    mapping = getattr(self.data, "_mmap", None)
                    if mapping is not None:
                        mapping.close()
                    self.data = None
            finally:
                if self.path is not None:
                    with suppress(FileNotFoundError):
                        os.unlink(self.path)
                    self.path = None
