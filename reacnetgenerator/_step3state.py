# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compact temporary state used while collecting reaction paths."""

from __future__ import annotations

import operator
import os
import tempfile
from typing import Literal

import numpy as np

_PACKED_BOOL_SPARSE_SELECTION_DIVISOR = 8


def _unsigned_dtype_for_maximum(maximum_value: int) -> np.dtype:
    """Return the narrowest unsigned dtype containing ``maximum_value``."""
    maximum = max(0, int(maximum_value))
    for dtype in (np.uint8, np.uint16, np.uint32, np.uint64):
        if maximum <= np.iinfo(dtype).max:
            return np.dtype(dtype)
    raise OverflowError("Value exceeds unsigned 64-bit integer storage")


class _MoleculeNameTable:
    """Map each molecule ID to one deduplicated species name."""

    def __init__(self, ids, names, *, validate: bool = True) -> None:
        self.ids = np.asarray(ids)
        self.names = np.asarray(names)
        if self.ids.ndim != 1 or self.names.ndim != 1:
            raise ValueError("Molecule name IDs and values must be one-dimensional")
        if not np.issubdtype(self.ids.dtype, np.unsignedinteger):
            raise TypeError("Molecule name IDs must use an unsigned integer dtype")
        if validate and self.ids.size and int(self.ids.max()) >= len(self.names):
            raise ValueError("Molecule name ID is outside the species-name table")

    @classmethod
    def from_names(cls, values):
        """Build a compact table from an iterable of molecule names."""
        values = list(values)
        builder = _MoleculeNameBuilder(len(values))
        for value in values:
            builder.append(value)
        return builder.finish()

    def __len__(self) -> int:
        return len(self.ids)

    def __iter__(self):
        return (self.names[int(name_id)] for name_id in self.ids)

    def __getitem__(self, index):
        return self.names[self.ids[index]]

    def __array__(self, dtype=None, copy=None):
        values = self.names[self.ids]
        if dtype is not None:
            values = values.astype(dtype, copy=False)
        if copy:
            values = values.copy()
        return values


class _MoleculeNameBuilder:
    """Build a fixed-size name table without retaining duplicate strings."""

    def __init__(self, molecule_count: int) -> None:
        molecule_count = int(molecule_count)
        if molecule_count < 0:
            raise ValueError("Molecule count must be non-negative")
        self.ids = np.empty(
            molecule_count,
            dtype=_unsigned_dtype_for_maximum(max(0, molecule_count - 1)),
        )
        self.names: list[str] = []
        self._name_ids: dict[str, int] = {}
        self._count = 0

    def append(self, name: str) -> None:
        """Append one molecule name in molecule-ID order."""
        if self._count >= len(self.ids):
            raise RuntimeError("More molecule names than the declared count")
        name = str(name)
        name_id = self._name_ids.get(name)
        if name_id is None:
            name_id = len(self.names)
            self._name_ids[name] = name_id
            self.names.append(name)
        self.ids[self._count] = name_id
        self._count += 1

    def finish(self) -> _MoleculeNameTable:
        """Return the completed compact table."""
        if self._count != len(self.ids):
            raise RuntimeError("Fewer molecule names than the declared count")
        names = np.asarray(self.names, dtype=str)
        target_dtype = _unsigned_dtype_for_maximum(max(0, len(names) - 1))
        ids = (
            self.ids
            if self.ids.dtype == target_dtype
            else self.ids.astype(target_dtype)
        )
        return _MoleculeNameTable(ids, names, validate=False)


def _is_full_slice(index) -> bool:
    return (
        isinstance(index, slice)
        and index.start is None
        and index.stop is None
        and index.step is None
    )


class _PackedBoolFrameView:
    """Expose packed atom-by-frame data as a lazy frame-by-atom view."""

    def __init__(self, matrix: _PackedBoolMatrix, frames=None) -> None:
        self._matrix = matrix
        self._frames = range(matrix.shape[1]) if frames is None else frames
        self.shape = (len(self._frames), matrix.shape[0])
        self.dtype = np.dtype(np.bool_)

    def __len__(self) -> int:
        return len(self._frames)

    def __iter__(self):
        for frame in self._frames:
            yield self._matrix._read_column(frame, slice(None))

    def __getitem__(self, index):
        if isinstance(index, slice):
            return _PackedBoolFrameView(self._matrix, self._frames[index])
        frame = self._frames[operator.index(index)]
        return self._matrix._read_column(frame, slice(None))

    def __array__(self, dtype=None, copy=None):
        if len(self) == 0:
            values = np.empty(self.shape, dtype=np.bool_)
        else:
            values = np.stack(tuple(self), axis=0)
        if dtype is not None:
            values = values.astype(dtype, copy=False)
        if copy:
            values = values.copy()
        return values


class _PackedBoolMatrix:
    """Store a logical row-major boolean matrix using one bit per cell."""

    def __init__(
        self,
        path: str,
        shape,
        *,
        mode: Literal["r", "r+", "w+", "c"] = "w+",
    ) -> None:
        self.path = path
        self.shape = tuple(int(value) for value in shape)
        if len(self.shape) != 2 or min(self.shape) <= 0:
            raise ValueError("Packed boolean matrix shape must be two positive values")
        self.dtype = np.dtype(np.bool_)
        self.packed_shape = (self.shape[0], (self.shape[1] + 7) // 8)
        self.data = np.memmap(
            path,
            mode=mode,
            dtype=np.uint8,
            shape=self.packed_shape,
        )
        if mode.startswith("w"):
            self.data.fill(0)
        self._closed = False

    @property
    def T(self) -> _PackedBoolFrameView:
        """Return a lazy frame-by-atom transpose view."""
        return _PackedBoolFrameView(self)

    @property
    def nbytes(self) -> int:
        """Return the size of the packed payload in bytes."""
        return int(np.prod(self.packed_shape, dtype=np.int64))

    def __len__(self) -> int:
        return self.shape[0]

    def _read_column(self, column, rows):
        column = operator.index(column)
        if column < 0:
            column += self.shape[1]
        if column < 0 or column >= self.shape[1]:
            raise IndexError("Packed boolean matrix column is out of range")
        values = self.data[rows, column // 8]
        return np.not_equal(
            np.bitwise_and(values, np.uint8(1 << (column % 8))),
            0,
        )

    def _unpack_rows(self, rows):
        packed = np.asarray(self.data[rows])
        values = np.unpackbits(packed, axis=-1, bitorder="little")
        return values[..., : self.shape[1]].astype(np.bool_, copy=False)

    def __getitem__(self, index):
        if isinstance(index, tuple):
            if len(index) != 2:
                raise IndexError("Packed boolean matrix requires two indices")
            rows, columns = index
            try:
                column = operator.index(columns)
            except TypeError:
                values = self._unpack_rows(rows)
                return values[..., columns]
            return self._read_column(column, rows)
        return self._unpack_rows(index)

    def __setitem__(self, index, value) -> None:
        if not _is_full_slice(index):
            raise TypeError("Packed boolean matrix only supports whole-matrix fill")
        if not np.isscalar(value):
            raise TypeError("Packed boolean matrix fill value must be a scalar")
        self.data.fill(255 if bool(value) else 0)
        if value and self.shape[1] % 8:
            self.data[:, -1] &= np.uint8((1 << (self.shape[1] % 8)) - 1)

    def __array__(self, dtype=None, copy=None):
        values = self._unpack_rows(slice(None))
        if dtype is not None:
            values = values.astype(dtype, copy=False)
        if copy:
            values = values.copy()
        return values

    def mark(self, atoms, frames, overlap) -> None:
        """Set cells selected by an atom-by-frame overlap mask."""
        atoms = np.asarray(atoms, dtype=np.int64).reshape((-1,))
        frames = np.asarray(frames, dtype=np.int64).reshape((-1,))
        overlap = np.asarray(overlap, dtype=np.bool_)
        if overlap.shape != (len(atoms), len(frames)):
            raise ValueError("Conflict overlap shape does not match selected cells")
        if len(atoms) == 0 or len(frames) == 0:
            return
        if np.any((frames < 0) | (frames >= self.shape[1])):
            raise IndexError("Packed boolean matrix frame is out of range")
        active_count = int(np.count_nonzero(overlap))
        if active_count == 0:
            return
        if active_count * _PACKED_BOOL_SPARSE_SELECTION_DIVISOR <= overlap.size:
            active_positions = np.flatnonzero(overlap)
            atom_positions, frame_positions = np.divmod(
                active_positions,
                len(frames),
            )
            selected_frames = frames[frame_positions]
            np.bitwise_or.at(
                self.data,
                (atoms[atom_positions], selected_frames // 8),
                np.asarray(1 << (selected_frames % 8), dtype=np.uint8),
            )
            return
        bit_masks = np.asarray(1 << (frames % 8), dtype=np.uint8)
        write_masks = overlap.view(np.uint8).copy()
        np.multiply(write_masks, bit_masks, out=write_masks)
        np.bitwise_or.at(
            self.data,
            (atoms[:, np.newaxis], (frames // 8)[np.newaxis, :]),
            write_masks,
        )

    def flush(self) -> None:
        """Flush the packed mapping to its temporary file."""
        if not self._closed:
            self.data.flush()

    def close(self) -> None:
        """Close the mapping and remove its temporary file."""
        if self._closed:
            return
        self._closed = True
        try:
            try:
                self.data.flush()
            finally:
                mmap = getattr(self.data, "_mmap", None)
                if mmap is not None:
                    mmap.close()
        finally:
            try:
                os.unlink(self.path)
            except FileNotFoundError:
                pass


class _AtomFrameStore:
    """Own compact disk-backed atom-by-frame matrices for Step 3."""

    def __init__(self, shape, maximum_molecule_id: int, *, directory=None) -> None:
        self.shape = tuple(int(value) for value in shape)
        if len(self.shape) != 2 or min(self.shape) <= 0:
            raise ValueError("Atom-frame matrix shape must be two positive values")
        self.molecule_dtype = _unsigned_dtype_for_maximum(maximum_molecule_id)
        atom_handle, self.atomeach_path = tempfile.mkstemp(
            prefix="reacnetgenerator-atomeach-",
            suffix=".mmap",
            dir=directory,
        )
        os.close(atom_handle)
        conflict_handle, conflict_path = tempfile.mkstemp(
            prefix="reacnetgenerator-conflict-",
            suffix=".mmap",
            dir=directory,
        )
        os.close(conflict_handle)
        try:
            self.atomeach = np.memmap(
                self.atomeach_path,
                mode="w+",
                dtype=self.molecule_dtype,
                shape=self.shape,
            )
            self.atomeach.fill(0)
            self.conflict = _PackedBoolMatrix(conflict_path, self.shape)
            self._closed = False
        except BaseException:
            atomeach = getattr(self, "atomeach", None)
            if atomeach is not None:
                mmap = getattr(atomeach, "_mmap", None)
                if mmap is not None:
                    mmap.close()
            conflict = getattr(self, "conflict", None)
            if conflict is not None:
                conflict.close()
            for path in (self.atomeach_path, conflict_path):
                try:
                    os.unlink(path)
                except FileNotFoundError:
                    pass
            self._closed = True
            raise

    @property
    def nbytes(self) -> int:
        """Return the total disk-backed payload size in bytes."""
        return int(self.atomeach.nbytes) + self.conflict.nbytes

    def flush(self) -> None:
        """Flush both mappings before they are consumed."""
        self.atomeach.flush()
        self.conflict.flush()

    def close(self) -> None:
        """Close both mappings and remove their temporary files."""
        if getattr(self, "_closed", False):
            return
        self._closed = True
        try:
            try:
                self.atomeach.flush()
            finally:
                mmap = getattr(self.atomeach, "_mmap", None)
                if mmap is not None:
                    mmap.close()
        finally:
            try:
                os.unlink(self.atomeach_path)
            except FileNotFoundError:
                pass
            finally:
                self.conflict.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()

    def __del__(self):
        try:
            self.close()
        except BaseException:
            pass
