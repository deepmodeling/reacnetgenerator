# SPDX-License-Identifier: LGPL-3.0-or-later
"""Test compact temporary state used by Step 3."""

import os
import tempfile

import numpy as np
import pytest

from reacnetgenerator._step3state import (
    _AtomFrameStore,
    _MoleculeNameBuilder,
    _unsigned_dtype_for_maximum,
)


@pytest.mark.parametrize(
    ("maximum", "expected"),
    [
        (0, np.uint8),
        (255, np.uint8),
        (256, np.uint16),
        (65535, np.uint16),
        (65536, np.uint32),
        (2**32, np.uint64),
    ],
)
def test_unsigned_dtype_uses_narrowest_width(maximum, expected):
    """Molecule IDs should use the narrowest lossless unsigned dtype."""
    assert _unsigned_dtype_for_maximum(maximum) == np.dtype(expected)


def test_unsigned_dtype_rejects_values_above_uint64():
    """Unsupported molecule ID ranges should fail instead of wrapping."""
    with pytest.raises(OverflowError, match="unsigned 64-bit"):
        _unsigned_dtype_for_maximum(int(np.iinfo(np.uint64).max) + 1)


def test_molecule_name_builder_deduplicates_species():
    """Per-molecule IDs should index one copy of each species name."""
    builder = _MoleculeNameBuilder(4)
    for name in ("A", "B", "A", "C"):
        builder.append(name)

    table = builder.finish()

    assert table.ids.dtype == np.dtype(np.uint8)
    assert table.ids.tolist() == [0, 1, 0, 2]
    assert table.names.tolist() == ["A", "B", "C"]
    assert list(table) == ["A", "B", "A", "C"]
    np.testing.assert_array_equal(table[np.array([3, 0, 2])], ["C", "A", "A"])
    np.testing.assert_array_equal(np.asarray(table), ["A", "B", "A", "C"])


def test_molecule_name_builder_checks_declared_count():
    """Name-table construction should detect incomplete or excess records."""
    incomplete = _MoleculeNameBuilder(1)
    with pytest.raises(RuntimeError, match="Fewer molecule names"):
        incomplete.finish()

    excess = _MoleculeNameBuilder(1)
    excess.append("A")
    with pytest.raises(RuntimeError, match="More molecule names"):
        excess.append("B")


def test_atom_frame_store_is_compact_transposable_and_recoverable(tmp_path):
    """Disk-backed matrices should retain ndarray semantics and clean up."""
    with _AtomFrameStore((3, 10), 256, directory=tmp_path) as store:
        atomeach_path = store.atomeach_path
        conflict_path = store.conflict.path
        assert store.atomeach.dtype == np.dtype(np.uint16)
        assert store.atomeach.nbytes == 60
        assert store.conflict.nbytes == 6
        assert store.nbytes == 66

        store.atomeach[0, [0, 8]] = 1
        store.atomeach[1, [8, 9]] = 256
        store.conflict.mark(
            np.array([0, 1]),
            np.array([0, 8, 9]),
            np.array(
                [
                    [True, True, False],
                    [False, True, True],
                ]
            ),
        )
        store.flush()

        expected = np.zeros((3, 10), dtype=np.bool_)
        expected[0, [0, 8]] = True
        expected[1, [8, 9]] = True
        np.testing.assert_array_equal(store.conflict, expected)
        np.testing.assert_array_equal(store.conflict.T, expected.T)
        np.testing.assert_array_equal(
            tuple(store.conflict.T[8:10]),
            tuple(expected.T[8:10]),
        )

    assert not os.path.exists(atomeach_path)
    assert not os.path.exists(conflict_path)


def test_packed_conflict_marks_sparse_selection(tmp_path):
    """Sparse overlap masks should mark only their selected cells."""
    with _AtomFrameStore((3, 10), 4, directory=tmp_path) as store:
        overlap = np.zeros((3, 10), dtype=np.bool_)
        overlap[2, 9] = True
        store.conflict.mark(np.arange(3), np.arange(10), overlap)

        expected = np.zeros((3, 10), dtype=np.bool_)
        expected[2, 9] = True
        np.testing.assert_array_equal(store.conflict, expected)


def test_atom_frame_store_cleans_up_after_second_allocation_failure(
    tmp_path, monkeypatch
):
    """A failed conflict allocation should remove the atom mapping file."""
    original_mkstemp = tempfile.mkstemp
    created_paths = []
    allocation_count = 0

    def fail_second_mkstemp(*args, **kwargs):
        nonlocal allocation_count
        allocation_count += 1
        if allocation_count == 2:
            raise OSError("simulated conflict allocation failure")
        handle, path = original_mkstemp(*args, **kwargs)
        created_paths.append(path)
        return handle, path

    monkeypatch.setattr(
        "reacnetgenerator._step3state.tempfile.mkstemp",
        fail_second_mkstemp,
    )

    with pytest.raises(OSError, match="simulated conflict allocation failure"):
        _AtomFrameStore((3, 10), 4, directory=tmp_path)

    assert allocation_count == 2
    assert len(created_paths) == 1
    assert not os.path.exists(created_paths[0])
