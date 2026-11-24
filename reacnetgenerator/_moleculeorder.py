# SPDX-License-Identifier: LGPL-3.0-or-later
"""Disk-backed structural ordering of aligned HMM results."""

import pickle
import sqlite3
import tempfile
from contextlib import ExitStack, closing
from pathlib import Path

import numpy as np

from .utils import bytestolist

_ORDER_CACHE_KIB = 8192


def _structure_key(record):
    """Encode detection's sorted, nonnegative atom/bond/order tuples numerically.

    Detection deduplicates these three fields before filtering. Big-endian
    integers preserve their numeric lexicographic order in SQLite BLOB keys,
    including proper-prefix tuples, without depending on pickle/compression.
    Only one molecule's structure is decoded; observed frames and full-length
    signals are never part of the sort key.
    """
    return tuple(
        np.asarray(bytestolist(field), dtype=">i8").tobytes() for field in record[:3]
    )


class _OrderedMoleculeResults:
    """Spool results and stream them in structural order with a bounded cache.

    The primary-key B-tree is built incrementally on disk, so iteration needs
    neither an in-memory list of keys nor an external SQL sort. Payloads live
    in a separate file, avoiding SQLite's row-size limit for long signals.
    Memory still scales with one indivisible molecule result, as in the HMM
    runner. All scratch files are removed on success and on failure.
    """

    def __enter__(self):
        self._stack = ExitStack()
        try:
            directory = Path(
                self._stack.enter_context(
                    tempfile.TemporaryDirectory(prefix="reacnetgenerator-order-")
                )
            )
            self._payload = self._stack.enter_context(
                (directory / "results.bin").open("w+b")
            )
            self._index = self._stack.enter_context(
                closing(sqlite3.connect(directory / "index.sqlite"))
            )
            # This disposable index has no durability requirement. Its B-tree
            # pages spill to disk instead of retaining every structural key.
            self._index.execute(f"PRAGMA cache_size = -{_ORDER_CACHE_KIB}")
            self._index.execute("PRAGMA journal_mode = OFF")
            self._index.execute("PRAGMA synchronous = OFF")
            self._index.execute(
                "CREATE TABLE results ("
                "atoms BLOB, bonds BLOB, orders BLOB, offset INTEGER, "
                "PRIMARY KEY (atoms, bonds, orders)) WITHOUT ROWID"
            )
        except BaseException:
            self._stack.close()
            raise
        return self

    def add(self, origin, filtered, record):
        """Retain one complete result without separating its aligned fields."""
        key = _structure_key(record)
        offset = self._payload.tell()
        pickle.dump(
            (origin, filtered, record), self._payload, protocol=pickle.HIGHEST_PROTOCOL
        )
        self._index.execute("INSERT INTO results VALUES (?, ?, ?, ?)", (*key, offset))

    def atom_set_count(self):
        """Count adjacent atom sets for ordered, bounded HMM worker dispatch."""
        return self._index.execute(
            "SELECT COUNT(DISTINCT atoms) FROM results"
        ).fetchone()[0]

    def __iter__(self):
        """Yield at most one compressed molecule and its signals at a time."""
        self._index.commit()
        self._payload.flush()
        with closing(
            self._index.execute(
                "SELECT offset FROM results ORDER BY atoms, bonds, orders"
            )
        ) as rows:
            for (offset,) in rows:
                self._payload.seek(offset)
                yield pickle.load(self._payload)

    def __exit__(self, *exc):
        return self._stack.__exit__(*exc)
