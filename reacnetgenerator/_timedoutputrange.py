# SPDX-License-Identifier: LGPL-3.0-or-later
"""Random-access presence checks for molecule range tables."""

from collections import OrderedDict

_RANGE_BLOCK_ROWS = 8192
_RANGE_CACHE_BLOCKS = 16


class _RangeColumn:
    """Read scalar integers through a shared, bounded cache of array blocks."""

    def __init__(self, dataset, column, cache):
        self.dataset = dataset
        self.column = column
        self.cache = cache
        self.length = len(dataset)

    def __len__(self):
        return self.length

    def __getitem__(self, index):
        if not 0 <= index < self.length:
            raise IndexError("Range column index is out of bounds")
        block, offset = divmod(index, _RANGE_BLOCK_ROWS)
        key = self.column, block
        if key not in self.cache:
            # Evict before reading so a miss never retains an extra full block.
            if len(self.cache) >= _RANGE_CACHE_BLOCKS:
                self.cache.popitem(last=False)
            start = block * _RANGE_BLOCK_ROWS
            self.cache[key] = self.dataset[
                start : min(start + _RANGE_BLOCK_ROWS, self.length)
            ]
        self.cache.move_to_end(key)
        return self.cache[key][offset]


class _MoleculeRangeLookup:
    """Check an immutable range table with at most 1 MiB of cached int64 data.

    All three columns share 16 blocks of at most 8192 rows. Small array and
    dictionary headers are additional; this is not a total process RSS bound.
    Create the lookup only after all ranges are flushed, and discard it before
    the backing file closes. The original scalar search remains the oracle.
    """

    def __init__(self, molecule_ids, starts, ends):
        cache = OrderedDict()
        self.columns = tuple(
            _RangeColumn(dataset, column, cache)
            for column, dataset in enumerate((molecule_ids, starts, ends))
        )

    def contains(self, molecule_id, frame):
        """Retain exact inclusive endpoints, gaps and missing-ID behavior."""
        return molecule_present(molecule_id, frame, *self.columns)


def _lower_bound(dataset, value, start=0, stop=None):
    """Return the first position whose integer value is at least ``value``."""
    right = len(dataset) if stop is None else stop
    left = start
    while left < right:
        middle = (left + right) // 2
        if int(dataset[middle]) < value:
            left = middle + 1
        else:
            right = middle
    return left


def molecule_present(molecule_id, frame, molecule_ids, starts, ends):
    """Return whether one molecule's maximal ranges contain ``frame``."""
    first = _lower_bound(molecule_ids, molecule_id)
    stop = _lower_bound(molecule_ids, molecule_id + 1, first)
    candidate = _lower_bound(starts, frame + 1, first, stop) - 1
    return candidate >= first and int(ends[candidate]) >= frame
