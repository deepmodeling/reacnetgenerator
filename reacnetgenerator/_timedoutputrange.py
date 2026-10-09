# SPDX-License-Identifier: LGPL-3.0-or-later
"""Random-access presence checks for molecule range tables."""


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
