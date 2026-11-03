# SPDX-License-Identifier: LGPL-3.0-or-later
"""Print reaction counts selected from a saved timeline."""

import argparse

from .timedoutput import iter_reaction_types


def main(argv=None):
    """Write a .reactionabcd-style summary to standard output."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("timeline", help="ReacNetGenerator HDF5 timeline")
    parser.add_argument(
        "--inferred-order",
        choices=("include", "exclude", "only"),
        default="include",
        help="how to count coordinate-input events with only bond-order changes",
    )
    args = parser.parse_args(argv)
    try:
        counts = iter_reaction_types(args.timeline, inferred_order=args.inferred_order)
        for _, reactant, product, count in sorted(counts, key=lambda row: -row[3]):
            print(f"{count} {reactant}->{product}")
    except (OSError, ValueError) as exc:
        parser.exit(1, f"{parser.prog}: {exc}\n")


if __name__ == "__main__":
    main()
