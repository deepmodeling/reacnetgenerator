# SPDX-License-Identifier: LGPL-3.0-or-later
"""Configuration shared by bond detection and timeline validation."""

import sys
from numbers import Real


def normalize_hysteresis(config, inputfiletype):
    """Check the two thresholds in the units and direction of the input."""
    if not isinstance(config, dict) or set(config) != {"form", "break"}:
        raise ValueError("bond_hysteresis requires form and break thresholds")
    for value in config.values():
        if (
            isinstance(value, bool)
            or not isinstance(value, Real)
            or not 0 <= value <= sys.float_info.max
        ):
            raise ValueError(
                "bond_hysteresis thresholds must be finite nonnegative numbers"
            )
    form, break_threshold = config["form"], config["break"]
    if inputfiletype in ("bond", "lammpsbondfile"):
        if form <= break_threshold:
            raise ValueError("bond_hysteresis for bond orders requires form > break")
    elif inputfiletype in ("dump", "lammpsdumpfile", "xyz", "extxyz"):
        if not 0 < form < break_threshold:
            raise ValueError("bond_hysteresis for distances requires 0 < form < break")
    else:
        raise ValueError("bond_hysteresis requires a supported inputfiletype")
    return {"form": float(form), "break": float(break_threshold)}


def molecule_range_basis(run_hmm, hysteresis):
    """Name the last operation defining molecule occupancy in this timeline."""
    if run_hmm:
        return "HMM signal"
    return "hysteresis signal" if hysteresis is not None else "observed signal"
