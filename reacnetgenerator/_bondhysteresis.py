# SPDX-License-Identifier: LGPL-3.0-or-later
"""Configuration shared by bond detection and timeline validation."""

import math
from numbers import Real


def normalize_hysteresis(config, inputfiletype):
    """Check the two thresholds in the units and direction of the input."""
    if not isinstance(config, dict) or set(config) != {"form", "break"}:
        raise ValueError("bond_hysteresis requires form and break thresholds")
    parsed = {}
    for key in ("form", "break"):
        value = config[key]
        if isinstance(value, bool) or not isinstance(value, Real):
            raise ValueError(
                "bond_hysteresis thresholds must be finite nonnegative numbers"
            )
        try:
            numeric = float(value)
        except OverflowError as exc:
            raise ValueError(
                "bond_hysteresis thresholds must be finite nonnegative numbers"
            ) from exc
        # Validate after float conversion so underflow or rounding collisions are rejected.
        if not math.isfinite(numeric) or numeric < 0:
            raise ValueError(
                "bond_hysteresis thresholds must be finite nonnegative numbers"
            )
        parsed[key] = numeric
    form, break_threshold = parsed["form"], parsed["break"]
    if inputfiletype in ("bond", "lammpsbondfile"):
        if form <= break_threshold:
            raise ValueError("bond_hysteresis for bond orders requires form > break")
    elif inputfiletype in ("dump", "lammpsdumpfile", "xyz", "extxyz"):
        if not 0 < form < break_threshold:
            raise ValueError("bond_hysteresis for distances requires 0 < form < break")
    else:
        raise ValueError("bond_hysteresis requires a supported inputfiletype")
    return parsed


def molecule_range_basis(run_hmm, hysteresis):
    """Name the last operation defining molecule occupancy in this timeline."""
    if run_hmm:
        return "HMM signal"
    return "hysteresis signal" if hysteresis is not None else "observed signal"
