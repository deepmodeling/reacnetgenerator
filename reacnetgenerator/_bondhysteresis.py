# SPDX-License-Identifier: LGPL-3.0-or-later
"""Configuration shared by bond detection and timeline validation."""

import math
from numbers import Real
from typing import TypedDict

# The optional map is defined separately because 'break' is a Python keyword.
_BaseThresholds = TypedDict("_BaseThresholds", {"form": float, "break": float})


class BondHysteresis(_BaseThresholds, total=False):
    """Normalized bond detection parameters recorded in the timeline."""

    pair_break: dict[str, float]


def _threshold(value) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(
            "bond_hysteresis thresholds must be finite nonnegative numbers"
        )
    try:
        numeric = float(value)
    except OverflowError as exc:
        raise ValueError("bond_hysteresis threshold is out of range") from exc
    if not math.isfinite(numeric) or numeric < 0:
        raise ValueError(
            "bond_hysteresis thresholds must be finite nonnegative numbers"
        )
    return numeric


def normalize_hysteresis(config, inputfiletype) -> BondHysteresis:
    """Validate units, finite float boundaries and unordered element pairs."""
    if (
        not isinstance(config, dict)
        or config.keys() - {"form", "break", "pair_break"}
        or not {"form", "break"} <= config.keys()
    ):
        raise ValueError("bond_hysteresis requires form and break thresholds")
    parsed = BondHysteresis(
        {
            "form": _threshold(config["form"]),
            "break": _threshold(config["break"]),
        }
    )
    form, break_threshold = parsed["form"], parsed["break"]
    if inputfiletype in ("bond", "lammpsbondfile"):
        if form <= break_threshold:
            raise ValueError("bond_hysteresis for bond orders requires form > break")
    elif inputfiletype in ("dump", "lammpsdumpfile", "xyz", "extxyz"):
        if not 0 < form < break_threshold:
            raise ValueError("bond_hysteresis for distances requires 0 < form < break")
    else:
        raise ValueError("bond_hysteresis requires a supported inputfiletype")
    if "pair_break" in config:
        from ase.data import atomic_numbers

        overrides = config["pair_break"]
        if (
            inputfiletype in ("bond", "lammpsbondfile")
            or not isinstance(overrides, dict)
            or not overrides
        ):
            raise ValueError(
                "bond_hysteresis pair_break requires a nonempty mapping for coordinate inputs"
            )
        pairs = {}
        for key, value in overrides.items():
            parts = key.split("-") if isinstance(key, str) else []
            if len(parts) != 2 or not all(atomic_numbers.get(el) for el in parts):
                raise ValueError(
                    "bond_hysteresis pair_break requires element pairs such as H-H"
                )
            pair = "-".join(sorted(parts))
            if pair in pairs:
                raise ValueError(f"duplicate bond_hysteresis pair_break for {pair}")
            limit = _threshold(value)
            if limit <= form:
                raise ValueError(
                    "bond_hysteresis pair_break must exceed the form multiplier"
                )
            pairs[pair] = limit
        parsed["pair_break"] = dict(sorted(pairs.items()))
    return parsed


def molecule_range_basis(run_hmm, hysteresis):
    """Name the last operation defining molecule occupancy in this timeline."""
    if run_hmm:
        return "HMM signal"
    return "hysteresis signal" if hysteresis is not None else "observed signal"
