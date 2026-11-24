# SPDX-License-Identifier: LGPL-3.0-or-later
"""Dependency-free public contract types for timeline validation."""

from dataclasses import dataclass

SCHEMA_VERSION = "1.1"
SUPPORTED_SCHEMA_VERSIONS = ("1.0", SCHEMA_VERSION)
TIME_UNITS = ("fs", "ps", "ns")
LENGTH_UNITS = ("angstrom", "bohr", "nm")


def _normalize_unit(value: object) -> str:
    """Use the same unit spelling in stored parameters and generated commands."""
    return str(value).strip().lower()


@dataclass(frozen=True)
class ValidationSummary:
    """Row counts established by successful timeline validation."""

    schema_version: str
    sources: int
    frames: int
    atoms: int
    atom_types: int
    species: int
    molecules: int
    molecule_ranges: int
    reaction_types: int
    reaction_events: int


class TimedOutputValidationError(ValueError):
    """A timeline violates the declared structural or semantic contract."""
