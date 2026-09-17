# SPDX-License-Identifier: LGPL-3.0-or-later
"""Dependency-free public contract types for timeline validation."""

from dataclasses import dataclass

SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class ValidationSummary:
    """Row counts established by a successful schema 1.0 validation."""

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
