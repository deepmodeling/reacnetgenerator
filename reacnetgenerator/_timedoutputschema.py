# SPDX-License-Identifier: LGPL-3.0-or-later
"""Installed machine-readable contract for the timeline format."""

import json
from importlib.resources import files


def read_schema_descriptor():
    """Return a fresh mapping so callers cannot mutate shared schema state."""
    resource = files("reacnetgenerator.schemas").joinpath("timed-output-schema.json")
    with resource.open(encoding="utf-8") as file:
        return json.load(file)
