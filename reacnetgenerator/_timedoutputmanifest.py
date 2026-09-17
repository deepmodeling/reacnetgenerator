# SPDX-License-Identifier: LGPL-3.0-or-later
"""Canonical semantic manifests for validated timeline artifacts."""

import hashlib
import hmac
import json
import struct
from collections.abc import Mapping
from dataclasses import asdict

import h5py
import numpy as np

from ._timedoutputvalidate import _root_metadata, validate_timed_output

_MANIFEST_FORMAT = "reacnetgenerator-timeline-semantic-manifest"
_MANIFEST_VERSION = "1.0"
_PATH_PARAMETERS = {
    "atomroutefilename",
    "imagefilename",
    "inputfilename",
    "jsonfilename",
    "moleculefilename",
    "moleculetimelinefilename",
    "output_dir",
    "reactionabcdfilename",
    "reactioneventfilename",
    "reactionfilename",
    "resultfilename",
    "speciesfilename",
    "tablefilename",
    "timed_output",
}
_DATASETS = (
    "sources/path",
    "sources/size_bytes",
    "sources/mtime_ns",
    "frames/source_id",
    "frames/source_frame",
    "frames/timestep",
    "atoms/type",
    "atoms/type_name",
    "species/name",
    "molecules/species_id",
    "molecules/atom_offsets",
    "molecules/atom_index",
    "molecules/bond_offsets",
    "molecules/bond_atom_index_1",
    "molecules/bond_atom_index_2",
    "molecules/bond_order",
    "molecule_ranges/molecule_id",
    "molecule_ranges/start_frame",
    "molecule_ranges/end_frame",
    "reaction_types/reactant",
    "reaction_types/product",
    "reaction_types/total_count",
    "reaction_events/transition",
    "reaction_events/reaction_type_id",
    "reaction_events/count",
)
_PROVENANCE_DATASETS = {
    "sources/path",
    "sources/size_bytes",
    "sources/mtime_ns",
}


def _normalized_configuration(configuration):
    """Remove location-only values while retaining interpretation settings."""
    parameters = {
        key: value
        for key, value in configuration["parameters"].items()
        if key not in _PATH_PARAMETERS
    }
    explicit = sorted(
        key
        for key in configuration["explicit_parameters"]
        if key not in _PATH_PARAMETERS
    )
    return {"parameters": parameters, "explicit_parameters": explicit}


def _dataset_digest(dataset, block_rows):
    digest = hashlib.sha256()
    string = h5py.check_string_dtype(dataset.dtype)
    for start in range(0, len(dataset), block_rows):
        if string is None:
            values = np.asarray(dataset[start : start + block_rows], dtype="<i8")
            digest.update(values.tobytes(order="C"))
        else:
            values = dataset.asstr()[start : start + block_rows]
            for value in values:
                encoded = value.encode("utf-8")
                digest.update(struct.pack("<Q", len(encoded)))
                digest.update(encoded)
    return digest.hexdigest()


def _canonical_hash(value):
    payload = json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def semantic_manifest(filename, *, include_provenance=False, block_rows=8192):
    """Build a storage-independent manifest after full contract validation."""
    summary = validate_timed_output(filename, block_rows=block_rows)
    with h5py.File(filename, "r") as file:
        configuration = _root_metadata(file)
        metadata = {
            "atom_index_convention": file.attrs["atom_index_convention"],
            "capabilities": json.loads(file.attrs["capabilities"]),
            "configuration": configuration
            if include_provenance
            else _normalized_configuration(configuration),
            "molecule_range_basis": file.attrs["molecule_range_basis"],
        }
        if include_provenance:
            metadata.update(
                created_utc=file.attrs["created_utc"],
                rng_version=file.attrs["rng_version"],
            )
        datasets = {}
        for path in _DATASETS:
            if not include_provenance and path in _PROVENANCE_DATASETS:
                continue
            dataset = file[path]
            datasets[path] = {
                "rows": len(dataset),
                "sha256": _dataset_digest(dataset, block_rows),
            }

    manifest = {
        "manifest_format": _MANIFEST_FORMAT,
        "manifest_version": _MANIFEST_VERSION,
        "timeline_schema_version": summary.schema_version,
        "provenance_included": bool(include_provenance),
        "metadata": metadata,
        "summary": asdict(summary),
        "datasets": datasets,
    }
    manifest["semantic_sha256"] = _canonical_hash(manifest)
    return manifest


def _pointer_token(value):
    return str(value).replace("~", "~0").replace("/", "~1")


def _differences(left, right, path):
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        for key in sorted(set(left) | set(right)):
            child = f"{path}/{_pointer_token(key)}"
            if key not in left or key not in right:
                yield child
            else:
                yield from _differences(left[key], right[key], child)
        return
    if isinstance(left, list) and isinstance(right, list):
        for position in range(max(len(left), len(right))):
            child = f"{path}/{position}"
            if position >= len(left) or position >= len(right):
                yield child
            else:
                yield from _differences(left[position], right[position], child)
        return
    if type(left) is not type(right) or left != right:
        yield path or "/"


def compare_semantic_manifests(left, right):
    """Validate manifest envelopes and locate all value differences."""
    for name, manifest in (("left", left), ("right", right)):
        if not isinstance(manifest, Mapping):
            raise ValueError(f"The {name} manifest must be a mapping")
        if manifest.get("manifest_format") != _MANIFEST_FORMAT:
            raise ValueError(f"The {name} manifest has an unsupported format")
        if manifest.get("manifest_version") != _MANIFEST_VERSION:
            raise ValueError(f"The {name} manifest has an unsupported version")
        checksum = manifest.get("semantic_sha256")
        content = dict(manifest)
        content.pop("semantic_sha256", None)
        expected = _canonical_hash(content)
        if not isinstance(checksum, str) or not hmac.compare_digest(checksum, expected):
            raise ValueError(f"The {name} manifest checksum is invalid")
    return tuple(_differences(left, right, ""))
