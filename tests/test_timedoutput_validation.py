# SPDX-License-Identifier: LGPL-3.0-or-later
"""Public validation contracts for persisted timeline artifacts."""

import json
import shutil
from pathlib import Path

import h5py
import pytest

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._detect import _Detect
from reacnetgenerator._hmmfilter import _HMMFilter
from reacnetgenerator._path import _CollectPaths
from reacnetgenerator._timedoutput import _TimedOutputWriter
from reacnetgenerator.timedoutput import (
    TimedOutputValidationError,
    ValidationSummary,
    compare_semantic_manifests,
    iter_transition_evidence,
    read_schema_descriptor,
    semantic_manifest,
    validate_timed_output,
)
from reacnetgenerator.timedoutputcheck import main as check_timed_output

INPUT = Path(__file__).parent / "inputs" / "reaction.bond"


@pytest.fixture(scope="module")
def valid_timeline(tmp_path_factory):
    """Build one small artifact through the real analysis and writer pipeline."""
    directory = tmp_path_factory.mktemp("timed-output-validation")
    path = directory / "timeline.h5"
    rng = ReacNetGenerator(
        inputfilename=str(INPUT),
        inputfiletype="bond",
        atomname=["H", "He"],
        nproc=1,
        runHMM=False,
        output_dir=directory,
        timed_output=path,
    )
    _Detect.gettype(rng).detect()
    _HMMFilter(rng).filter()
    rng.timed_output = directory / "empty-tables.h5"
    with _TimedOutputWriter(rng):
        pass
    rng.timed_output = path
    _CollectPaths.getstype(rng).collect()
    return path


def test_validate_real_artifact_and_reject_broken_offsets(valid_timeline, tmp_path):
    """Accept writer output and identify the malformed public column by name."""
    summary = validate_timed_output(valid_timeline, block_rows=2)

    assert isinstance(summary, ValidationSummary)
    assert summary.schema_version == "1.1"
    assert summary.sources == 1
    assert summary.frames == 8
    assert summary.atoms == 3
    assert summary.atom_types == 2
    assert summary.molecules > 0

    malformed = tmp_path / "malformed.h5"
    shutil.copyfile(valid_timeline, malformed)
    with h5py.File(malformed, "r+") as file:
        file["molecules/atom_offsets"][-1] += 1

    with pytest.raises(
        TimedOutputValidationError, match=r"molecules/atom_offsets.*atom_index"
    ):
        validate_timed_output(malformed, block_rows=2)


def test_validate_accepts_declared_empty_tables(valid_timeline):
    """Keep zero-row dictionaries and payloads valid without sentinel records."""
    summary = validate_timed_output(
        valid_timeline.with_name("empty-tables.h5"), block_rows=2
    )

    assert summary.molecules == 0
    assert summary.molecule_ranges == 0
    assert summary.reaction_types == 0
    assert summary.reaction_events == 0


def test_schema_1_0_remains_readable_without_transition_evidence(
    valid_timeline, tmp_path
):
    """Keep accepted aggregate-only artifacts readable after the 1.1 extension."""
    legacy = _copy_timeline(valid_timeline, tmp_path)
    with h5py.File(legacy, "r+") as file:
        file.attrs["schema_version"] = "1.0"
        file.attrs["capabilities"] = json.dumps(["molecule_ranges", "reaction_events"])
        del file["transition_evidence"]

    assert validate_timed_output(legacy, block_rows=2).schema_version == "1.0"
    with pytest.raises(ValueError, match="does not contain transition evidence"):
        list(iter_transition_evidence(legacy, block_rows=2))


@pytest.mark.parametrize(
    "run_hmm,basis",
    [
        (True, "observed signal"),
        (False, "HMM signal"),
        ("false", "observed signal"),
    ],
)
def test_validate_rejects_conflicting_signal_basis(
    valid_timeline, tmp_path, run_hmm, basis
):
    """Reject metadata that gives one range two incompatible meanings."""
    malformed = _copy_timeline(valid_timeline, tmp_path)
    with h5py.File(malformed, "r+") as file:
        configuration = json.loads(file.attrs["configuration"])
        configuration["parameters"]["runHMM"] = run_hmm
        file.attrs["configuration"] = json.dumps(configuration)
        file.attrs["molecule_range_basis"] = basis

    with pytest.raises(
        TimedOutputValidationError,
        match=r"configuration/parameters/runHMM.*molecule_range_basis",
    ):
        validate_timed_output(malformed, block_rows=2)


def test_validate_rejects_nonstandard_json_constants(valid_timeline, tmp_path):
    """Fail before consumers disagree about non-finite configuration values."""
    malformed = _copy_timeline(valid_timeline, tmp_path)
    with h5py.File(malformed, "r+") as file:
        file.attrs["configuration"] = (
            '{"parameters":{"runHMM":false,"cutoff":NaN},"explicit_parameters":[]}'
        )

    with pytest.raises(
        TimedOutputValidationError,
        match=r"configuration must be valid JSON",
    ):
        validate_timed_output(malformed, block_rows=2)


def test_validate_rejects_a_frame_stride_that_conflicts_with_configuration(
    valid_timeline, tmp_path
):
    """Reject frame mappings that cannot result from the declared sampling stride."""
    malformed = _copy_timeline(valid_timeline, tmp_path)
    with h5py.File(malformed, "r+") as file:
        configuration = json.loads(file.attrs["configuration"])
        configuration["parameters"]["stepinterval"] = 2
        file.attrs["configuration"] = json.dumps(configuration)

    with pytest.raises(
        TimedOutputValidationError,
        match=r"frames/source_frame.*stepinterval",
    ):
        validate_timed_output(malformed, block_rows=2)


def _copy_timeline(valid_timeline, tmp_path):
    path = tmp_path / "mutated.h5"
    shutil.copyfile(valid_timeline, path)
    return path


@pytest.mark.parametrize(
    "dataset",
    [
        "sources/path",
        "frames/source_frame",
        "molecules/bond_order",
        "reaction_types/product",
    ],
)
def test_validate_requires_declared_column_types_and_alignment(
    valid_timeline, tmp_path, dataset
):
    """Reject wrong public types and row alignment before interpreting values."""
    malformed = _copy_timeline(valid_timeline, tmp_path)
    with h5py.File(malformed, "r+") as file:
        values = file[dataset][:]
        del file[dataset]
        if dataset in {"sources/path", "reaction_types/product"}:
            file.create_dataset(dataset, data=list(range(len(values))), dtype="<i8")
        else:
            file.create_dataset(dataset, data=values[:-1], dtype="<i8")

    with pytest.raises(TimedOutputValidationError, match=dataset):
        validate_timed_output(malformed, block_rows=2)


@pytest.mark.parametrize(
    "dataset,value",
    [
        ("frames/source_id", 1),
        ("atoms/type", 2),
        ("molecules/species_id", 3),
        ("molecule_ranges/molecule_id", 0),
        ("reaction_events/reaction_type_id", 1),
        ("transition_evidence/reaction_type_id", 1),
        ("transition_evidence/participant_molecule_id", 0),
    ],
)
def test_validate_rejects_dangling_references(valid_timeline, tmp_path, dataset, value):
    """Reject IDs that cannot be joined to their declared dictionary table."""
    malformed = _copy_timeline(valid_timeline, tmp_path)
    with h5py.File(malformed, "r+") as file:
        file[dataset][0] = value

    with pytest.raises(TimedOutputValidationError, match=dataset):
        validate_timed_output(malformed, block_rows=2)


def test_validate_rejects_nonmaximal_ranges_and_wrong_event_totals(
    valid_timeline, tmp_path
):
    """Enforce compact-range and aggregate-count meaning across public tables."""
    malformed = _copy_timeline(valid_timeline, tmp_path)
    with h5py.File(malformed, "r+") as file:
        ranges = file["molecule_ranges"]
        for name in ("molecule_id", "start_frame", "end_frame"):
            column = ranges[name]
            column.resize((len(column) + 1,))
        ranges["molecule_id"][-1] = ranges["molecule_id"][-2]
        ranges["start_frame"][-1] = ranges["end_frame"][-2]
        ranges["end_frame"][-1] = ranges["start_frame"][-1]

    with pytest.raises(TimedOutputValidationError, match=r"molecule_ranges.*maximal"):
        validate_timed_output(malformed, block_rows=2)

    malformed = _copy_timeline(valid_timeline, tmp_path)
    with h5py.File(malformed, "r+") as file:
        file["reaction_types/total_count"][0] += 1

    with pytest.raises(
        TimedOutputValidationError,
        match=r"reaction_types/total_count.*reaction_events",
    ):
        validate_timed_output(malformed, block_rows=2)


def test_validate_requires_one_evidence_row_per_reaction_instance(
    valid_timeline, tmp_path
):
    """Reject a structurally valid evidence table that omits its final instance."""
    malformed = _copy_timeline(valid_timeline, tmp_path)
    with h5py.File(malformed, "r+") as file:
        evidence = file["transition_evidence"]
        rows = len(evidence["transition"])
        participant_stop = int(evidence["participant_offsets"][-2])
        change_stop = int(evidence["bond_change_offsets"][-2])
        for name in ("transition", "reaction_type_id"):
            evidence[name].resize((rows - 1,))
        for name in ("participant_offsets", "bond_change_offsets"):
            evidence[name].resize((rows,))
        for name in ("participant_molecule_id", "participant_side"):
            evidence[name].resize((participant_stop,))
        for name in (
            "bond_atom_index_1",
            "bond_atom_index_2",
            "before_order",
            "after_order",
        ):
            evidence[name].resize((change_stop,))

    with pytest.raises(
        TimedOutputValidationError,
        match=r"transition evidence.*(missing|counts disagree)",
    ):
        validate_timed_output(malformed, block_rows=2)


def test_validate_recomputes_evidence_bond_changes(valid_timeline, tmp_path):
    """Reject bond differences that disagree with referenced molecule graphs."""
    malformed = _copy_timeline(valid_timeline, tmp_path)
    with h5py.File(malformed, "r+") as file:
        before = file["transition_evidence/before_order"]
        assert len(before) > 0
        before[0] = int(before[0]) + 1

    with pytest.raises(
        TimedOutputValidationError,
        match=r"transition evidence bond changes disagree",
    ):
        validate_timed_output(malformed, block_rows=2)


def test_validate_requires_participants_at_their_reaction_frames(
    valid_timeline, tmp_path
):
    """Reject a participant whose definition exists outside its claimed frame."""
    malformed = _copy_timeline(valid_timeline, tmp_path)
    with h5py.File(malformed, "r+") as file:
        evidence = file["transition_evidence"]
        molecule_id = int(
            evidence["participant_molecule_id"][int(evidence["participant_offsets"][0])]
        )
        ranges = file["molecule_ranges"]
        keep = ranges["molecule_id"][:] != molecule_id
        for name in ("molecule_id", "start_frame", "end_frame"):
            column = ranges[name]
            values = column[:][keep]
            column.resize((len(values),))
            column[:] = values

    with pytest.raises(
        TimedOutputValidationError,
        match=r"transition evidence participant is absent",
    ):
        validate_timed_output(malformed, block_rows=2)


def _repack_timeline(source, destination):
    """Copy values through h5py while deliberately changing storage layout."""
    with h5py.File(source) as old, h5py.File(destination, "w") as new:
        new.attrs.update(old.attrs)

        def copy(name, item):
            if isinstance(item, h5py.Group):
                new.require_group(name)
            else:
                parent, _, leaf = name.rpartition("/")
                group = new.require_group(parent) if parent else new
                string = h5py.check_string_dtype(item.dtype)
                if string is None:
                    group.create_dataset(
                        leaf,
                        data=item[:],
                        maxshape=(None,),
                        chunks=True,
                        compression="gzip",
                    )
                else:
                    group.create_dataset(
                        leaf,
                        data=item.asstr()[:],
                        maxshape=(None,),
                        dtype=h5py.string_dtype("utf-8"),
                        chunks=True,
                    )

        old.visititems(copy)


def test_semantic_manifest_ignores_layout_and_relocation(valid_timeline, tmp_path):
    """Identify equivalent science despite HDF5 layout and path provenance changes."""
    relocated = tmp_path / "relocated.h5"
    _repack_timeline(valid_timeline, relocated)
    with h5py.File(relocated, "r+") as file:
        file.attrs["created_utc"] = "2025-01-01T00:00:00+00:00"
        file.attrs["rng_version"] = "different-build"
        file["sources/path"][0] = "/relocated/input.bond"
        file["sources/size_bytes"][0] += 10
        file["sources/mtime_ns"][0] += 10
        configuration = json.loads(file.attrs["configuration"])
        for key in ("inputfilename", "output_dir", "timed_output"):
            configuration["parameters"][key] = f"/relocated/{key}"
        file.attrs["configuration"] = json.dumps(configuration)

    original = semantic_manifest(valid_timeline, block_rows=2)
    moved = semantic_manifest(relocated, block_rows=3)

    assert original["manifest_format"] == (
        "reacnetgenerator-timeline-semantic-manifest"
    )
    assert original["manifest_version"] == "1.0"
    assert original["semantic_sha256"] == moved["semantic_sha256"]
    assert compare_semantic_manifests(original, moved) == ()

    original_provenance = semantic_manifest(
        valid_timeline, include_provenance=True, block_rows=2
    )
    moved_provenance = semantic_manifest(
        relocated, include_provenance=True, block_rows=3
    )
    assert compare_semantic_manifests(original_provenance, moved_provenance)


def test_semantic_manifest_locates_result_changes(valid_timeline, tmp_path):
    """Report the public dataset whose values changed, independent of HDF5 bytes."""
    changed = _copy_timeline(valid_timeline, tmp_path)
    with h5py.File(changed, "r+") as file:
        file["frames/timestep"][0] += 1

    differences = compare_semantic_manifests(
        semantic_manifest(valid_timeline, block_rows=2),
        semantic_manifest(changed, block_rows=2),
    )

    assert "/datasets/frames~1timestep/sha256" in differences
    assert "/semantic_sha256" in differences


def test_semantic_manifest_preserves_analysis_settings(valid_timeline, tmp_path):
    """Report interpretation-relevant configuration changes by JSON pointer."""
    changed = _copy_timeline(valid_timeline, tmp_path)
    with h5py.File(changed, "r+") as file:
        configuration = json.loads(file.attrs["configuration"])
        configuration["parameters"]["max_component_atoms"] += 1
        file.attrs["configuration"] = json.dumps(configuration)

    differences = compare_semantic_manifests(
        semantic_manifest(valid_timeline, block_rows=2),
        semantic_manifest(changed, block_rows=2),
    )

    assert "/metadata/configuration/parameters/max_component_atoms" in differences
    assert "/semantic_sha256" in differences


def test_compare_rejects_a_manifest_with_a_stale_checksum(valid_timeline):
    """Detect manual or partial edits to a stored reference before comparison."""
    manifest = semantic_manifest(valid_timeline)
    stale = dict(manifest)
    stale["semantic_sha256"] = "0" * 64

    with pytest.raises(ValueError, match="checksum"):
        compare_semantic_manifests(manifest, stale)


def test_compare_rejects_a_non_ascii_checksum():
    """Treat malformed external checksum text as an invalid manifest."""
    malformed = {
        "manifest_format": "reacnetgenerator-timeline-semantic-manifest",
        "manifest_version": "1.0",
        "semantic_sha256": "é",
    }

    with pytest.raises(ValueError, match="The left manifest checksum is invalid"):
        compare_semantic_manifests(malformed, malformed)


def test_schema_descriptor_lists_the_public_contract(valid_timeline):
    """Expose the installed schema to consumers without requiring source docs."""
    descriptor = read_schema_descriptor()

    assert descriptor["descriptor_format"] == "reacnetgenerator-timeline-schema"
    assert descriptor["descriptor_version"] == "1.0"
    assert descriptor["timeline_format"] == "reacnetgenerator-timeline"
    assert descriptor["schema_version"] == "1.1"
    assert descriptor["datasets"]["molecules/atom_offsets"]["dtype"] == "int64-le"
    assert descriptor["datasets"]["species/name"]["dtype"] == "utf8-vlen"
    assert {item["id"] for item in descriptor["constraints"]} >= {
        "offset-endpoints",
        "frame-stride-consistency",
        "reference-integrity",
        "reaction-totals",
        "transition-evidence-totals",
        "signal-basis-consistency",
    }
    manifest = semantic_manifest(valid_timeline, include_provenance=True)
    assert set(descriptor["datasets"]) == set(manifest["datasets"])


def test_validation_and_manifest_scan_numeric_columns_in_blocks(
    valid_timeline, monkeypatch
):
    """Keep ordinary numeric scans within the caller's row budget."""
    unbounded_per_record = {
        "/molecules/atom_index",
        "/molecules/bond_atom_index_1",
        "/molecules/bond_atom_index_2",
        "/transition_evidence/participant_molecule_id",
        "/transition_evidence/participant_side",
        "/transition_evidence/bond_atom_index_1",
        "/transition_evidence/bond_atom_index_2",
        "/transition_evidence/before_order",
        "/transition_evidence/after_order",
    }
    checked = set()
    getitem = h5py.Dataset.__getitem__

    def bounded(dataset, key, *args, **kwargs):
        if (
            dataset.dtype.kind == "i"
            and dataset.name not in unbounded_per_record
            and isinstance(key, slice)
        ):
            assert key.start is not None and key.stop is not None
            assert key.stop - key.start <= 2
            checked.add(dataset.name)
        return getitem(dataset, key, *args, **kwargs)

    monkeypatch.setattr(h5py.Dataset, "__getitem__", bounded)

    validate_timed_output(valid_timeline, block_rows=2)
    semantic_manifest(valid_timeline, block_rows=2)

    assert {"/frames/timestep", "/reaction_events/count"} <= checked


def test_cli_writes_and_compares_manifests(valid_timeline, tmp_path, capsys):
    """Provide CI-friendly exit codes and atomic reference-manifest output."""
    reference = tmp_path / "reference.json"
    assert (
        check_timed_output(
            [
                str(valid_timeline),
                "--write-manifest",
                str(reference),
                "--block-rows",
                "2",
            ]
        )
        == 0
    )
    result = json.loads(capsys.readouterr().out)
    assert result["status"] == "valid"
    assert reference.exists()
    assert not list(tmp_path.glob("*.incomplete"))

    changed = _copy_timeline(valid_timeline, tmp_path)
    with h5py.File(changed, "r+") as file:
        file["frames/timestep"][0] += 1

    assert (
        check_timed_output(
            [str(changed), "--compare-manifest", str(reference), "--block-rows", "2"]
        )
        == 1
    )
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "/datasets/frames~1timestep/sha256" in captured.err


def test_cli_removes_an_incomplete_manifest_after_publish_failure(
    valid_timeline, tmp_path, monkeypatch, capsys
):
    """Clean the atomic sibling even after fdopen has closed its descriptor."""
    reference = tmp_path / "reference.json"

    def fail_replace(*args):
        raise OSError("publish failed")

    monkeypatch.setattr("reacnetgenerator.timedoutputcheck.os.replace", fail_replace)

    assert (
        check_timed_output(
            [
                str(valid_timeline),
                "--write-manifest",
                str(reference),
                "--block-rows",
                "2",
            ]
        )
        == 1
    )
    assert "publish failed" in capsys.readouterr().err
    assert not reference.exists()
    assert not list(tmp_path.glob("*.incomplete"))


def test_cli_never_replaces_an_input_artifact(valid_timeline, tmp_path, capsys):
    """Reject an output alias before an atomic JSON publish can destroy HDF5."""
    artifact = _copy_timeline(valid_timeline, tmp_path)

    assert check_timed_output([str(artifact), "--write-manifest", str(artifact)]) == 1
    assert "must not overwrite" in capsys.readouterr().err
    assert validate_timed_output(artifact).frames == 8
