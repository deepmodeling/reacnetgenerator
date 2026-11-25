# SPDX-License-Identifier: LGPL-3.0-or-later
"""Unit tests for optional physical-unit metadata in timeline configuration."""

import json
from pathlib import Path

import h5py
import numpy as np
import pytest

from reacnetgenerator._timedoutputcontract import TimedOutputValidationError
from reacnetgenerator._timedoutputvalidate import _root_metadata
from reacnetgenerator.commandline import main_parser, parm2cmd
from reacnetgenerator.reacnetgen import ReacNetGenerator
from reacnetgenerator.timedoutput import validate_timed_output


def test_reacnetgenerator_records_optional_physical_units(tmp_path):
    """Record units only when supplied by the caller."""
    traj = tmp_path / "run.lammpstrj"
    traj.write_text("", encoding="utf-8")
    rng = ReacNetGenerator(
        inputfilename=str(traj),
        atomname=["C", "H", "O", "Cl"],
        inputfiletype="lammpsdumpfile",
        timestep_value=0.0001,
        time_unit=" PS ",
        length_unit="ANGSTROM",
    )
    prov = rng.parameter_provenance()
    assert prov["parameters"]["timestep_value"] == pytest.approx(0.0001)
    assert prov["parameters"]["time_unit"] == "ps"
    assert prov["parameters"]["length_unit"] == "angstrom"
    assert {"timestep_value", "time_unit", "length_unit"} <= set(
        prov["explicit_parameters"]
    )
    default = ReacNetGenerator(
        inputfilename=str(traj), atomname=["C"], inputfiletype="lammpsdumpfile"
    ).parameter_provenance()
    assert (
        not {"timestep_value", "time_unit", "length_unit"}
        & default["parameters"].keys()
    )


@pytest.mark.parametrize(
    "value",
    [
        -0.1,
        0.0,
        float("inf"),
        float("nan"),
        True,
        np.bool_(True),
        pytest.param(10**309, id="integer-overflow"),
    ],
)
def test_reacnetgenerator_rejects_invalid_timestep_value(tmp_path, value):
    """Report invalid durations through the public validation error."""
    traj = tmp_path / "run.lammpstrj"
    with pytest.raises(ValueError, match="timestep_value"):
        ReacNetGenerator(
            inputfilename=str(traj),
            atomname=["C"],
            inputfiletype="lammpsdumpfile",
            timestep_value=value,
        )


def test_reacnetgenerator_rejects_timestep_value_that_underflows(tmp_path):
    """Reject durations that cannot be stored as positive Python floats."""
    value = np.finfo(np.longdouble).tiny
    if float(value) != 0.0:
        pytest.skip("longdouble has no normalized positive value below float range")
    with pytest.raises(ValueError, match="timestep_value"):
        ReacNetGenerator(
            inputfilename=str(tmp_path / "run.lammpstrj"),
            atomname=["C"],
            inputfiletype="lammpsdumpfile",
            timestep_value=value,
        )


def test_reacnetgenerator_rejects_invalid_units(tmp_path):
    """Reject unsupported time and length units."""
    traj = tmp_path / "run.lammpstrj"
    with pytest.raises(ValueError, match="time_unit"):
        ReacNetGenerator(
            inputfilename=str(traj),
            atomname=["C"],
            inputfiletype="lammpsdumpfile",
            time_unit="years",
        )
    with pytest.raises(ValueError, match="length_unit"):
        ReacNetGenerator(
            inputfilename=str(traj),
            atomname=["C"],
            inputfiletype="lammpsdumpfile",
            length_unit="miles",
        )


@pytest.mark.parametrize("timestep_value", [0.0001, np.float32(0.1)])
def test_cli_parses_and_serializes_physical_units(timestep_value):
    """Expose the optional metadata through both CLI directions."""
    parser = main_parser()
    args = parser.parse_args(
        [
            "-i",
            "run.lammpstrj",
            "-a",
            "C",
            "H",
            "--timestep-value",
            "0.0001",
            "--time-unit",
            "ps",
            "--length-unit",
            "angstrom",
        ]
    )
    assert args.timestep_value == pytest.approx(0.0001)
    assert args.time_unit == "ps"
    assert args.length_unit == "angstrom"

    cmd = parm2cmd(
        {
            "inputfilename": "run.lammpstrj",
            "atomname": ["C", "H"],
            "inputfiletype": "lammpsdumpfile",
            "timestep_value": timestep_value,
            "time_unit": " PS ",
            "length_unit": "ANGSTROM",
        }
    )
    roundtrip = parser.parse_args(cmd[1:])
    assert roundtrip.timestep_value == float(timestep_value)
    assert roundtrip.time_unit == "ps"
    assert roundtrip.length_unit == "angstrom"


@pytest.mark.parametrize("option", ["--ti", "--tim", "--time", "--timed-output"])
def test_cli_preserves_timed_output_option_prefixes(option):
    """Preserve CLI prefixes that resolved to --timed-output before unit flags."""
    args = main_parser().parse_args(
        ["-i", "run.lammpstrj", "-a", "C", option, "timeline.h5", "--time-unit", "fs"]
    )
    assert args.timed_output == "timeline.h5"
    assert args.time_unit == "fs"


def test_root_metadata_validator_accepts_valid_and_rejects_bad_units(tmp_path):
    """Maintain compatibility with timelines created without unit metadata."""

    def _make_h5(name, params):
        path = tmp_path / name
        with h5py.File(path, "w") as f:
            f.attrs.update(
                format="reacnetgenerator-timeline",
                schema_version="1.1",
                status="complete",
                capabilities=json.dumps(
                    ["molecule_ranges", "reaction_events", "transition_evidence"]
                ),
                rng_version="2.0.0",
                created_utc="2026-11-24T00:00:00+00:00",
                atom_index_convention="zero-based RNG canonical atom order",
                molecule_range_basis="observed signal",
                configuration=json.dumps(
                    {"parameters": params, "explicit_parameters": []}
                ),
            )
        return path

    good = _make_h5(
        "good.h5",
        {
            "runHMM": False,
            "stepinterval": 1,
            "timestep_value": 0.0001,
            "time_unit": "ps",
            "length_unit": "angstrom",
        },
    )
    with h5py.File(good, "r") as f:
        cfg = _root_metadata(f)
    assert cfg["parameters"]["time_unit"] == "ps"

    # Legacy file without unit keys must continue to pass.
    legacy = _make_h5("legacy.h5", {"runHMM": False, "stepinterval": 1})
    with h5py.File(legacy, "r") as f:
        _root_metadata(f)

    for key, value in (
        ("timestep_value", 0),
        ("time_unit", "invalid_unit"),
        ("length_unit", "invalid_unit"),
    ):
        bad = _make_h5(
            f"bad-{key}.h5", {"runHMM": False, "stepinterval": 1, key: value}
        )
        with h5py.File(bad, "r") as f:
            with pytest.raises(TimedOutputValidationError, match=key):
                _root_metadata(f)


def test_real_timeline_persists_physical_units(tmp_path):
    """Persist explicit units without changing the stored timestep axis."""
    trajectory = Path(__file__).parent / "inputs" / "reaction.bond"
    path = tmp_path / "timeline.h5"
    rng = ReacNetGenerator(
        inputfilename=str(trajectory),
        inputfiletype="bond",
        atomname=["H", "He"],
        nproc=1,
        runHMM=False,
        output_dir=tmp_path,
        timed_output=path,
        timestep_value=0.1,
        time_unit="fs",
        length_unit="angstrom",
    )
    rng.run_items(["reactions"])

    assert validate_timed_output(path).schema_version == "1.2"
    with h5py.File(path, "r") as file:
        parameters = json.loads(file.attrs["configuration"])["parameters"]
        assert parameters["timestep_value"] == 0.1
        assert parameters["time_unit"] == "fs"
        assert parameters["length_unit"] == "angstrom"
        np.testing.assert_array_equal(file["frames/timestep"][:], np.arange(8))
