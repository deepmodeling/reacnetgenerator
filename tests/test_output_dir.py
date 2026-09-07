# SPDX-License-Identifier: LGPL-3.0-or-later
"""Tests for explicit output-directory contracts."""

from pathlib import Path

import pytest

from reacnetgenerator import ReacNetGenerator, run
from reacnetgenerator.commandline import main_parser, parm2cmd


def test_output_dir_maps_default_artifacts(tmp_path):
    """Map default artifact names below an explicit output directory."""
    output_dir = tmp_path / "artifacts"
    generator = ReacNetGenerator(
        inputfilename=str(tmp_path / "trajectory.dump"),
        inputfiletype="dump",
        atomname=["H"],
        output_dir=output_dir,
    )

    assert output_dir.is_dir()
    assert generator.artifacts["species"] == str(output_dir / "species")
    assert generator.artifacts["reactions"] == str(output_dir / "reaction")
    assert generator.artifacts["network"] == str(output_dir / "svg")
    assert all(Path(path).parent == output_dir for path in generator.artifacts.values())


def test_output_dir_is_optional_and_preserves_legacy_names(tmp_path):
    """Keep basename-derived artifact paths when output_dir is omitted."""
    input_path = tmp_path / "trajectory.dump"
    generator = ReacNetGenerator(
        inputfilename=str(input_path),
        inputfiletype="dump",
        atomname=["H"],
    )

    assert generator.speciesfilename == f"{input_path}.species"
    assert Path(generator.artifacts["species"]).parent == tmp_path


def test_cli_and_parm2cmd_expose_output_dir(tmp_path):
    """Round-trip output_dir through CLI parsing and command rebuilding."""
    output_dir = tmp_path / "artifacts"
    parser = main_parser()
    args = parser.parse_args(
        ["-i", "trajectory.dump", "-a", "H", "--output-dir", str(output_dir)]
    )
    assert args.output_dir == str(output_dir)
    command = parm2cmd(
        {
            "inputfilename": "trajectory.dump",
            "atomname": ["H"],
            "inputfiletype": "dump",
            "output_dir": output_dir,
        }
    )
    assert parser.parse_args(command[1:]).output_dir == str(output_dir)


def test_run_rejects_unknown_items(tmp_path):
    """Reject unsupported semantic output items."""
    with pytest.raises(ValueError, match="Unsupported output items"):
        run(
            input_path="trajectory.dump",
            output_dir=tmp_path / "artifacts",
            input_type="dump",
            atomname=["H"],
            items=("species", "unknown"),
        )


def test_run_rejects_empty_items():
    """Reject an empty semantic output request."""
    with pytest.raises(ValueError, match="at least one output stage"):
        run(
            input_path="trajectory.dump",
            input_type="dump",
            atomname=["H"],
            items=(),
        )


def test_run_wrapper_returns_artifacts_without_running(tmp_path, monkeypatch):
    """Return the generator artifact mapping from the run convenience API."""

    def fake_runanddraw(self, *, run, draw, report):
        assert (run, draw, report) == (True, True, False)
        return dict(self.artifacts)

    monkeypatch.setattr(
        "reacnetgenerator.reacnetgen.ReacNetGenerator.runanddraw", fake_runanddraw
    )
    artifacts = run(
        input_path="trajectory.dump",
        output_dir=tmp_path / "artifacts",
        input_type="dump",
        atomname=["H"],
        items=("species", "network"),
    )
    assert artifacts["report"].endswith("/html") or artifacts["report"].endswith(
        "\\html"
    )


def test_run_wrapper_preserves_legacy_paths_without_output_dir(tmp_path, monkeypatch):
    """Let the convenience API omit output_dir without changing legacy paths."""

    def fake_runanddraw(self, *, run, draw, report):
        assert (run, draw, report) == (True, False, False)
        return dict(self.artifacts)

    monkeypatch.setattr(
        "reacnetgenerator.reacnetgen.ReacNetGenerator.runanddraw", fake_runanddraw
    )
    input_path = tmp_path / "trajectory.dump"
    artifacts = run(
        input_path=input_path,
        input_type="dump",
        atomname=["H"],
        items=("species",),
    )
    assert artifacts["species"] == f"{input_path}.species"


def test_run_report_includes_network_dependency(tmp_path, monkeypatch):
    """Request network generation whenever an HTML report is requested."""

    def fake_runanddraw(self, *, run, draw, report):
        assert (run, draw, report) == (True, True, True)
        return dict(self.artifacts)

    monkeypatch.setattr(
        "reacnetgenerator.reacnetgen.ReacNetGenerator.runanddraw", fake_runanddraw
    )
    run(
        input_path="trajectory.dump",
        output_dir=tmp_path / "artifacts",
        input_type="dump",
        atomname=["H"],
        items=("report",),
    )


def test_nolabel_network_stays_inside_output_dir(tmp_path, monkeypatch):
    """Prefix only the image basename for unlabeled network output."""
    from reacnetgenerator._draw import _DrawNetwork

    output_dir = tmp_path / "artifacts"
    generator = ReacNetGenerator(
        inputfilename="trajectory.dump",
        inputfiletype="dump",
        atomname=["H"],
        output_dir=output_dir,
        nolabel=True,
    )
    drawer = _DrawNetwork(generator)
    monkeypatch.setattr(drawer, "_readtable", lambda _: ([[0]], ["H"]))
    monkeypatch.setattr(drawer, "_handlespecies", lambda names: (names, {"H": "H"}))

    drawer._draw()

    assert (output_dir / "svg").is_file()
    assert (output_dir / "nolabel_svg").is_file()
    assert not (tmp_path / "nolabel_artifacts").exists()
