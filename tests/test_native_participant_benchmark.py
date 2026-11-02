# SPDX-License-Identifier: LGPL-3.0-or-later
"""Check benchmark windows, platform defaults and validated result publication."""

import hashlib
import importlib.util
import json
import multiprocessing
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def benchmark():
    """Load the standalone benchmark without executing its CLI."""
    spec = importlib.util.spec_from_file_location(
        "native_participant_benchmark",
        ROOT / "benchmarks/native_participant_metadata.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("start,frames", [(0, 2), (1, 2), (2, 1)])
def test_window_offset_preserves_bytes_and_hash(tmp_path, benchmark, start, frames):
    """Select frame ordinals rather than timestep values, including the last frame."""
    blocks = [f"ITEM: TIMESTEP\n{100 * i}\npayload {i}\n".encode() for i in range(3)]
    source = tmp_path / "input.dump"
    source.write_bytes(b"".join(blocks))
    target = tmp_path / "window.dump"
    digest = benchmark.extract_window(source, target, frames, start)
    expected = b"".join(blocks[start : start + frames])
    assert target.read_bytes() == expected
    assert digest == hashlib.sha256(expected).hexdigest()


@pytest.mark.parametrize("start,available", [(1, 1), (3, 0)])
def test_short_window_is_rejected(tmp_path, benchmark, start, available):
    """Do not publish a partial requested window as a complete input."""
    source = tmp_path / "input.dump"
    source.write_bytes(b"ITEM: TIMESTEP\n0\nITEM: TIMESTEP\n100\n")
    with pytest.raises(ValueError, match=f"found only {available}"):
        benchmark.extract_window(source, tmp_path / "window.dump", 2, start)


@pytest.mark.parametrize(
    "methods,expected", [(["spawn"], "spawn"), (["spawn", "fork"], "fork")]
)
def test_supported_platform_default(benchmark, monkeypatch, methods, expected):
    """Prefer fork only when the platform supports it."""
    monkeypatch.setattr(multiprocessing, "get_all_start_methods", lambda: methods)
    args = benchmark.parse_args(
        ["--trajectory", "input.dump", "--output-dir", "out", "--atom-names", "H"]
    )
    assert args.start_method == expected


def test_negative_offset_is_rejected(benchmark):
    """Reject an invalid offset before creating an output directory."""
    with pytest.raises(SystemExit):
        benchmark.parse_args(
            [
                "--trajectory",
                "input.dump",
                "--output-dir",
                "out",
                "--atom-names",
                "H",
                "--start-frame",
                "-1",
            ]
        )


@pytest.mark.parametrize("valid", [True, False])
def test_publish_only_validated_effective_parameters(
    tmp_path, benchmark, monkeypatch, capsys, valid
):
    """Run a real tiny dump; rejected native output cannot leave a measurement."""
    folder = tmp_path / "result"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark",
            "--trajectory",
            str(ROOT / "tests/inputs/reaction.dump"),
            "--output-dir",
            str(folder),
            "--atom-names",
            "H",
            "O",
            "--frames",
            "2",
            "--nproc",
            "1",
        ],
    )
    if not valid:

        def reject(path):
            raise ValueError("native validation rejected output")

        monkeypatch.setattr("reacnetgenerator.timedoutput.semantic_manifest", reject)
    previous_method = multiprocessing.get_start_method(allow_none=True)
    try:
        if valid:
            benchmark.main()
        else:
            with pytest.raises(ValueError, match="native validation rejected"):
                benchmark.main()
    finally:
        multiprocessing.set_start_method(previous_method, force=True)
    output = capsys.readouterr().out
    if valid:
        measurement = json.loads((folder / "measurement.json").read_text())
        assert measurement["requested_parameters"]["getoriginfile"] is False
        assert measurement["parameters"]["getoriginfile"] is True
        assert measurement["parameters"]["runHMM"] is False
        assert measurement["frames"] == 2
        assert (folder / "semantic-manifest.json").is_file()
        assert '"pipeline_seconds"' in output
    else:
        assert not (folder / "measurement.json").exists()
        assert not (folder / "semantic-manifest.json").exists()
        assert '"pipeline_seconds"' not in output
