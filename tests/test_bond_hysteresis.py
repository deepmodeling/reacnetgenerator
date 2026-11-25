# SPDX-License-Identifier: LGPL-3.0-or-later
"""Check stateful edge detection through the public analysis and timeline APIs."""

import json

import h5py
import pytest

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator.commandline import main_parser, parm2cmd
from reacnetgenerator.timedoutput import (
    TimedOutputValidationError,
    iter_molecule_ranges,
    iter_molecules,
    validate_timed_output,
)


def _analyze(tmp_path, sources, kind, thresholds, **options):
    target = tmp_path / "run.h5"
    rng = ReacNetGenerator(
        inputfilename=[str(p) for p in sources],
        inputfiletype=kind,
        atomname=["C"],
        bond_hysteresis=thresholds,
        output_dir=tmp_path,
        timed_output=target,
        nproc=2,
        **options,
    )
    rng.run_items(["species", "reactions"])
    validate_timed_output(target)
    with h5py.File(target) as file:
        params = json.loads(file.attrs["configuration"])["parameters"]
        assert params["bond_hysteresis"] == thresholds
        assert file.attrs["molecule_range_basis"] == (
            "HMM signal" if params["runHMM"] else "hysteresis signal"
        )
    return target


def _pair_orders(target):
    molecules = {m.molecule_id: m.bonds for m in iter_molecules(target)}
    with h5py.File(target) as file:
        orders = [0] * len(file["frames/timestep"])
    for span in iter_molecule_ranges(target):
        bonds = molecules[span.molecule_id]
        if bonds:
            ((i, j, order),) = bonds
            assert (i, j) == (0, 1)
            for frame in range(span.start_frame, span.end_frame + 1):
                orders[frame] = order
    return orders


def _reaxff_frame(step, order):
    # Some writers report slightly different orders in each direction.
    rows = (
        ["1 1 0 0 0 0 0", "2 1 0 0 0 0 0"]
        if order is None
        else [f"1 1 1 2 0 {order} 0 0 0", f"2 1 1 1 0 {order - 0.01} 0 0 0"]
    )
    return "\n".join(
        [f"# Timestep {step}", "#", "# Number of particles 2", "#", *rows, "#", ""]
    )


def test_reaxff_raw_orders_and_split(tmp_path):
    """Exercise equality, both deadband states, absent reports and re-entry."""
    sequence = [0.45, 0.6, 0.45, 0.30, 0.45, 1.8, None, 0.45]
    parts = [tmp_path / "first.bond", tmp_path / "second.bond"]
    frames = [_reaxff_frame(i, value) for i, value in enumerate(sequence)]
    parts[0].write_text("".join(frames[:2]))
    parts[1].write_text("".join(frames[2:]))
    target = _analyze(
        tmp_path, parts, "bond", {"form": 0.6, "break": 0.3}, runHMM=False, split=2
    )
    assert _pair_orders(target) == [0, 1, 1, 0, 0, 2, 0, 0]
    with h5py.File(target, "r+") as file:
        file.attrs["schema_version"] = "1.1"
    with pytest.raises(TimedOutputValidationError, match=r"1\.2"):
        validate_timed_output(target)


@pytest.mark.parametrize("use_ase", [False, True])
def test_coordinate_orders_and_periodic_images(tmp_path, use_ase):
    """Keep backend order assignment after periodic edge selection."""
    # C radii sum to 1.52 A. Only the periodic image qualifies at formation.
    distances = [1.7, 1.2, 1.7, 2.0, 1.7]
    path = tmp_path / "periodic.xyz"
    path.write_text(
        "".join(
            f'2\nLattice="8 0 0 0 8 0 0 0 8"\nC 0 0 0\nC {8 - d} 0 0\n'
            for d in distances
        )
    )
    target = _analyze(
        tmp_path,
        [path],
        "extxyz",
        {"form": 1.0, "break": 1.25},
        runHMM=False,
        use_ase=use_ase,
    )
    orders = _pair_orders(target)
    assert [bool(v) for v in orders] == [False, True, True, False, False]
    if use_ase:
        assert orders[1] == 1
    else:
        assert orders[1] > 1


def test_skipped_frames_and_hmm_provenance(tmp_path):
    """Only sampled frames feed hysteresis, independently of the downstream HMM."""
    path = tmp_path / "two.bond"
    path.write_text(
        "".join(_reaxff_frame(i, v) for i, v in enumerate([0.4, 0.8, 0.4, 0.8, 0.4]))
    )
    args = ([path], "bond", {"form": 0.6, "break": 0.3})
    target = _analyze(tmp_path, *args, runHMM=False, stepinterval=2)
    assert _pair_orders(target) == [0, 0, 0]
    _analyze(tmp_path, *args)  # Default HMM runs after per-bond detection.


@pytest.mark.parametrize(
    ("kind", "config"),
    [
        ("bond", {"form": 0.3, "break": 0.6}),
        ("bond", {"form": True, "break": 0.3}),
        ("xyz", {"form": 1, "break": float("inf")}),
        ("dump", {"form": 1, "break": 1}),
        ("extxyz", {"form": 0, "break": 2}),
        ("bond", {"form": 0.6}),
    ],
)
def test_invalid_hysteresis(kind, config):
    """Reject conflicting thresholds before starting detection."""
    with pytest.raises(ValueError, match="bond_hysteresis"):
        ReacNetGenerator(
            inputfilename="unused",
            inputfiletype=kind,
            atomname=["C"],
            bond_hysteresis=config,
        )


def test_hysteresis_cli_roundtrip():
    """Use the same pair of thresholds when reconstructing a command."""
    args = main_parser().parse_args(
        parm2cmd(
            {
                "inputfilename": "dump",
                "atomname": ["C"],
                "inputfiletype": "dump",
                "bond_hysteresis": {"form": 1.25, "break": 1.35},
            }
        )[1:]
    )
    assert args.bond_hysteresis == [1.25, 1.35]


def test_distance_equalities(tmp_path):
    """A pair forms at the inner boundary and breaks at the outer one."""
    path = tmp_path / "pair.xyz"
    path.write_text(
        "".join(f"2\nframe\nC 0 0 0\nC {d} 0 0\n" for d in [1.6, 1.52, 1.6, 1.9, 1.6])
    )
    target = _analyze(
        tmp_path,
        [path],
        "xyz",
        {"form": 1.0, "break": 1.25},
        pbc=False,
        runHMM=False,
        use_ase=True,
    )
    assert _pair_orders(target) == [0, 1, 1, 0, 0]
