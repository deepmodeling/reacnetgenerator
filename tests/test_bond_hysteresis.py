# SPDX-License-Identifier: LGPL-3.0-or-later
"""Check stateful edge detection through the public analysis and timeline APIs."""

import json
import math
import os
import sys
from contextlib import nullcontext
from fractions import Fraction

import h5py
import numpy as np
import pytest
from ase import Atoms
from ase.data import atomic_numbers, covalent_radii
from ase.neighborlist import neighbor_list

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator import _detect as detect_module
from reacnetgenerator._detect import _Detect, _DetectLAMMPSdump
from reacnetgenerator.commandline import _commandline, main_parser, parm2cmd
from reacnetgenerator.timedoutput import (
    TimedOutputValidationError,
    iter_molecule_ranges,
    iter_molecules,
    iter_reaction_events,
    validate_timed_output,
)
from reacnetgenerator.utils import bytestolist, read_compressed_block


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
        ("xyz", {"form": Fraction(1), "break": Fraction(1) + Fraction(1, 10**30)}),
        ("xyz", {"form": Fraction(1, 10**400), "break": 1.2}),
        ("bond", {"form": 10**400, "break": 0.3}),
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


@pytest.mark.parametrize(
    "ase_option",
    [
        {"ase_cutoff_mult": 1.3},
        {"custom_cutoffs": "C-C:1.8"},
    ],
)
def test_hysteresis_rejects_conflicting_ase_cutoff_settings(ase_option):
    """Reject combining coordinate hysteresis with single-cutoff ASE parameters."""
    with pytest.raises(ValueError, match="bond_hysteresis"):
        ReacNetGenerator(
            inputfilename="unused",
            inputfiletype="dump",
            atomname=["C"],
            bond_hysteresis={"form": 1.2, "break": 1.4},
            **ase_option,
        )


def test_dump_oscillating_bond_suppresses_spurious_reactions(tmp_path):
    """A thermally oscillating C-O bond stays bound in the deadband without spurious events."""
    # C (0.76 A) + O (0.66 A) = 1.42 A; form <= 1.704 A (1.2x), break >= 1.988 A (1.4x).
    distances = [1.45, 1.85, 1.45, 1.85, 1.45]
    frames = []
    for step, dist in enumerate(distances):
        frames.append(
            "\n".join(
                [
                    "ITEM: TIMESTEP",
                    str(step),
                    "ITEM: NUMBER OF ATOMS",
                    "2",
                    "ITEM: BOX BOUNDS pp pp pp",
                    "0.0 10.0",
                    "0.0 10.0",
                    "0.0 10.0",
                    "ITEM: ATOMS id type x y z",
                    "1 1 0.0 0.0 0.0",
                    f"2 2 {dist} 0.0 0.0",
                    "",
                ]
            )
        )
    dump = tmp_path / "oscillate.dump"
    dump.write_text("".join(frames))

    hyst_target = tmp_path / "hyst.h5"
    rng_hyst = ReacNetGenerator(
        inputfilename=str(dump),
        inputfiletype="dump",
        atomname=["C", "O"],
        bond_hysteresis={"form": 1.2, "break": 1.4},
        runHMM=False,
        use_ase=True,
        output_dir=tmp_path / "hyst",
        timed_output=hyst_target,
        nproc=2,
    )
    rng_hyst.run_items(["species", "reactions"])
    validate_timed_output(hyst_target)
    assert list(iter_reaction_events(hyst_target)) == []

    base_target = tmp_path / "base.h5"
    rng_base = ReacNetGenerator(
        inputfilename=str(dump),
        inputfiletype="dump",
        atomname=["C", "O"],
        runHMM=False,
        use_ase=True,
        output_dir=tmp_path / "base",
        timed_output=base_target,
        nproc=2,
    )
    rng_base.run_items(["species", "reactions"])
    assert len(list(iter_reaction_events(base_target))) > 0

    with h5py.File(hyst_target, "r+") as file:
        file.attrs["molecule_range_basis"] = "observed signal"
    with pytest.raises(TimedOutputValidationError, match="molecule_range_basis"):
        validate_timed_output(hyst_target)


def test_hysteresis_cli_roundtrip():
    """Use the same pair of thresholds when reconstructing a command."""
    args = main_parser().parse_args(
        parm2cmd(
            {
                "inputfilename": "dump",
                "atomname": ["C"],
                "inputfiletype": "dump",
                "bond_hysteresis": {
                    "form": 1.25,
                    "break": 1.35,
                    "pair_break": {"H-O": 1.3, "H-H": 1.75},
                },
            }
        )[1:]
    )
    assert args.bond_hysteresis == [1.25, 1.35]
    assert args.bond_break_mult == [["H-H", "1.75"], ["H-O", "1.3"]]


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


@pytest.mark.parametrize("ase_orders", [False, True])
def test_pair_break_multipliers(tmp_path, ase_orders):
    """Extend H-H, shorten H-O and keep the global O-O cutoff across frames."""
    # The first pair crosses the periodic boundary and all pairs start unbound.
    distances = [
        (0.7, 1.1, 1.5),
        (0.95, 1.30, 1.9),
        (1.10, 1.24, 2.0),
        (0.85, 1.1, 1.7),
    ]
    trajectory = tmp_path / "pairs.xyz"
    trajectory.write_text(
        "".join(
            f'6\nLattice="20 0 0 0 20 0 0 0 20"\nH 0 0 0\nH {20 - hh} 0 0\n'
            f"O 5 0 0\nH {5 + ho} 0 0\nO 10 0 0\nO {10 + oo} 0 0\n"
            for hh, ho, oo in distances
        )
    )
    target = tmp_path / "pairs.h5"
    config = {"form": 1.25, "break": 1.5, "pair_break": {"H-H": 1.75, "O-H": 1.3}}
    rng = ReacNetGenerator(
        inputfilename=str(trajectory),
        inputfiletype="extxyz",
        atomname=["H", "O"],
        output_dir=tmp_path,
        timed_output=target,
        bond_hysteresis=config,
        use_ase=ase_orders,
        runHMM=False,
        nproc=1,
    )
    rng.run_items(["species", "reactions"])
    validate_timed_output(target)
    with h5py.File(target) as file:
        assert json.loads(file.attrs["configuration"])["parameters"][
            "bond_hysteresis"
        ] == {"form": 1.25, "break": 1.5, "pair_break": {"H-H": 1.75, "H-O": 1.3}}
    edges = [set() for _ in distances]
    definitions = {m.molecule_id: m.bonds for m in iter_molecules(target)}
    for span in iter_molecule_ranges(target):
        for frame in range(span.start_frame, span.end_frame + 1):
            edges[frame].update((a, b) for a, b, _ in definitions[span.molecule_id])
    assert edges == [{(0, 1), (2, 3), (4, 5)}, {(0, 1), (4, 5)}, set(), {(2, 3)}]


@pytest.mark.parametrize(
    ("kind", "pairs"),
    [
        ("bond", {"H-H": 2.0}),
        ("dump", {}),
        ("dump", None),
        ("dump", {1: 1.6}),
        ("dump", {"H-H": True}),
        ("dump", {"H-Qq": 1.6}),
        ("dump", {"X-H": 1.75}),
        ("dump", {"H-H": 1.25}),
        ("dump", {"H-H": -1.5}),
        ("dump", {"H-H": float("inf")}),
        ("dump", {"H-H": float("nan")}),
        ("dump", {"H-H": 10**400}),
        ("dump", {"O-H": 1.6, "H-O": 1.7}),
        ("dump", {"H": 1.7}),
    ],
)
def test_invalid_pair_break_multipliers(kind, pairs):
    """Reject invalid boundary values, non-elements and ambiguous declarations."""
    with pytest.raises(ValueError, match="bond_hysteresis"):
        ReacNetGenerator(
            inputfilename="unused",
            inputfiletype=kind,
            atomname=["H", "O"],
            bond_hysteresis={
                "form": 1.25,
                "break": 0.3 if kind == "bond" else 1.5,
                "pair_break": pairs,
            },
        )


def test_pair_break_boundary_canonical_order_and_search_scope(tmp_path):
    """Break at exact pair_break equality, sort keys canonically, and scope search radii."""
    # H-H covalent radius sum is 0.62 A; form=1.25 -> 0.775 A, break=1.75 -> 1.085 A.
    path = tmp_path / "hh.xyz"
    path.write_text(
        "".join(f"2\nframe\nH 0 0 0\nH {d} 0 0\n" for d in [0.7, 1.084, 1.085, 0.9])
    )
    target = tmp_path / "hh.h5"
    rng = ReacNetGenerator(
        inputfilename=str(path),
        inputfiletype="xyz",
        atomname=["H", "O", "C"],
        pbc=False,
        runHMM=False,
        use_ase=True,
        output_dir=tmp_path,
        timed_output=target,
        nproc=1,
        bond_hysteresis={
            "form": 1.25,
            "break": 1.35,
            "pair_break": {"O-H": 1.3, "H-H": 1.75},
        },
    )
    assert rng.bond_hysteresis is not None
    assert list(rng.bond_hysteresis["pair_break"]) == ["H-H", "H-O"]
    detector = _Detect.gettype(rng)
    assert detector._atom_outer is not None
    assert detector._atom_outer[atomic_numbers["H"]] == pytest.approx(1.75)
    assert detector._atom_outer[atomic_numbers["O"]] == pytest.approx(1.35)
    assert detector._atom_outer[atomic_numbers["C"]] == pytest.approx(1.35)

    rng.run_items(["species", "reactions"])
    validate_timed_output(target)
    assert _pair_orders(target) == [1, 1, 0, 0]

    # Validator rejects tampered pair_break metadata in stored schema 1.2 timelines.
    with h5py.File(target, "r+") as file:
        config = json.loads(file.attrs["configuration"])
        config["parameters"]["bond_hysteresis"]["pair_break"] = {"H-H": 1.25}
        file.attrs["configuration"] = json.dumps(config)
    with pytest.raises(TimedOutputValidationError, match="pair_break"):
        validate_timed_output(target)


def test_bond_break_mult_cli_execution(tmp_path, monkeypatch):
    """Exercise _commandline parsing of --bond-hysteresis and --bond-break-mult."""
    path = tmp_path / "cli_pairs.xyz"
    path.write_text(
        "".join(f"2\nframe\nH 0 0 0\nH {d} 0 0\n" for d in [0.7, 0.95, 1.1])
    )
    target = tmp_path / "cli_pairs.h5"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "reacnetgenerator",
            "--type",
            "xyz",
            "-i",
            str(path),
            "-a",
            "H",
            "O",
            "--nopbc",
            "--nohmm",
            "--use-ase",
            "-n",
            "1",
            "--output-dir",
            str(tmp_path),
            "--timed-output",
            str(target),
            "--items",
            "species,reactions",
            "--bond-hysteresis",
            "1.25",
            "1.35",
            "--bond-break-mult",
            "O-H",
            "1.3",
            "--bond-break-mult",
            "H-H",
            "1.75",
        ],
    )
    _commandline()
    validate_timed_output(target)
    with h5py.File(target) as file:
        hysteresis = json.loads(file.attrs["configuration"])["parameters"][
            "bond_hysteresis"
        ]
    assert hysteresis == {
        "form": 1.25,
        "break": 1.35,
        "pair_break": {"H-H": 1.75, "H-O": 1.3},
    }
    assert list(hysteresis["pair_break"]) == ["H-H", "H-O"]
    assert _pair_orders(target) == [1, 1, 0]


@pytest.mark.parametrize(
    ("extra_args", "error_pattern"),
    [
        (["--bond-break-mult", "H-H", "1.75"], "--bond-break-mult requires"),
        (
            [
                "--bond-hysteresis",
                "1.25",
                "1.35",
                "--bond-break-mult",
                "H-H",
                "1.6",
                "--bond-break-mult",
                "H-H",
                "1.75",
            ],
            "duplicate --bond-break-mult for H-H",
        ),
        (
            [
                "--bond-hysteresis",
                "1.25",
                "1.35",
                "--bond-break-mult",
                "O-H",
                "1.6",
                "--bond-break-mult",
                "H-O",
                "1.75",
            ],
            "duplicate --bond-break-mult for H-O",
        ),
        (
            [
                "--bond-hysteresis",
                "1.25",
                "1.35",
                "--bond-break-mult",
                "H-H",
                "not-a-number",
            ],
            "invalid --bond-break-mult multiplier",
        ),
    ],
)
def test_bond_break_mult_cli_errors(monkeypatch, capsys, extra_args, error_pattern):
    """Reject missing base hysteresis, duplicates in either order, and bad numbers."""
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "reacnetgenerator",
            "--type",
            "dump",
            "-i",
            "unused.dump",
            "-a",
            "H",
            "O",
            *extra_args,
        ],
    )
    with pytest.raises(SystemExit) as exc:
        _commandline()
    assert exc.value.code == 2
    assert error_pattern in capsys.readouterr().err


def _run_hysteresis_sequence(
    atomname,
    frames,
    *,
    pbc,
    use_ase,
    thresholds,
    force_neighbor_list=False,
    monkeypatch=None,
):
    rng = ReacNetGenerator(
        inputfilename="unused",
        inputfiletype="dump",
        atomname=atomname,
        pbc=pbc,
        use_ase=use_ase,
        bond_hysteresis=thresholds,
    )
    detector = _DetectLAMMPSdump(rng)
    detector._active_bonds = set()
    if force_neighbor_list and monkeypatch is not None:
        monkeypatch.setattr(
            detector, "_get_kdtree_hysteresis_pairs", lambda *args, **kwargs: None
        )
    return [
        detector._get_hysteresis_bonds(atoms.copy(), cell) for atoms, cell in frames
    ]


@pytest.mark.parametrize("use_ase", [False, True])
def test_kdtree_hysteresis_matches_neighbor_list_periodic_and_boundary(
    monkeypatch, use_ase
):
    """Match ASE neighbor_list across periodic boundaries and negative-epsilon wraps."""
    cell = np.diag([10.0, 10.0, 10.0])
    # Includes -1e-16 wrapping to box edge and cross-boundary H-H/C-O pairs.
    frames = [
        (
            Atoms(
                "H2CO",
                positions=[
                    [-1e-16, 1.0, 1.0],
                    [9.25, 1.0, 1.0],
                    [5.0, 5.0, 5.0],
                    [6.30, 5.0, 5.0],
                ],
            ),
            cell,
        ),
        (
            Atoms(
                "H2CO",
                positions=[
                    [0.05, 1.0, 1.0],
                    [9.02, 1.0, 1.0],
                    [5.0, 5.0, 5.0],
                    [6.90, 5.0, 5.0],
                ],
            ),
            cell,
        ),
        (
            Atoms(
                "H2CO",
                positions=[
                    [0.05, 1.0, 1.0],
                    [8.80, 1.0, 1.0],
                    [5.0, 5.0, 5.0],
                    [7.20, 5.0, 5.0],
                ],
            ),
            cell,
        ),
    ]
    thresholds = {"form": 1.25, "break": 1.35, "pair_break": {"H-H": 1.75, "C-O": 1.45}}
    fast = _run_hysteresis_sequence(
        ["H", "C", "O"], frames, pbc=True, use_ase=use_ase, thresholds=thresholds
    )
    ref = _run_hysteresis_sequence(
        ["H", "C", "O"],
        frames,
        pbc=True,
        use_ase=use_ase,
        thresholds=thresholds,
        force_neighbor_list=True,
        monkeypatch=monkeypatch,
    )
    assert fast == ref


@pytest.mark.parametrize("use_ase", [False, True])
def test_kdtree_hysteresis_matches_neighbor_list_nonperiodic_pair_break(
    monkeypatch, use_ase
):
    """Match ASE neighbor_list for non-periodic trajectories with pair overrides."""
    cell = np.diag([20.0, 20.0, 20.0])
    frames = [
        (
            Atoms("H2O2", positions=[[0, 0, 0], [0.7, 0, 0], [5, 0, 0], [6.4, 0, 0]]),
            cell,
        ),
        (
            Atoms("H2O2", positions=[[0, 0, 0], [1.0, 0, 0], [5, 0, 0], [6.85, 0, 0]]),
            cell,
        ),
        (
            Atoms("H2O2", positions=[[0, 0, 0], [1.15, 0, 0], [5, 0, 0], [6.85, 0, 0]]),
            cell,
        ),
    ]
    thresholds = {"form": 1.25, "break": 1.35, "pair_break": {"H-H": 1.75}}
    fast = _run_hysteresis_sequence(
        ["H", "O"], frames, pbc=False, use_ase=use_ase, thresholds=thresholds
    )
    ref = _run_hysteresis_sequence(
        ["H", "O"],
        frames,
        pbc=False,
        use_ase=use_ase,
        thresholds=thresholds,
        force_neighbor_list=True,
        monkeypatch=monkeypatch,
    )
    assert fast == ref


def test_kdtree_hysteresis_falls_back_for_triclinic_and_narrow_cells(monkeypatch):
    """Fall back to neighbor_list on triclinic or narrow periodic cells."""
    triclinic = np.array([[10.0, 0.0, 0.0], [1.5, 10.0, 0.0], [0.0, 0.0, 10.0]])
    narrow = np.diag([2.0, 10.0, 10.0])
    atoms = Atoms("HH", positions=[[0.0, 0.0, 0.0], [0.72, 0.0, 0.0]])
    thresholds = {"form": 1.25, "break": 1.35, "pair_break": {"H-H": 1.75}}

    rng = ReacNetGenerator(
        inputfilename="unused",
        inputfiletype="dump",
        atomname=["H"],
        pbc=True,
        use_ase=True,
        bond_hysteresis=thresholds,
    )
    detector = _DetectLAMMPSdump(rng)
    detector._active_bonds = set()
    calls = []
    orig_nl = __import__(
        "reacnetgenerator._detect", fromlist=["neighbor_list"]
    ).neighbor_list

    def spy_nl(*args, **kwargs):
        calls.append(True)
        return orig_nl(*args, **kwargs)

    monkeypatch.setattr("reacnetgenerator._detect.neighbor_list", spy_nl)
    res_tri = detector._get_hysteresis_bonds(atoms.copy(), triclinic)
    res_narrow = detector._get_hysteresis_bonds(atoms.copy(), narrow)
    assert len(calls) == 2
    assert res_tri == ([[1], [0]], [[1], [1]])
    assert res_narrow == ([[1], [0]], [[1], [1]])


def test_hysteresis_frame_caches_and_65536_widening(tmp_path, monkeypatch):
    """Reuse repeated hysteresis frames safely across the uint16 step boundary."""
    detector_rng = ReacNetGenerator(
        inputfilename=str(tmp_path / "dummy.xyz"),
        inputfiletype="xyz",
        atomname=["C", "O"],
        pbc=False,
        use_ase=True,
        bond_hysteresis={"form": 1.2, "break": 1.4},
        output_dir=tmp_path,
        timed_output=tmp_path / "t.h5",
        nproc=1,
    )
    path = tmp_path / "dummy.xyz"
    # Exercise widening both on a changed graph (frame 2) and on an unchanged
    # frame with partial array rollback (frame 4).
    frame_defs = [
        (1.45, 1.45),
        (1.45, 2.20),
        (1.45, 1.45),
        (1.45, 2.20),
        (1.45, 2.20),
        (1.45, 2.20),
    ]
    path.write_text(
        "".join(
            f"4\nf{idx}\nC 0 0 0\nO {d1} 0 0\nC 8 0 0\nO {8 + d2} 0 0\n"
            for idx, (d1, d2) in enumerate(frame_defs)
        )
    )
    detector = _Detect.gettype(detector_rng)
    orig_readstep = detector._readstepfunc
    step_map = [0, 1, 65536, 2, 65537, 65538]
    monkeypatch.setattr(
        detector,
        "_readstepfunc",
        lambda item: (
            orig_readstep(item)[0],
            (step_map[item[0]], step_map[item[0]]),
        ),
    )
    detector.detect()
    try:
        with open(detector_rng.moleculetempfilename, "rb") as handle:
            blocks = list(read_compressed_block(handle))
        # 4 definitions: (0,1), (2,3), (2,), (3,)
        assert detector_rng.temp1it == 4
        observed_histories = {
            tuple(bytestolist(blocks[i])): bytestolist(blocks[i + 3]).tolist()
            for i in range(0, len(blocks), 4)
        }
        assert observed_histories == {
            (0, 1): [0, 1, 65536, 2, 65537, 65538],
            (2, 3): [0, 65536],
            (2,): [1, 2, 65537, 65538],
            (3,): [1, 2, 65537, 65538],
        }
    finally:
        os.remove(detector_rng.moleculetempfilename)


def test_fast_dump_parser_hysteresis_parity_and_fallback(tmp_path):
    """Match line-by-line parsing across sorted, shuffled, triclinic and non-numeric dumps."""
    frames_spec = [
        # (ts, box, columns_header, rows)
        (
            0,
            ("0.0 15.0", "0.0 15.0", "0.0 15.0"),
            "ITEM: ATOMS id type x y z",
            [
                "1 1 0.0 0.0 0.0",
                "2 2 1.45 0.0 0.0",
                "3 1 8.0 0.0 0.0",
                "4 2 9.45 0.0 0.0",
            ],
        ),
        (
            10,
            ("0.0 16.0 1.2", "0.0 15.0 -0.5", "0.0 15.0 0.3"),
            "ITEM: ATOMS id type x y z",
            [
                "4 2 10.20 0.0 0.0",
                "2 2 1.85 0.0 0.0",
                "1 1 0.0 0.0 0.0",
                "3 1 8.0 0.0 0.0",
            ],
        ),
        (
            20,
            ("0.0 15.0", "0.0 15.0", "0.0 15.0"),
            "ITEM: ATOMS id type x y z",
            [
                "1 1 0.0 0.0 0.0",
                "2 2 2.20 0.0 0.0",
                "3 1 8.0 0.0 0.0",
                "4 2 9.45 0.0 0.0",
            ],
        ),
    ]
    dump_path = tmp_path / "mixed.dump"
    dump_path.write_text(
        "".join(
            "\n".join(
                [
                    "ITEM: TIMESTEP",
                    str(ts),
                    "ITEM: NUMBER OF ATOMS",
                    "4",
                    "ITEM: BOX BOUNDS xy xz yz pp pp pp"
                    if len(box[0].split()) == 3
                    else "ITEM: BOX BOUNDS pp pp pp",
                    *box,
                    hdr,
                    *rows,
                    "",
                ]
            )
            for ts, box, hdr, rows in frames_spec
        )
    )
    # Non-numeric extra column triggers clean fallback to line-by-line parser.
    fallback_path = tmp_path / "fallback.dump"
    fallback_path.write_text(
        "".join(
            "\n".join(
                [
                    "ITEM: TIMESTEP",
                    str(ts),
                    "ITEM: NUMBER OF ATOMS",
                    "4",
                    "ITEM: BOX BOUNDS xy xz yz pp pp pp"
                    if len(box[0].split()) == 3
                    else "ITEM: BOX BOUNDS pp pp pp",
                    *box,
                    "ITEM: ATOMS id type x y z tag",
                    *(f"{row} ok" for row in rows),
                    "",
                ]
            )
            for ts, box, _hdr, rows in frames_spec
        )
    )

    def run_dump(src):
        target = tmp_path / f"{src.stem}.h5"
        rng = ReacNetGenerator(
            inputfilename=str(src),
            inputfiletype="dump",
            atomname=["C", "O"],
            pbc=True,
            use_ase=True,
            runHMM=False,
            bond_hysteresis={"form": 1.2, "break": 1.4},
            output_dir=tmp_path / src.stem,
            timed_output=target,
            nproc=1,
        )
        rng.run_items(["species", "reactions"])
        validate_timed_output(target)
        mols = {m.molecule_id: (m.atom_index, m.bonds) for m in iter_molecules(target)}
        spans = [
            (s.start_frame, s.end_frame, mols[s.molecule_id])
            for s in iter_molecule_ranges(target)
        ]
        return sorted(spans)

    assert run_dump(dump_path) == run_dump(fallback_path)


def test_hysteresis_openbabel_aromatic_and_kdtree_edge_fallbacks(monkeypatch):
    """Cover Open Babel aromatic order 5->12 and kdtree/radii edge fallbacks."""
    angles = np.arange(6) * np.pi / 3.0
    ring = np.column_stack((np.cos(angles), np.sin(angles), np.zeros(6)))
    benzene = Atoms("C6H6", positions=np.vstack((1.39 * ring, 2.48 * ring)) + 10.0)
    cell = np.diag([20.0, 20.0, 20.0])

    rng_ob = ReacNetGenerator(
        inputfilename="unused",
        inputfiletype="dump",
        atomname=["C", "H"],
        pbc=True,
        use_ase=False,
        bond_hysteresis={"form": 1.25, "break": 1.35},
    )
    det_ob = _DetectLAMMPSdump(rng_ob)
    det_ob._active_bonds = {(0, 1)}
    bonds, orders = det_ob._get_hysteresis_bonds(benzene, cell)
    assert any(2 in atom_orders for atom_orders in orders)
    assert len(bonds) == 12

    class _OrderFiveOBMol(detect_module.openbabel.OBMol):
        def GetBond(self, idx):
            bond = super().GetBond(idx)
            bond.SetBondOrder(5)
            return bond

    monkeypatch.setattr(detect_module.openbabel, "OBMol", _OrderFiveOBMol)
    det_ob_five = _DetectLAMMPSdump(
        ReacNetGenerator(
            inputfilename="unused",
            inputfiletype="dump",
            atomname=["C"],
            pbc=False,
            use_ase=False,
            bond_hysteresis={"form": 1.25, "break": 1.35},
        )
    )
    det_ob_five._active_bonds = set()
    bonds5, orders5 = det_ob_five._get_hysteresis_bonds(
        Atoms("CC", positions=[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]), cell
    )
    assert bonds5 == [[1], [0]]
    assert orders5 == [[12], [12]]
    monkeypatch.undo()

    # Single atom, unknown atomic number 0, and non-finite coordinates.
    det_ase = _DetectLAMMPSdump(
        ReacNetGenerator(
            inputfilename="unused",
            inputfiletype="dump",
            atomname=["H"],
            pbc=False,
            use_ase=True,
            bond_hysteresis={"form": 1.25, "break": 1.35},
        )
    )
    det_ase._active_bonds = set()
    assert det_ase._get_hysteresis_bonds(Atoms("H", positions=[[0, 0, 0]]), cell) == (
        [[]],
        [[]],
    )
    assert det_ase._hysteresis_radii_and_cutoff(np.array([0, 1])) is None
    assert (
        det_ase._get_kdtree_hysteresis_pairs(
            np.array([1, 1]), np.array([[0.0, 0.0, 0.0], [np.nan, 0.0, 0.0]]), cell, 1.5
        )
        is None
    )


def test_hysteresis_cache_eviction_source_errors_and_dump_fallbacks(
    tmp_path, monkeypatch
):
    """Cover bounded record cache eviction, source-frame errors, and dump parser fallbacks."""
    rng = ReacNetGenerator(
        inputfilename="unused",
        inputfiletype="dump",
        atomname=["C", "O"],
        pbc=False,
        use_ase=True,
        bond_hysteresis={"form": 1.2, "break": 1.4},
        output_dir=tmp_path,
        timed_output=tmp_path / "t.h5",
        nproc=1,
    )
    det = _DetectLAMMPSdump(rng)
    det.N = 3
    det.atomtype = np.array([0, 0, 1])
    det._active_bonds = set()
    monkeypatch.setattr(detect_module, "_DETECT_HYSTERESIS_RECORD_CACHE_MAX", 2)
    recs = det._connectmolecule([[], [], []], [[], [], []], frame=0, timestep=0)
    assert len(recs) == 3 and len(det._record_bytes_cache) == 2

    # Byte-budget eviction and oversized component bypass.
    det._clear_hysteresis_state()
    det._active_bonds = set()
    det._record_bytes_cache = {}
    monkeypatch.undo()
    # Two small single-atom entries fit, third evicts oldest when byte budget is tight;
    # a 3-atom ringed component whose key_cost exceeds max_entry bypasses without evicting.
    monkeypatch.setattr(
        detect_module, "_DETECT_HYSTERESIS_RECORD_CACHE_MAX_ENTRY_BYTES", 700
    )
    monkeypatch.setattr(detect_module, "_DETECT_HYSTERESIS_RECORD_CACHE_MAX_BYTES", 900)
    small_recs = det._connectmolecule([[], [], []], [[], [], []], frame=0, timestep=0)
    assert len(small_recs) == 3 and len(det._record_bytes_cache) == 1
    assert det._record_bytes_cache_bytes <= 900
    prev_keys = list(det._record_bytes_cache.keys())
    big_recs = det._connectmolecule(
        [[1, 2], [0, 2], [0, 1]],
        [[1, 1], [1, 1], [1, 1]],
        frame=1,
        timestep=1,
    )
    assert len(big_recs) == 1
    assert list(det._record_bytes_cache.keys()) == prev_keys
    monkeypatch.setattr(
        detect_module, "_DETECT_HYSTERESIS_RECORD_CACHE_MAX_ENTRY_BYTES", 400
    )
    # 0-bond single atom has key_cost=320 (<400) but key_cost + len(encoded) > 400.
    det._connectmolecule([[]], [[]], frame=2, timestep=2)
    assert list(det._record_bytes_cache.keys()) == prev_keys
    monkeypatch.undo()

    # Stdin helper delegation and _source_frames incomplete-frame errors.
    stdin_sentinel = object()
    monkeypatch.setattr(
        detect_module.fileinput, "FileInput", lambda files: nullcontext(stdin_sentinel)
    )
    with det._open_source_file("-") as opened:
        assert opened is stdin_sentinel
    monkeypatch.undo()

    p1 = tmp_path / "a.txt"
    p2 = tmp_path / "b.txt"
    p1.write_text("1\n2\n")
    p2.write_text("3\n")
    det.inputfilename = [str(p1)]
    with pytest.raises(ValueError, match="Incomplete final source frame"):
        list(det._source_frames(3))
    p1.write_text("1\n2\n3\n")
    det.inputfilename = [str(p1), str(p2)]
    with pytest.raises(ValueError, match="Incomplete source frame at file boundary"):
        list(det._source_frames(2))

    # Instance-level _source_lines override retains zip_longest fallback.
    xyz = tmp_path / "two.xyz"
    xyz.write_text("2\nf\nC 0 0 0\nC 1.4 0 0\n")
    rng_xyz = ReacNetGenerator(
        inputfilename=str(xyz),
        inputfiletype="xyz",
        atomname=["C"],
        pbc=False,
        use_ase=True,
        bond_hysteresis={"form": 1.2, "break": 1.4},
        output_dir=tmp_path / "xyz",
        timed_output=tmp_path / "xyz.h5",
        nproc=1,
    )
    det_xyz = _Detect.gettype(rng_xyz)
    orig_lines = det_xyz._source_lines
    called_lines = []
    monkeypatch.setattr(
        det_xyz,
        "_source_lines",
        lambda n: called_lines.append(n) or orig_lines(n),
    )
    det_xyz.detect()
    assert called_lines == [4]
    os.remove(rng_xyz.moleculetempfilename)

    # Non-contiguous x/y/z columns and fast-parser malformed fallback paths.
    noncontig = tmp_path / "noncontig.dump"
    noncontig.write_text(
        "\n".join(
            [
                "ITEM: TIMESTEP",
                "0",
                "ITEM: NUMBER OF ATOMS",
                "2",
                "ITEM: BOX BOUNDS pp pp pp",
                "0.0 10.0",
                "0.0 10.0",
                "0.0 10.0",
                "ITEM: ATOMS x id y type z",
                "0.0 1 0.0 1 0.0",
                "1.45 2 0.0 2 0.0",
                "",
            ]
        )
    )
    det_nc = _DetectLAMMPSdump(
        ReacNetGenerator(
            inputfilename=str(noncontig),
            inputfiletype="dump",
            atomname=["C", "O"],
            pbc=False,
            use_ase=True,
            bond_hysteresis={"form": 1.2, "break": 1.4},
        )
    )
    with open(noncontig) as fh:
        stride = det_nc._readNfunc(fh)
    assert stride == 11
    raw_lines = tuple(noncontig.read_text().splitlines(keepends=True))
    assert det_nc._try_fast_parse_frame(0, raw_lines) is not None
    # Malformed timestep, mismatched box bounds, NaN ID, out-of-range ID, duplicate ID.
    bad_ts = ("ITEM: TIMESTEP\n", "bad\n", *raw_lines[2:])
    assert det_nc._try_fast_parse_frame(0, bad_ts) is None
    bad_box = (*raw_lines[:5], "0.0\n", *raw_lines[6:])
    assert det_nc._try_fast_parse_frame(0, bad_box) is None
    nan_id = (*raw_lines[:9], "0.0 nan 0.0 1 0.0\n", raw_lines[10])
    assert det_nc._try_fast_parse_frame(0, nan_id) is None
    oob_id = (*raw_lines[:9], "0.0 0 0.0 1 0.0\n", raw_lines[10])
    assert det_nc._try_fast_parse_frame(0, oob_id) is None
    dup_id = (*raw_lines[:9], "0.0 2 0.0 1 0.0\n", raw_lines[10])
    assert det_nc._try_fast_parse_frame(0, dup_id) is None
    with pytest.raises(ValueError, match="LAMMPS dump file format error"):
        det_nc._readstepfunc((0, ("bad line\n",)))

    # extxyz short frame guard.
    det_ext = _Detect.gettype(
        ReacNetGenerator(
            inputfilename="unused",
            inputfiletype="extxyz",
            atomname=["H"],
        )
    )
    with pytest.raises(RuntimeError, match="No cell information"):
        det_ext._readstepfunc((0, ("1\n",)))


def test_kdtree_hysteresis_exact_ase_threshold_and_strict_search_cutoff(tmp_path):
    """Preserve bit-for-bit ASE distance rounding and strict cutoff semantics."""
    r_c = covalent_radii[6]
    form_mult = 1.25
    break_mult = 1.35
    # Construct a 3D displacement right on the C-C formation threshold where
    # np.sum(dv * dv) <= thr while einsum("ij,ij->i", dv, dv) rounds 1 ULP higher.
    dv = np.array([1.8137306394406951, 0.5308986975642098, -0.19628484525982356])
    thr = 2.0 * r_c * form_mult
    assert float(np.sqrt(np.sum(dv * dv))) <= thr
    assert float(np.sqrt(np.einsum("ij,ij->i", dv[None, :], dv[None, :])[0])) > thr

    for periodic in (False, True):
        traj = tmp_path / f"cc_two_frames_{periodic}.xyz"
        p0 = (0.0, 0.0, 0.0)
        p1 = tuple(float(x) for x in dv)
        frame_txt = (
            f"2\nf\nC {p0[0]!r} {p0[1]!r} {p0[2]!r}\nC {p1[0]!r} {p1[1]!r} {p1[2]!r}\n"
        )
        traj.write_text(frame_txt * 2)
        target = tmp_path / f"cc_{periodic}.h5"
        rng = ReacNetGenerator(
            inputfilename=str(traj),
            inputfiletype="xyz",
            atomname=["C"],
            pbc=periodic,
            cell=np.diag([10.0, 10.0, 10.0]) if periodic else None,
            use_ase=True,
            runHMM=False,
            bond_hysteresis={"form": form_mult, "break": break_mult},
            output_dir=tmp_path / f"out_{periodic}",
            timed_output=target,
            nproc=1,
        )
        rng.run_items(["species", "reactions"])
        validate_timed_output(target)
        spans = list(iter_molecule_ranges(target))
        assert [(s.start_frame, s.end_frame) for s in spans] == [(0, 1)]

    # Strict ASE search cutoff: scaled formation distance equals sum of ASE search radii.
    form_edge = math.nextafter(1.75, math.inf)
    break_edge = math.nextafter(form_edge, math.inf)
    r_h = covalent_radii[1]
    d_exact = (2.0 * r_h) * form_edge
    hh = Atoms("HH", positions=[[0.0, 0.0, 0.0], [d_exact, 0.0, 0.0]], pbc=False)
    ase_i, _, _ = neighbor_list(
        "ijd", hh, np.array([r_h * break_edge, r_h * break_edge])
    )
    assert len(ase_i) == 0

    det_hh = _DetectLAMMPSdump(
        ReacNetGenerator(
            inputfilename="unused",
            inputfiletype="dump",
            atomname=["H"],
            pbc=False,
            use_ase=True,
            bond_hysteresis={"form": form_edge, "break": break_edge},
        )
    )
    det_hh._active_bonds = set()
    assert det_hh._get_hysteresis_bonds(hh, np.diag([20.0, 20.0, 20.0])) == (
        [[], []],
        [[], []],
    )


def test_kdtree_hysteresis_strict_pbc_eligibility_and_active_bonds_sync():
    """Require exact diagonal cells and full boolean PBC for KD-tree; sync _active_bonds."""
    # Nonzero off-diagonal tilt (1e-12) must fall back to ASE.
    tilt_cell = np.array([[10.0, 0.0, 0.0], [1e-12, 10.0, 0.0], [0.0, 0.0, 10.0]])
    hh_tilt = Atoms("HH", positions=[[0.0, 0.0, 0.0], [0.7750000000005, 10.0, 0.0]])
    det_tilt = _DetectLAMMPSdump(
        ReacNetGenerator(
            inputfilename="unused",
            inputfiletype="dump",
            atomname=["H"],
            pbc=True,
            use_ase=True,
            bond_hysteresis={"form": 1.25, "break": 1.35},
        )
    )
    det_tilt._active_bonds = set()
    assert (
        det_tilt._get_kdtree_hysteresis_pairs(
            np.array([1, 1]), hh_tilt.positions, tilt_cell, 1.5
        )
        is None
    )
    assert det_tilt._get_hysteresis_bonds(hh_tilt, tilt_cell) == (
        [[1], [0]],
        [[1], [1]],
    )
    assert det_tilt._active_bonds == {(0, 1)}

    # External _active_bonds reset must invalidate retained-edge reuse.
    hh_hold = Atoms("HH", positions=[[0.0, 0.0, 0.0], [0.81, 0.0, 0.0]])
    det_tilt._active_bonds = set()
    assert det_tilt._get_hysteresis_bonds(hh_hold, np.diag([10.0, 10.0, 10.0])) == (
        [[], []],
        [[], []],
    )

    # Mixed/partial/invalid PBC and wrapped boundary edge cases fall back to ASE cleanly.
    cell = np.diag([10.0, 10.0, 10.0])
    det_bad_pbc = _DetectLAMMPSdump(
        ReacNetGenerator(
            inputfilename="unused",
            inputfiletype="dump",
            atomname=["H"],
            pbc=[True, True],
            use_ase=False,
            bond_hysteresis={"form": 1.25, "break": 1.35},
        )
    )
    assert (
        det_bad_pbc._get_kdtree_hysteresis_pairs(
            np.array([1, 1]), np.array([[0.0, 0.0, 0.2], [0.0, 0.0, 9.8]]), cell, 1.5
        )
        is None
    )

    for pbc_flag, expect_kdtree, expect_bonds in (
        ([True, True, False], False, [[], []]),
        ([False, False, False], True, [[], []]),
        ([True, True, True], True, [[1], [0]]),
    ):
        det_mix = _DetectLAMMPSdump(
            ReacNetGenerator(
                inputfilename="unused",
                inputfiletype="dump",
                atomname=["H"],
                pbc=pbc_flag,
                use_ase=False,
                bond_hysteresis={"form": 1.25, "break": 1.35},
            )
        )
        det_mix._active_bonds = set()
        kdtree_res = det_mix._get_kdtree_hysteresis_pairs(
            np.array([1, 1]), np.array([[0.0, 0.0, 0.2], [0.0, 0.0, 9.8]]), cell, 1.5
        )
        assert (kdtree_res is not None) is expect_kdtree
        pair_atoms = Atoms("HH", positions=[[0.0, 0.0, 0.2], [0.0, 0.0, 9.8]])
        bonds, _ = det_mix._get_hysteresis_bonds(pair_atoms, cell)
        assert bonds == expect_bonds
        # Second identical call hits Open Babel unchanged-graph reuse of retained pairs.
        bonds_repeat, _ = det_mix._get_hysteresis_bonds(pair_atoms, cell)
        assert bonds_repeat == expect_bonds

    # Tiny negative coordinate wrapping to box edge and extreme coordinates fall back to ASE.
    assert (
        det_tilt._get_kdtree_hysteresis_pairs(
            np.array([1, 1]), np.array([[-1e-16, 0.0, 0.0], [0.5, 0.0, 0.0]]), cell, 1.5
        )
        is None
    )
    assert (
        det_tilt._get_kdtree_hysteresis_pairs(
            np.array([1, 1]),
            np.array([[1e8, 0.0, 0.0], [1e8 + 0.5, 0.0, 0.0]]),
            cell,
            1.5,
        )
        is None
    )


def test_fast_dump_parser_rejects_ragged_rows_and_noncanonical_integers(tmp_path):
    """Preserve row boundaries and int()-strict ID/type parsing in LAMMPS dumps."""
    valid_dump = tmp_path / "valid.dump"
    valid_dump.write_text(
        "\n".join(
            [
                "ITEM: TIMESTEP",
                "0",
                "ITEM: NUMBER OF ATOMS",
                "2",
                "ITEM: BOX BOUNDS pp pp pp",
                "0.0 10.0",
                "0.0 10.0",
                "0.0 10.0",
                "ITEM: ATOMS id type x y z",
                "1 1 0.0 0.0 0.0",
                "2 1 1.4 0.0 0.0",
                "",
            ]
        )
    )
    det = _DetectLAMMPSdump(
        ReacNetGenerator(
            inputfilename=str(valid_dump),
            inputfiletype="dump",
            atomname=["C"],
            pbc=False,
            use_ase=True,
            bond_hysteresis={"form": 1.2, "break": 1.4},
        )
    )
    with open(valid_dump) as fh:
        assert det._readNfunc(fh) == 11
    base_lines = tuple(valid_dump.read_text().splitlines(keepends=True))
    # Lines without trailing newline still parse via the fast path.
    no_nl = (*base_lines[:10], base_lines[10].rstrip("\n"))
    assert det._try_fast_parse_frame(0, no_nl) is not None

    # Ragged rows: total float count matches N*n_cols, but shifts trailing '2' as ID 2.
    ragged = (*base_lines[:9], "1 1 0 0 0 2\n", "1 1 0 0\n")
    assert det._try_fast_parse_frame(0, ragged) is None
    with pytest.raises(ValueError, match="duplicate atom ID"):
        det._readstepfunc((0, ragged))

    # Noncanonical integer tokens ('1.0', '1e0') and overflow-sized int strings must fall back.
    for bad_row in ("1.0 1 0 0 0\n", "1 1e0 0 0 0\n"):
        frame = (*base_lines[:9], bad_row, base_lines[10])
        assert det._try_fast_parse_frame(0, frame) is None
        with pytest.raises(ValueError):
            det._readstepfunc((0, frame))

    huge_id = (*base_lines[:9], "100000000000000000000 1 0 0 0\n", base_lines[10])
    assert det._try_fast_parse_frame(0, huge_id) is None
    with pytest.raises(ValueError, match="out of range"):
        det._readstepfunc((0, huge_id))

    # Non-string line element cleanly triggers TypeError fallback in _try_fast_parse_frame.
    assert det._try_fast_parse_frame(0, (*base_lines[:9], None, base_lines[10])) is None
