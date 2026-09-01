# SPDX-License-Identifier: LGPL-3.0-or-later
"""Test oversized connected-component protection."""

import numpy as np
import pytest
from ase import Atoms

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._detect import _DetectLAMMPSbond, _DetectLAMMPSdump
from reacnetgenerator.commandline import main_parser, parm2cmd


def make_detector(
    detector_class=_DetectLAMMPSdump,
    atom_count=300,
    **kwargs,
):
    """Build a detector with initialized synthetic atom metadata."""
    rng = ReacNetGenerator(
        inputfiletype=(
            "lammpsbondfile"
            if detector_class is _DetectLAMMPSbond
            else "lammpsdumpfile"
        ),
        inputfilename="synthetic",
        atomname=["H"],
        pbc=False,
        **kwargs,
    )
    detector = detector_class(rng)
    detector.N = atom_count
    detector.atomtype = np.zeros(atom_count, dtype=int)
    return detector


def chain_graph(atom_count, component_size):
    """Build one linear component and leave all remaining atoms isolated."""
    bond = [[] for _ in range(atom_count)]
    level = [[] for _ in range(atom_count)]
    for atom in range(component_size - 1):
        bond[atom].append(atom + 1)
        bond[atom + 1].append(atom)
        level[atom].append(1)
        level[atom + 1].append(1)
    return bond, level


def test_default_component_limits():
    """Use the documented atom and fractional limits by default."""
    rng = ReacNetGenerator(
        inputfiletype="lammpsdumpfile",
        inputfilename="synthetic",
        atomname=["H"],
    )

    assert rng.max_component_atoms == 256
    assert rng.max_component_fraction == 0.1
    args = main_parser().parse_args(["-i", "synthetic", "-a", "H"])
    assert args.max_component_atoms == 256
    assert args.max_component_fraction == 0.1


def test_component_at_exact_limit_is_allowed():
    """Allow a component whose size equals the effective limit."""
    detector = make_detector(max_component_fraction=0)
    bond, level = chain_graph(300, 256)

    detector._connectmolecule(bond, level, frame=0, timestep=0)


def test_component_above_atom_limit_fails_with_context():
    """Report actionable context when an ASE component exceeds the limit."""
    detector = make_detector(use_ase=True, max_component_fraction=0)
    bond, level = chain_graph(300, 257)

    with pytest.raises(RuntimeError) as exc_info:
        detector._connectmolecule(bond, level, frame=2, timestep=50)

    message = str(exc_info.value)
    assert "frame 2 (timestep 50)" in message
    assert "257 of 300 atoms" in message
    assert "limit of 256" in message
    assert "H:257" in message
    assert "Bond detection backend: ASE" in message
    assert "--ase-pair-cutoffs El1-El2:0" in message


def test_fractional_limit_scales_with_system_size():
    """Scale the effective limit with the total system size."""
    detector = make_detector(atom_count=3000)
    bond, level = chain_graph(3000, 301)

    with pytest.raises(RuntimeError, match="limit of 300"):
        detector._connectmolecule(bond, level, frame=0, timestep=0)


def test_fractional_limit_rounds_up():
    """Round a fractional component limit up to the next atom."""
    detector = make_detector(atom_count=2561)
    boundary_bond, boundary_level = chain_graph(2561, 257)
    oversized_bond, oversized_level = chain_graph(2561, 258)

    detector._connectmolecule(boundary_bond, boundary_level, frame=0, timestep=0)
    with pytest.raises(RuntimeError, match="limit of 257"):
        detector._connectmolecule(oversized_bond, oversized_level, frame=0, timestep=0)


def test_component_limits_can_be_overridden():
    """Accept public API overrides for both component limits."""
    detector = make_detector(
        max_component_atoms=400,
        max_component_fraction=0,
    )
    bond, level = chain_graph(300, 300)

    detector._connectmolecule(bond, level, frame=0, timestep=0)


def test_zero_atom_limit_disables_guard():
    """Disable the guard when the atom limit is explicitly zero."""
    detector = make_detector(
        max_component_atoms=0,
        max_component_fraction=1,
    )
    bond, level = chain_graph(300, 300)

    detector._connectmolecule(bond, level, frame=0, timestep=0)


def test_empty_component_list_is_allowed():
    """Allow an empty component collection without raising an error."""
    detector = make_detector(max_component_fraction=0)

    detector._validate_component_sizes([], frame=0, timestep=0)


@pytest.mark.parametrize(
    ("key", "value", "message"),
    [
        ("max_component_atoms", -1, "non-negative integer"),
        ("max_component_atoms", 1.5, "non-negative integer"),
        ("max_component_atoms", True, "non-negative integer"),
        ("max_component_fraction", -0.1, "between 0 and 1"),
        ("max_component_fraction", 1.1, "between 0 and 1"),
        ("max_component_fraction", True, "between 0 and 1"),
    ],
)
def test_python_api_rejects_invalid_component_limits(key, value, message):
    """Reject invalid component limits passed through the Python API."""
    with pytest.raises(ValueError, match=message):
        ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename="synthetic",
            atomname=["H"],
            **{key: value},
        )


@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("--max-component-atoms", "-1"),
        ("--max-component-atoms", "1.5"),
        ("--max-component-fraction", "-0.1"),
        ("--max-component-fraction", "1.1"),
    ],
)
def test_cli_rejects_invalid_component_limits(option, value):
    """Reject invalid component limits passed through the CLI."""
    parser = main_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["-i", "synthetic", "-a", "H", option, value])


def test_parm2cmd_preserves_zero_component_limits():
    """Preserve explicit zero values when rebuilding a CLI command."""
    command = parm2cmd(
        {
            "inputfilename": "synthetic",
            "inputfiletype": "lammpsdumpfile",
            "atomname": ["H"],
            "max_component_atoms": 0,
            "max_component_fraction": 0,
        }
    )

    assert command[command.index("--max-component-atoms") + 1] == "0"
    assert command[command.index("--max-component-fraction") + 1] == "0"
    args = main_parser().parse_args(command[1:])
    assert args.max_component_atoms == 0
    assert args.max_component_fraction == 0


def test_openbabel_error_recommends_ase_controls():
    """Recommend ASE controls for oversized coordinate-based components."""
    detector = make_detector(use_ase=False, max_component_fraction=0)
    bond, level = chain_graph(300, 257)

    with pytest.raises(RuntimeError, match="consider --use-ase"):
        detector._connectmolecule(bond, level, frame=0, timestep=0)


def test_explicit_bond_error_recommends_upstream_cutoff():
    """Recommend upstream bond controls for explicit bond trajectories."""
    detector = make_detector(
        detector_class=_DetectLAMMPSbond,
        max_component_fraction=0,
    )
    bond, level = chain_graph(300, 257)

    with pytest.raises(RuntimeError, match="upstream bond-order cutoff"):
        detector._connectmolecule(bond, level, frame=0, timestep=0)


def test_detect_stops_before_writing_molecule_data(tmp_path):
    """Propagate the guard through trajectory detection before serialization."""
    atom_count = 300
    trajectory = tmp_path / "synthetic.dump"
    atom_lines = "\n".join(
        f"{atom + 1} 1 {0.6 * atom} 0.0 0.0" for atom in range(atom_count)
    )
    trajectory.write_text(
        f"""ITEM: TIMESTEP
50
ITEM: NUMBER OF ATOMS
{atom_count}
ITEM: BOX BOUNDS pp pp pp
0.0 200.0
0.0 10.0
0.0 10.0
ITEM: ATOMS id type x y z
{atom_lines}
"""
    )
    rng = ReacNetGenerator(
        inputfiletype="lammpsdumpfile",
        inputfilename=trajectory,
        atomname=["H"],
        pbc=False,
        use_ase=True,
        nproc=1,
    )
    detector = _DetectLAMMPSdump(rng)

    with pytest.raises(RuntimeError, match=r"frame 0 \(timestep 50\)"):
        detector.detect()
    assert detector.moleculetempfilename is None


def test_ase_pair_cutoff_can_remove_oversized_component():
    """Allow an explicit ASE pair cutoff to remove an unintended component."""
    atom_count = 300
    atoms = Atoms(
        ["H"] * atom_count,
        positions=[(0.6 * atom, 0.0, 0.0) for atom in range(atom_count)],
        cell=[200.0, 10.0, 10.0],
        pbc=False,
    )

    detector = make_detector(use_ase=True)
    bond, level = detector._getbondfromase(atoms, atoms.cell.array)
    with pytest.raises(RuntimeError, match="ASE natural cutoffs"):
        detector._connectmolecule(bond, level, frame=0, timestep=0)

    detector = make_detector(use_ase=True, custom_cutoffs="H-H:0")
    bond, level = detector._getbondfromase(atoms, atoms.cell.array)
    molecules = detector._connectmolecule(bond, level, frame=0, timestep=0)
    assert len(molecules) == atom_count
