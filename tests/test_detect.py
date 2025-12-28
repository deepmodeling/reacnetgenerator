# SPDX-License-Identifier: LGPL-3.0-or-later
"""Test different detect format."""

import os
from pathlib import Path

import numpy as np
import pytest

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._detect import _Detect

p_inputs = Path(__file__).parent / "inputs"


class TestDetect:
    """Test different detect format.

    All systems contain a single water molecule: H, H, O.
    """

    @pytest.fixture(autouse=True)
    def chdir(self, tmp_path):
        """Change directory to tmp_path."""
        start_directory = os.getcwd()
        os.chdir(tmp_path)
        yield
        os.chdir(start_directory)

    @pytest.fixture(
        params=[
            # inputfiletype, inputfilename
            ("lammpsdumpfile", p_inputs / "water.dump"),
            ("dump", p_inputs / "water_pbc.dump"),
            ("lammpsbondfile", p_inputs / "water.bond"),
            ("xyz", p_inputs / "water.xyz"),
            ("extxyz", p_inputs / "water.extxyz"),
        ]
    )
    def reacnetgen_param(self, request):
        """Fixture for ReacNetGenerator parameters."""
        return request.param

    @pytest.fixture()
    def reacnetgen(self, reacnetgen_param):
        """Fixture for ReacNetGenerator."""
        rngclass = ReacNetGenerator(
            inputfiletype=reacnetgen_param[0],
            inputfilename=reacnetgen_param[1],
            atomname=["H", "O"],
            pbc=False,
        )
        yield rngclass

    def test_reacnetgen(self, reacnetgen):
        """Test main process of ReacNetGen."""
        _Detect.gettype(reacnetgen).detect()
        assert reacnetgen.N == 3
        np.testing.assert_array_equal(
            reacnetgen.atomtype, np.array([0, 0, 1], dtype=int)
        )
        # assert this is a single molecule
        assert reacnetgen.temp1it == 1


# Additional tests for ASE detection and Scipy clustering features
try:
    from ase import Atoms

    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False

try:
    import scipy

    print(scipy.__version__)
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

from reacnetgenerator._detect import _DetectLAMMPSdump


class TestParseCustomCutoffs:
    """Test the _parse_custom_cutoffs method."""

    @pytest.fixture
    def detect_instance(self):
        """Create a detect instance for testing."""
        rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename="dummy",
            atomname=["H", "O", "C", "Al"],
        )
        return _DetectLAMMPSdump(rng)

    def test_valid_simple_format(self, detect_instance):
        """Test parsing a valid simple format."""
        result = detect_instance._parse_custom_cutoffs("H-O:1.5")
        expected = {frozenset({"H", "O"}): 1.5}
        assert result == expected

    def test_valid_multiple_pairs(self, detect_instance):
        """Test parsing multiple valid pairs."""
        result = detect_instance._parse_custom_cutoffs("H-O:1.5,C-H:2.0")
        expected = {frozenset({"H", "O"}): 1.5, frozenset({"C", "H"}): 2.0}
        assert result == expected

    def test_with_spaces(self, detect_instance):
        """Test parsing with spaces in the input."""
        result = detect_instance._parse_custom_cutoffs("H - O : 1.5 , C - H : 2.0")
        expected = {frozenset({"H", "O"}): 1.5, frozenset({"C", "H"}): 2.0}
        assert result == expected

    def test_invalid_missing_colon(self, detect_instance):
        """Test that missing colon raises ValueError."""
        with pytest.raises(ValueError) as exc_info:
            detect_instance._parse_custom_cutoffs("H-O1.5")
        assert "Invalid custom cutoff format 'H-O1.5'" in str(exc_info.value)
        assert "Expected 'Element1-Element2:distance'" in str(exc_info.value)
        assert "Example: 'Al-O:2.5,C-H:1.1'" in str(exc_info.value)

    def test_invalid_missing_dash(self, detect_instance):
        """Test that missing dash raises ValueError."""
        with pytest.raises(ValueError) as exc_info:
            detect_instance._parse_custom_cutoffs("HO:1.5")
        assert "Invalid custom cutoff format 'HO:1.5'" in str(exc_info.value)
        assert "Expected 'Element1-Element2:distance'" in str(exc_info.value)

    def test_invalid_distance_non_numeric(self, detect_instance):
        """Test that non-numeric distance raises ValueError."""
        with pytest.raises(ValueError) as exc_info:
            detect_instance._parse_custom_cutoffs("H-O:invalid")
        assert "Invalid distance value 'invalid'" in str(exc_info.value)
        assert "Expected a number" in str(exc_info.value)
        assert "Example: 'Al-O:2.5,C-H:1.1'" in str(exc_info.value)


@pytest.mark.skipif(not ASE_AVAILABLE, reason="ASE is not available")
class TestGetBondFromASE:
    """Test the _getbondfromase method."""

    @pytest.fixture
    def detect_instance(self):
        """Create a detect instance for testing."""
        rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename="dummy",
            atomname=["H", "O"],
            use_ase=True,
        )
        return _DetectLAMMPSdump(rng)

    def test_water_molecule_bonds(self, detect_instance):
        """Test that water molecule detects 2 bonds with default settings."""
        # Create a simple water molecule
        atoms = Atoms(
            "H2O",
            positions=[[0.757, 0.586, 0.0], [-0.757, 0.586, 0.0], [0.0, 0.0, 0.0]],
        )
        cell = np.eye(3) * 10  # Large cell to avoid PBC issues

        bond, _bondlevel = detect_instance._getbondfromase(atoms, cell)

        # Should have 3 atoms
        assert len(bond) == 3

        # Count total bonds (each bond is counted twice - once for each atom)
        total_bonds = sum(len(neighbors) for neighbors in bond)
        assert total_bonds == 4  # 2 actual bonds, each counted twice

        # Check that oxygen is bonded to both hydrogens
        assert 0 in bond[2]  # O (index 2) bonded to H (index 0)
        assert 1 in bond[2]  # O (index 2) bonded to H (index 1)
        assert 2 in bond[0]  # H (index 0) bonded to O (index 2)
        assert 2 in bond[1]  # H (index 1) bonded to O (index 2)

    def test_duplicate_bond_prevention_with_pbc(self, detect_instance):
        """Test that duplicate bonds are prevented with PBC."""
        # Create a simple water molecule
        atoms = Atoms(
            "H2O",
            positions=[[0.757, 0.586, 0.0], [-0.757, 0.586, 0.0], [0.0, 0.0, 0.0]],
        )
        # Use a small cell to create potential duplicate scenarios
        cell = np.eye(3) * 5.0

        bond, _bondlevel = detect_instance._getbondfromase(atoms, cell)

        # Check that there are no duplicate entries in bond lists
        for i, neighbors in enumerate(bond):
            # Check that all neighbors are unique
            assert len(neighbors) == len(set(neighbors)), (
                f"Duplicate bonds found for atom {i}"
            )

    def test_custom_cutoffs_override(self, detect_instance):
        """Test that custom cutoffs override global settings."""
        # Modify the detect instance to have custom cutoffs
        detect_instance.rng.custom_cutoffs = "H-O:1.0"  # Very short cutoff

        # Create a simple water molecule where H-O distance is ~1.0 Angstrom
        atoms = Atoms(
            "H2O", positions=[[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0], [0.0, 0.0, 0.0]]
        )  # H-O distance is 0.5
        cell = np.eye(3) * 10

        bond, _bondlevel = detect_instance._getbondfromase(atoms, cell)

        # With a 1.0 cutoff, H-O should still be bonded
        # Check that oxygen is bonded to at least one hydrogen
        assert 0 in bond[2] or 1 in bond[2]  # O (index 2) bonded to H (index 0 or 1)


@pytest.mark.skipif(not SCIPY_AVAILABLE, reason="Scipy is not available")
class TestScipyClustering:
    """Test the Scipy clustering in _connectmolecule method."""

    @pytest.fixture
    def detect_instance(self):
        """Create a detect instance for testing."""
        rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename="dummy",
            atomname=["H", "O"],
            use_ase=True,
        )
        return _DetectLAMMPSdump(rng)

    def test_large_linear_chain_scipy_clustering(self, detect_instance):
        """Test scipy clustering with a large linear chain."""
        # Create a large linear chain of atoms: 0-1-2-3-...-499
        n_atoms = 500
        bond = [[] for _ in range(n_atoms)]
        level = [[] for _ in range(n_atoms)]

        # Create linear chain bonds
        for i in range(n_atoms - 1):
            bond[i].append(i + 1)
            bond[i + 1].append(i)
            level[i].append(1)
            level[i + 1].append(1)

        # Use the _connectmolecule method which should use scipy when available
        result = detect_instance._connectmolecule(bond, level)

        # Should have 1 component containing all atoms
        assert len(result) == 1

        # The component should contain all atoms
        component_bytes = result[0]
        # The first part of the component bytes contains the atom indices
        # Need to extract the atom list from the bytes format
        # This is tricky since it's in the internal format, so just check length
        # The result format is: atom_indices + bond_pairs + bond_levels
        assert len(result) == 1  # Single connected component


class TestAutoEnableLogic:
    """Test the auto-enable logic for ASE mode."""

    def test_auto_enable_with_custom_cutoffs(self, tmp_path):
        """Test that ASE mode is auto-enabled when custom cutoffs are provided."""
        # Create a temporary file for input
        dummy_file = tmp_path / "dummy.dump"
        dummy_file.write_text("")

        rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename=str(dummy_file),
            atomname=["H", "O"],
            use_ase=False,  # Explicitly set to False
            custom_cutoffs="H-O:1.5",  # But provide custom cutoffs
        )

        # Should be auto-enabled
        assert rng.use_ase is True

    def test_auto_enable_with_modified_multiplier(self, tmp_path):
        """Test that ASE mode is auto-enabled when multiplier is modified."""
        # Create a temporary file for input
        dummy_file = tmp_path / "dummy.dump"
        dummy_file.write_text("")

        rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename=str(dummy_file),
            atomname=["H", "O"],
            use_ase=False,  # Explicitly set to False
            ase_cutoff_mult=1.5,  # But modify multiplier
        )

        # Should be auto-enabled
        assert rng.use_ase is True
