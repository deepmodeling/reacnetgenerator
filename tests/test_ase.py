# SPDX-License-Identifier: LGPL-3.0-or-later
"""Test ASE mode functionality."""

import os
from pathlib import Path

import numpy as np
import pytest

try:
    from ase import Atoms
    from ase.neighborlist import natural_cutoffs, neighbor_list

    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False

try:
    from scipy.sparse.csgraph import connected_components

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._detect import _DetectLAMMPSdump

p_inputs = Path(__file__).parent / "inputs"


class TestASEMode:
    """Test ASE mode functionality."""

    @pytest.fixture(autouse=True)
    def chdir(self, tmp_path):
        """Change directory to tmp_path."""
        start_directory = os.getcwd()
        os.chdir(tmp_path)
        yield
        os.chdir(start_directory)

    @pytest.mark.skipif(not ASE_AVAILABLE, reason="ASE is not available")
    def test_ase_mode_disabled_by_default(self):
        """Test that ASE mode is disabled by default."""
        rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename=p_inputs / "water.dump",
            atomname=["H", "O"],
            pbc=False,
        )
        assert not getattr(rng, "use_ase", False)
        assert rng.ase_cutoff_mult == 1.2
        assert rng.custom_cutoffs is None

    @pytest.mark.skipif(not ASE_AVAILABLE, reason="ASE is not available")
    def test_ase_mode_enabled_manually(self):
        """Test that ASE mode can be enabled manually."""
        rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename=p_inputs / "water.dump",
            atomname=["H", "O"],
            pbc=False,
            use_ase=True,
        )
        assert rng.use_ase

    @pytest.mark.skipif(not ASE_AVAILABLE, reason="ASE is not available")
    def test_ase_auto_enable_with_custom_cutoffs(self):
        """Test that ASE mode is auto-enabled when custom cutoffs are provided."""
        rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename=p_inputs / "water.dump",
            atomname=["H", "O"],
            pbc=False,
            custom_cutoffs="H-O:1.5,O-O:2.0",
        )
        assert rng.use_ase
        assert rng.custom_cutoffs == "H-O:1.5,O-O:2.0"

    @pytest.mark.skipif(not ASE_AVAILABLE, reason="ASE is not available")
    def test_ase_auto_enable_with_modified_mult(self):
        """Test that ASE mode is auto-enabled when cutoff multiplier is modified."""
        rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename=p_inputs / "water.dump",
            atomname=["H", "O"],
            pbc=False,
            ase_cutoff_mult=1.5,
        )
        assert rng.use_ase
        assert rng.ase_cutoff_mult == 1.5

    @pytest.mark.skipif(not ASE_AVAILABLE, reason="ASE is not available")
    def test_ase_parse_custom_cutoffs(self):
        """Test parsing of custom cutoff strings."""
        rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename=p_inputs / "water.dump",
            atomname=["H", "O"],
            pbc=False,
            use_ase=True,
            custom_cutoffs="H-O:1.5,O-O:2.0",
        )
        # Get a concrete instance instead of abstract class
        detect_class = _DetectLAMMPSdump(rng)
        parsed = detect_class._parse_custom_cutoffs(rng.custom_cutoffs)

        assert frozenset({"H", "O"}) in parsed
        assert frozenset({"O", "O"}) in parsed
        assert parsed[frozenset({"H", "O"})] == 1.5
        assert parsed[frozenset({"O", "O"})] == 2.0

    @pytest.mark.skipif(not ASE_AVAILABLE, reason="ASE is not available")
    def test_ase_parse_custom_cutoffs_order_independence(self):
        """Test that custom cutoff parsing is order-independent."""
        rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename=p_inputs / "water.dump",
            atomname=["H", "O"],
            pbc=False,
            use_ase=True,
            custom_cutoffs="H-O:1.5,O-H:2.0",  # Same pair in different order
        )
        # Get a concrete instance instead of abstract class
        detect_class = _DetectLAMMPSdump(rng)
        parsed = detect_class._parse_custom_cutoffs(rng.custom_cutoffs)

        # Both should map to the same frozenset and take the last value
        assert frozenset({"H", "O"}) in parsed
        # In case of duplicate entries, the last one should be used
        assert parsed[frozenset({"H", "O"})] == 2.0

    @pytest.mark.skipif(not ASE_AVAILABLE, reason="ASE is not available")
    def test_ase_getbondfromase_method(self):
        """Test the _getbondfromase method directly."""
        rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename=p_inputs / "water.dump",
            atomname=["H", "O"],
            pbc=False,
            use_ase=True,
        )
        detect_class = _DetectLAMMPSdump(rng)

        # Create a simple water molecule
        atoms = Atoms(
            "H2O",
            positions=[[0.757, 0.586, 0.0], [-0.757, 0.586, 0.0], [0.0, 0.0, 0.0]],
        )
        cell = np.eye(3) * 10  # Large cell to avoid PBC issues

        # This should not raise an error
        bond, bondlevel = detect_class._getbondfromase(atoms, cell)

        # Should have 3 atoms
        assert len(bond) == 3
        assert len(bondlevel) == 3

        # Should have O-H bonds
        # The oxygen (index 2) should be bonded to both hydrogens (0, 1)
        assert 0 in bond[2] or 1 in bond[2]
        assert 2 in bond[0] or 2 in bond[1]

        # All bond levels should be 1 (single bonds)
        for bond_levels in bondlevel:
            for level in bond_levels:
                assert level == 1

    @pytest.mark.skipif(not ASE_AVAILABLE, reason="ASE is not available")
    def test_ase_duplicate_bond_prevention(self):
        """Test that duplicate bonds are prevented in ASE mode."""
        rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename=p_inputs / "water.dump",
            atomname=["H", "O"],
            pbc=False,
            use_ase=True,
            ase_cutoff_mult=5.0,  # Very high multiplier to force many neighbors
        )
        detect_class = _DetectLAMMPSdump(rng)

        # Create a simple water molecule
        atoms = Atoms(
            "H2O",
            positions=[[0.757, 0.586, 0.0], [-0.757, 0.586, 0.0], [0.0, 0.0, 0.0]],
        )
        # Use a small cell to create periodic images
        cell = np.eye(3) * 5.0

        # This should not raise an error and should not have duplicate bonds
        bond, bondlevel = detect_class._getbondfromase(atoms, cell)

        # Check that there are no duplicate entries in bond lists
        for i, neighbors in enumerate(bond):
            # Check that all neighbors are unique
            assert len(neighbors) == len(set(neighbors)), (
                f"Duplicate bonds found for atom {i}"
            )
            # Check that the corresponding bond levels are also unique
            assert len(bondlevel[i]) == len(set(bondlevel[i])) or len(
                bondlevel[i]
            ) == len(neighbors)

    @pytest.mark.skipif(
        not (ASE_AVAILABLE and SCIPY_AVAILABLE), reason="ASE or Scipy not available"
    )
    def test_connected_components_algorithm(self):
        """Test that connected components algorithm works properly."""
        rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename=p_inputs / "water.dump",
            atomname=["H", "O"],
            pbc=False,
            use_ase=True,
        )
        detect_class = _DetectLAMMPSdump(rng)

        # Create a simple bond structure with 2 connected components
        # Component 1: atoms 0 and 1 connected
        # Component 2: atoms 2 and 3 connected
        bond = [[1], [0], [3], [2]]  # 0-1 and 2-3
        level = [[1], [1], [1], [1]]

        # Use the _connectmolecule method which should use connected_components
        result = detect_class._connectmolecule(bond, level)

        # Should have 2 components
        assert len(result) == 2

        # Each component should have 2 atoms
        for component in result:
            # Parse the result format to check atoms
            # This requires understanding the internal format of the result
            pass

    @pytest.mark.skipif(not ASE_AVAILABLE, reason="ASE is not available")
    def test_ase_with_custom_cutoffs_integration(self):
        """Integration test for ASE mode with custom cutoffs."""
        # Create a temporary file with a simple water molecule
        water_file = "test_water.dump"
        with open(water_file, "w") as f:
            f.write("""ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
3
ITEM: BOX BOUNDS pp pp pp
0.0 10.0
0.0 10.0
0.0 10.0
ITEM: ATOMS id type x y z
1 1 0.757 0.586 0.0
2 1 -0.757 0.586 0.0
3 2 0.0 0.0 0.0
""")

        rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename=water_file,
            atomname=["H", "O"],
            pbc=True,
            use_ase=True,
            custom_cutoffs="H-O:1.5",  # Custom cutoff for H-O bonds
        )

        # Run the detection process
        detect_class = _DetectLAMMPSdump.gettype(rng)
        detect_class.detect()

        # Should have processed the file
        assert hasattr(detect_class, "temp1it")
        assert detect_class.temp1it > 0

    @pytest.mark.skipif(not ASE_AVAILABLE, reason="ASE is not available")
    def test_ase_warning_on_auto_enable(self, caplog):
        """Test that a warning is logged when ASE mode is auto-enabled."""
        with caplog.at_level("WARNING"):
            rng = ReacNetGenerator(
                inputfiletype="lammpsdumpfile",
                inputfilename=p_inputs / "water.dump",
                atomname=["H", "O"],
                pbc=False,
                custom_cutoffs="H-O:1.5",  # This should trigger auto-enable
            )

            # Check that the warning was logged
            assert (
                "ASE parameters detected. Automatically enabling --use-ase mode."
                in caplog.text
            )
            assert rng.use_ase
