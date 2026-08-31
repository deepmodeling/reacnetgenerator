# SPDX-License-Identifier: LGPL-3.0-or-later
"""Test different detect format."""

import os
from pathlib import Path

import numpy as np
import pytest

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator import _detect as detect_module
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

    def test_extxyz_uses_each_frames_lattice(self, tmp_path, mocker):
        """Variable-cell extxyz frames should use their own lattice."""
        rng = ReacNetGenerator(
            inputfiletype="extxyz",
            inputfilename=tmp_path / "variable-cell.extxyz",
            atomname=["H"],
            pbc=True,
        )
        detector = _Detect.gettype(rng)
        get_bonds = mocker.patch.object(
            detector, "_getbondfromcrd", return_value=([[]], [[]])
        )
        mocker.patch.object(detector, "_connectmolecule", return_value=[])
        frame_cells = [np.eye(3) * 5.0, np.eye(3) * 10.0]

        for step, cell in enumerate(frame_cells):
            lattice = " ".join(str(value) for value in cell.ravel())
            lines = (
                "1\n",
                f'Lattice="{lattice}" Properties=species:S:1:pos:R:3\n',
                "H 0.0 0.0 0.0\n",
            )
            detector._readstepfunc((step, lines))

        for call, expected_cell in zip(
            get_bonds.call_args_list, frame_cells, strict=True
        ):
            np.testing.assert_array_equal(call.args[1], expected_cell)


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


@pytest.mark.skipif(not ASE_AVAILABLE, reason="ASE is not available")
class TestPeriodicOpenBabelCandidates:
    """Keep periodic cKDTree candidates equivalent to Open Babel."""

    @pytest.fixture
    def detect_instance(self):
        """Create a periodic detector that retains Open Babel bond orders."""
        rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename="dummy",
            atomname=["H", "C", "N", "O", "F", "P", "Cl"],
            pbc=True,
            use_ase=False,
        )
        return _DetectLAMMPSdump(rng)

    @staticmethod
    def _molecule_records(detect_instance, result):
        assert result is not None
        return detect_instance._connectmolecule(*result)

    def _assert_matches_reference(self, detect_instance, atoms, cell):
        """Return the accelerated result after checking molecule semantics."""
        reference = detect_instance._getbondfromopenbabel(atoms, cell)
        accelerated = detect_instance._getbondfromperiodicneighborlist(atoms, cell)

        assert accelerated is not None
        assert self._molecule_records(
            detect_instance, accelerated
        ) == self._molecule_records(detect_instance, reference)
        return accelerated

    def test_matches_reference_across_periodic_boundary(self, detect_instance):
        """Find neighbors on opposite faces without changing molecule records."""
        cell = np.eye(3) * 12.0
        atoms = Atoms(
            "CHOH",
            positions=[
                [0.2, 1.0, 1.0],
                [11.3, 1.0, 1.01],
                [6.0, 6.0, 6.0],
                [6.8, 6.0, 6.01],
            ],
            cell=cell,
            pbc=True,
        )

        accelerated = self._assert_matches_reference(detect_instance, atoms, cell)
        assert 1 in accelerated[0][0]

    def test_matches_reference_when_openbabel_removes_excess_bonds(
        self,
        detect_instance,
    ):
        """Reuse Open Babel's valence and angle cleanup after candidate search."""
        cell = np.eye(3) * 12.0
        atoms = Atoms(
            "CHHHHH",
            positions=np.array(
                [
                    [6.00, 6.00, 6.000],
                    [6.90, 6.00, 6.001],
                    [5.05, 6.00, 6.002],
                    [6.00, 7.00, 6.003],
                    [6.00, 4.95, 6.004],
                    [6.00, 6.00, 7.080],
                ]
            ),
            cell=cell,
            pbc=True,
        )

        self._assert_matches_reference(detect_instance, atoms, cell)

    def test_matches_reference_for_equal_z_coordinates(self, detect_instance):
        """Keep rounded equal-z coordinates on the accelerated path."""
        cell = np.eye(3) * 12.0
        atoms = Atoms(
            "CHHHHH",
            positions=np.array(
                [
                    [6.00, 6.00, 6.00],
                    [6.90, 6.00, 6.00],
                    [5.05, 6.00, 6.00],
                    [6.00, 7.00, 6.00],
                    [6.00, 4.95, 6.00],
                    [6.00, 6.00, 7.08],
                ]
            ),
            cell=cell,
            pbc=True,
        )

        self._assert_matches_reference(detect_instance, atoms, cell)

    def test_matches_reference_for_phosphorus_sixth_bond_rule(
        self,
        detect_instance,
    ):
        """Preserve Open Babel's element-specific phosphorus insertion rule."""
        cell = np.eye(3) * 12.0
        angles = np.arange(6) * np.pi / 3.0
        positions = np.concatenate(
            (
                [[6.0, 6.0, 6.0]],
                np.column_stack(
                    (
                        6.0 + 1.6 * np.cos(angles),
                        6.0 + 1.6 * np.sin(angles),
                        6.01 + np.arange(6) * 0.01,
                    )
                ),
            )
        )
        atoms = Atoms("PFFFFFH", positions=positions, cell=cell, pbc=True)

        accelerated = self._assert_matches_reference(detect_instance, atoms, cell)
        assert accelerated[0][0] == [1, 2, 3, 4, 5]
        assert accelerated[0][6] == []

    @pytest.mark.parametrize("seed", range(3))
    def test_matches_reference_for_seeded_mixed_elements(
        self,
        detect_instance,
        seed,
    ):
        """Differentially cover insertion, cleanup, and bond-order perception."""
        random = np.random.default_rng(seed)
        atomic_numbers = random.choice(
            np.array([1, 6, 7, 8, 9, 15, 17]),
            size=96,
        )
        positions = random.uniform(0.0, 18.0, size=(96, 3))
        cell = np.eye(3) * 18.0
        atoms = Atoms(
            numbers=atomic_numbers,
            positions=positions,
            cell=cell,
            pbc=True,
        )

        self._assert_matches_reference(detect_instance, atoms, cell)

    @pytest.mark.parametrize(
        ("atoms", "cell"),
        [
            (
                Atoms("HH", positions=[[0, 0, 0], [0, 0, 1.07]]),
                np.eye(3) * 10.0,
            ),
            (
                Atoms("CH", positions=[[0, 0, 0], [1, 0, 0.1]]),
                np.eye(3) * 2.0,
            ),
            (
                Atoms("HH", positions=[[0, 0, 0], [np.nan, 0, 1.0]]),
                np.eye(3) * 10.0,
            ),
            (
                Atoms(numbers=[0, 1], positions=[[0, 0, 0], [0, 0, 1.0]]),
                np.eye(3) * 10.0,
            ),
        ],
        ids=[
            "cutoff-boundary",
            "narrow-cell",
            "invalid-coordinate",
            "invalid-radius",
        ],
    )
    def test_falls_back_for_ambiguous_or_unsupported_geometry(
        self,
        detect_instance,
        atoms,
        cell,
    ):
        """Leave unsafe geometries to the historical Open Babel path."""
        atoms.set_cell(cell)
        atoms.set_pbc(True)

        assert detect_instance._getbondfromperiodicneighborlist(atoms, cell) is None

    @pytest.mark.parametrize(
        "cell",
        [
            np.array(
                [
                    [10.0, 0.0, 0.0],
                    [2.0, 9.0, 0.0],
                    [1.0, 1.5, 8.0],
                ]
            ),
            np.array(
                [
                    [0.0, 10.0, 0.0],
                    [-9.0, 0.0, 0.0],
                    [0.0, 0.0, 8.0],
                ]
            ),
            np.diag([10.0, 0.0, 10.0]),
            np.diag([10.0, -9.0, 8.0]),
            np.array(
                [
                    [np.nan, 0.0, 0.0],
                    [0.0, 9.0, 0.0],
                    [0.0, 0.0, 8.0],
                ]
            ),
        ],
        ids=["triclinic", "rotated", "zero-length", "negative-length", "nonfinite"],
    )
    def test_falls_back_for_unsupported_cells(self, detect_instance, cell):
        """Leave unsupported cells to the reference implementation."""
        atoms = Atoms("CHOH", positions=np.zeros((4, 3)), pbc=True)

        assert detect_instance._getbondfromperiodicneighborlist(atoms, cell) is None

    def test_dispatches_only_safe_large_periodic_frames(
        self,
        detect_instance,
        monkeypatch,
    ):
        """Use cKDTree only above its crossover and preserve every fallback."""
        atoms = Atoms("HH", positions=[[0, 0, 0], [0, 0, 1]], cell=np.eye(3) * 10)
        accelerated = ([[1], [0]], [[1], [1]])
        reference = ([[], []], [[], []])
        calls = []
        monkeypatch.setattr(
            detect_module,
            "_OPENBABEL_PERIODIC_NEIGHBOR_MIN_ATOMS",
            2,
        )
        monkeypatch.setattr(
            detect_instance,
            "_getbondfromperiodicneighborlist",
            lambda *args: calls.append("fast") or accelerated,
        )
        monkeypatch.setattr(
            detect_instance,
            "_getbondfromopenbabel",
            lambda *args: calls.append("reference") or reference,
        )

        assert detect_instance._getbondfromcrd(atoms, np.eye(3) * 10) == accelerated
        assert calls == ["fast"]

        calls.clear()
        monkeypatch.setattr(
            detect_instance,
            "_getbondfromperiodicneighborlist",
            lambda *args: calls.append("fast") or None,
        )
        assert detect_instance._getbondfromcrd(atoms, np.eye(3) * 10) == reference
        assert calls == ["fast", "reference"]

        calls.clear()
        monkeypatch.setattr(
            detect_module,
            "_OPENBABEL_PERIODIC_NEIGHBOR_MIN_ATOMS",
            3,
        )
        assert detect_instance._getbondfromcrd(atoms, np.eye(3) * 10) == reference
        assert calls == ["reference"]

        calls.clear()
        monkeypatch.setattr(
            detect_module,
            "_OPENBABEL_PERIODIC_NEIGHBOR_MIN_ATOMS",
            2,
        )
        triclinic_cell = np.array(
            [[10.0, 0.0, 0.0], [1.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
        )
        assert detect_instance._getbondfromcrd(atoms, triclinic_cell) == reference
        assert calls == ["reference"]

        calls.clear()
        detect_instance.pbc = False
        assert detect_instance._getbondfromcrd(atoms, np.eye(3) * 10) == reference
        assert calls == ["reference"]

        calls.clear()
        detect_instance.pbc = True
        detect_instance.use_ase = True
        monkeypatch.setattr(
            detect_instance,
            "_getbondfromase",
            lambda *args: calls.append("ase") or accelerated,
        )
        assert detect_instance._getbondfromcrd(atoms, np.eye(3) * 10) == accelerated
        assert calls == ["ase"]


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

    def test_custom_cutoffs_change_bond_graph(self):
        """A pair cutoff should detect a bond excluded by natural cutoffs."""
        default_rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename="dummy",
            atomname=["H", "O"],
            use_ase=True,
        )
        custom_rng = ReacNetGenerator(
            inputfiletype="lammpsdumpfile",
            inputfilename="dummy",
            atomname=["H", "O"],
            use_ase=True,
            custom_cutoffs="H-O:1.5",
        )
        default_detector = _DetectLAMMPSdump(default_rng)
        custom_detector = _DetectLAMMPSdump(custom_rng)
        # At 1.3 A, H-O is outside the default natural cutoff but inside the
        # configured pair cutoff.
        atoms = Atoms("HO", positions=[[0.0, 0.0, 0.0], [1.3, 0.0, 0.0]])
        cell = np.eye(3) * 10

        default_bonds, _ = default_detector._getbondfromase(atoms.copy(), cell)
        custom_bonds, _ = custom_detector._getbondfromase(atoms.copy(), cell)

        assert default_bonds == [[], []]
        assert custom_bonds == [[1], [0]]


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
