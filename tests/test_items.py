# SPDX-License-Identifier: LGPL-3.0-or-later
"""Tests for --items selective execution."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from reacnetgenerator import ReacNetGenerator

INPUT_NAME = "water.bond"
INPUT_SOURCE = Path(__file__).parent / "inputs" / INPUT_NAME


@pytest.fixture(autouse=True)
def local_input(tmp_path, monkeypatch):
    """Run every case against a checked-in trajectory in an isolated directory."""
    shutil.copyfile(INPUT_SOURCE, tmp_path / INPUT_NAME)
    monkeypatch.chdir(tmp_path)


@pytest.fixture()
def rng():
    """Create a generator without downloading any external fixture."""
    return ReacNetGenerator(
        inputfilename=INPUT_NAME,
        atomname=["H", "O"],
        inputfiletype="lammpsbondfile",
        runHMM=False,
        SMILES=False,
        nproc=1,
    )


class TestRunItems:
    """Test the run_items() method and --items CLI."""

    def test_species_only(self, rng):
        """Species-only mode produces .species file without PATH step."""
        rng.run_items(["species"])
        assert os.path.exists(rng.speciesfilename)
        # Reaction file should NOT exist — PATH was skipped
        assert not os.path.exists(rng.reactionfilename)

    def test_species_reactions(self, rng):
        """species+reactions produces both .species and .reaction files."""
        rng.run_items(["species", "reactions"])
        assert os.path.exists(rng.speciesfilename)
        assert os.path.exists(rng.reactionfilename)
        assert rng.mname.dtype == object

    def test_smiles_molecule_names_use_object_dtype(self):
        """SMILES names should not expand to the longest fixed-width string."""
        rng = ReacNetGenerator(
            inputfilename=INPUT_NAME,
            atomname=["H", "O"],
            inputfiletype="lammpsbondfile",
            runHMM=False,
            SMILES=True,
            nproc=1,
        )
        rng.run_items(["species", "reactions"])
        assert rng.mname.dtype == object

    def test_explicit_species_overrides_legacy_flag_for_one_run(self, rng):
        """Honor a requested species file and then restore caller configuration."""
        rng.needprintspecies = False

        rng.run_items(["species", "reactions"])

        assert os.path.exists(rng.speciesfilename)
        assert rng.needprintspecies is False

    def test_all_items_matches_runanddraw(self, rng):
        """Passing all items is equivalent to runanddraw()."""
        rng.run_items(["species", "reactions", "network", "report"])
        assert os.path.exists(rng.speciesfilename)
        assert os.path.exists(rng.reactionfilename)
        assert os.path.exists(rng.resultfilename)

    def test_none_defaults_to_all(self, rng):
        """run_items(None) produces the full analysis."""
        rng.run_items(None)
        assert os.path.exists(rng.speciesfilename)
        assert os.path.exists(rng.reactionfilename)
        assert os.path.exists(rng.resultfilename)

    def test_invalid_item_raises(self, rng):
        """Unknown item name raises ValueError."""
        with pytest.raises(ValueError, match="Unknown item"):
            rng.run_items(["nonexistent"])

    def test_cli_items_species(self):
        """CLI --items species produces .species file."""
        subprocess.check_call(
            [
                sys.executable,
                "-m",
                "reacnetgenerator",
                "-i",
                INPUT_NAME,
                "-a",
                "H",
                "O",
                "--nohmm",
                "--items",
                "species",
            ]
        )
        assert os.path.exists(f"{INPUT_NAME}.species")
