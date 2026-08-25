# SPDX-License-Identifier: LGPL-3.0-or-later
"""Tests for --items selective execution."""

import json
import os
import subprocess
import sys

import pytest

from reacnetgenerator import ReacNetGenerator

with open(os.path.join(os.path.dirname(__file__), "test.json")) as f:
    test_data = json.load(f)


@pytest.fixture(autouse=True)
def chdir(tmp_path):
    """Change directory to tmp_path."""
    start = os.getcwd()
    os.chdir(tmp_path)
    yield
    os.chdir(start)


@pytest.fixture()
def rng():
    """Create a ReacNetGenerator instance with the first test param set."""
    return ReacNetGenerator(**test_data[0]["rngparams"])


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
                test_data[0]["rngparams"]["inputfilename"],
                "-a",
                *test_data[0]["rngparams"]["atomname"],
                "--nohmm",
                "--items",
                "species",
                "--urls",
                test_data[0]["rngparams"]["urls"][0]["fn"],
                test_data[0]["rngparams"]["urls"][0]["url"][0],
            ]
        )
        assert os.path.exists(
            f"{test_data[0]['rngparams']['inputfilename']}.species"
        )
