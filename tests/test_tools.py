# SPDX-License-Identifier: LGPL-3.0-or-later
"""Test the tools module."""
from pathlib import Path

import numpy as np

from reacnetgenerator.tools import calculate_rate, read_species

p_cwd = Path(__file__).parent


def test_species():
    """Test the species reading."""
    step_idx, n_species = read_species(p_cwd / "methane.species")
    np.testing.assert_array_equal(
        step_idx, np.array([0, 1000, 2000, 3000, 4000, 5000, 6000, 7000], dtype=int)
    )
    np.testing.assert_array_equal(
        n_species["[H]OO"], np.array([0, 0, 0, 1, 1, 1, 1, 1], dtype=int)
    )


def test_rates():
    """Test the rate calculation."""
    rates = calculate_rate(
        p_cwd / "methane.species",
        p_cwd / "methane.reactionabcd",
        np.eye(3) * 3.7601e1,
        0.1,
    )
    assert "OO->O=O" in rates
