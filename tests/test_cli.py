# SPDX-License-Identifier: LGPL-3.0-or-later
"""Tests for the command line interface of the reacnetgenerator package."""

import subprocess
import sys

import reacnetgenerator


def test_module_script():
    """Test the module script functionality by verifying the correct output when executed as a module."""
    expected_version = reacnetgenerator.__version__
    output = subprocess.check_output(
        [sys.executable, "-m", "reacnetgenerator", "--version"]
    ).decode("ascii")
    assert output.splitlines()[0] == f"ReacNetGenerator v{expected_version}"
