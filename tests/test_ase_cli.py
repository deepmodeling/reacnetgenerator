# SPDX-License-Identifier: LGPL-3.0-or-later
"""Test ASE mode command line interface."""

from pathlib import Path

from reacnetgenerator.commandline import main_parser

p_inputs = Path(__file__).parent / "inputs"


class TestASECommandLine:
    """Test ASE mode command line interface."""

    def test_ase_cli_arguments_exist(self):
        """Test that ASE-related command line arguments exist."""
        parser = main_parser()

        # Parse an example command with ASE args
        args = parser.parse_args(
            [
                "-i",
                str(p_inputs / "water.dump"),
                "-a",
                "H",
                "O",
                "--use-ase",
                "--ase-cutoff-mult",
                "1.5",
                "--ase-pair-cutoffs",
                "H-O:1.2,C-H:2.0",
            ]
        )

        # Check that the arguments are parsed correctly
        assert args.use_ase is True
        assert args.ase_cutoff_mult == 1.5
        assert args.ase_pair_cutoffs == "H-O:1.2,C-H:2.0"

    def test_ase_cli_default_values(self):
        """Test that ASE-related command line arguments have correct defaults."""
        parser = main_parser()

        # Parse a minimal command without ASE args
        args = parser.parse_args(["-i", str(p_inputs / "water.dump"), "-a", "H", "O"])

        # Check default values
        assert args.use_ase is False
        assert args.ase_cutoff_mult == 1.2
        assert args.ase_pair_cutoffs is None

    def test_ase_cli_help_contains_ase_args(self):
        """Test that ASE-related arguments exist in the parser."""
        parser = main_parser()

        # Check specific arguments exist in the parser
        actions = {action.dest: action for action in parser._actions}
        assert "use_ase" in actions
        assert "ase_cutoff_mult" in actions
        assert "ase_pair_cutoffs" in actions

        # Check that the help text contains expected content
        help_text = parser.format_help()
        assert "--use-ase" in help_text
        assert "--ase-cutoff-mult" in help_text
        assert "--ase-pair-cutoffs" in help_text
