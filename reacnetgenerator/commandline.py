# SPDX-License-Identifier: LGPL-3.0-or-later
# cython: language_level=3
# cython: linetrace=True
"""Command line interface for reacnetgenerator."""
import argparse
import textwrap
from typing import List

from . import __version__
from ._detect import _Detect


def main_parser() -> argparse.ArgumentParser:
    """Return main parser.

    Returns
    -------
    argparse.ArgumentParser
        reacnetgenerator cli parser
    """
    parser = argparse.ArgumentParser(
        description="ReacNetGenerator is an automatic reaction network generator for reactive molecular dynamics simulation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=textwrap.dedent(
            """\
        Examples:
            reacnetgenerator --type bond -i bonds.reaxc -a C H O --nohmm
        """
        ),
    )
    parser.add_argument(
        "-i",
        "--inputfilename",
        nargs="*",
        help="Input trajectory file(s), e.g. bonds.reaxc; If multiple file are given, they are concatenated.",
        required=True,
    )
    parser.add_argument(
        "-a",
        "--atomname",
        help="Atomic names in the trajectory, e.g. C H O",
        nargs="*",
        required=True,
    )
    parser.add_argument(
        "--nohmm",
        help=(
            "Process trajectory without Hidden Markov Model (HMM). If one wants to enable HMM, firstly "
            "read the related section in the paper"
        ),
        action="store_true",
    )
    parser.add_argument("--miso", help="Merge the isomers", type=int, default=0)
    parser_type = parser.add_mutually_exclusive_group()
    parser_type.add_argument(
        "--dump",
        help="(deprecated) Process the LAMMPS dump file. Please use `--type dump` instead.",
        action="store_true",
    )
    parser_type.add_argument(
        "--type",
        "-t",
        help="Input file type",
        choices=list(_Detect.subclasses.keys()),
        default="lammpsbondfile",
    )
    parser.add_argument("--nopbc", help="Disable PBC.", action="store_true")
    parser.add_argument(
        "--cell",
        "-c",
        nargs="+",
        type=float,
        help=(
            "Cell information if PBC is enabled and the input file does not contain cell information. "
            "The cell information should be in 3x3 matrix format."
        ),
    )
    parser.add_argument(
        "-n",
        "-np",
        "--nproc",
        help="Number of processes, by default using all processes",
        type=int,
    )
    parser.add_argument(
        "-s",
        "--selectatoms",
        help="Select element(s) in the reaction network, e.g. C",
        nargs="*",
    )
    parser.add_argument(
        "--stepinterval",
        help="Step interval, i.e. read every N step from the trajectory",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--split",
        help="Split number for the time axis; the whole trajectroy will be splited into N parts for analysis",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--maxspecies",
        help="Max number of nodes (species) in the network",
        type=int,
        default=20,
    )
    parser.add_argument(
        "--matrixa",
        help="Transition matrix A of HMM parameters",
        type=float,
        nargs=4,
        default=[0.999, 0.001, 0.001, 0.999],
    )
    parser.add_argument(
        "--matrixb",
        help="Emission matrix B of HMM parameters",
        type=float,
        nargs=4,
        default=[0.6, 0.4, 0.4, 0.6],
    )
    parser.add_argument(
        "--urls",
        action="append",
        nargs=2,
        type=str,
        help="Download files before analysis, in the format of `--urls filename url`",
    )
    parser.add_argument(
        "--version", action="version", version="ReacNetGenerator v%s" % __version__
    )
    return parser


def _commandline():
    parser = main_parser()
    args = parser.parse_args()
    import numpy as np

    from .reacnetgen import ReacNetGenerator

    ReacNetGenerator(
        inputfilename=args.inputfilename,
        atomname=args.atomname,
        miso=args.miso,
        runHMM=not args.nohmm,
        inputfiletype=("lammpsdumpfile" if args.dump else args.type),
        nproc=args.nproc,
        selectatoms=args.selectatoms,
        stepinterval=args.stepinterval,
        split=args.split,
        maxspecies=args.maxspecies,
        urls=[{"fn": url[0], "url": url[1]} for url in args.urls]
        if args.urls
        else None,
        a=np.array(args.matrixa).reshape((2, 2)),
        b=np.array(args.matrixb).reshape((2, 2)),
        pbc=not args.nopbc,
        cell=args.cell,
    ).runanddraw()


def parm2cmd(pp: dict) -> List[str]:
    """Convert a parameter dictionary to command line arguments.

    Parameters
    ----------
    pp : dict
        Parameter dictionary

    Returns
    -------
    List[str]
        Command line arguments
    """
    commands = ["reacnetgenerator", "-i", pp["inputfilename"], "-a", *pp["atomname"]]
    if not pp.get("runHMM", True):
        commands.append("--nohmm")
    if pp["inputfiletype"]:
        commands.extend(("--type", pp["inputfiletype"]))
    if pp["atomname"]:
        commands.extend(("-s", pp["atomname"][0]))
    if pp.get("urls", []):
        commands.extend(("--urls", pp["urls"][0]["fn"], pp["urls"][0]["url"][0]))
    if pp.get("a", []):
        commands.extend(
            (
                "--matrixa",
                str(pp["a"][0][0]),
                str(pp["a"][0][1]),
                str(pp["a"][1][0]),
                str(pp["a"][1][1]),
            )
        )
    if pp.get("b", []):
        commands.extend(
            (
                "--matrixb",
                str(pp["b"][0][0]),
                str(pp["b"][0][1]),
                str(pp["b"][1][0]),
                str(pp["b"][1][1]),
            )
        )
    if pp.get("cell", []):
        commands.extend(
            ("--cell", str(pp["cell"][0]), str(pp["cell"][1]), str(pp["cell"][2]))
        )
    for ii in ["nproc", "selectatoms", "stepinterval", "split", "maxspecies"]:
        if pp.get(ii, None):
            commands.extend((f"--{ii}", str(pp[ii])))
    return commands
