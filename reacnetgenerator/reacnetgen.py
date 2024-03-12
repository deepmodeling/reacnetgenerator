#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# SPDX-License-Identifier: LGPL-3.0-or-later
# cython: language_level=3
# cython: linetrace=True
"""The main module of ReacNetGenerator."""
doc_run = """ReacNetGenerator: an automatic reaction network generator for reactive
molecular dynamics simulation.

Please cite: ReacNetGenerator: an automatic reaction network generator
for reactive molecular dynamic simulations, Phys. Chem. Chem. Phys.,
2020, 22 (2): 683-691, doi: 10.1039/C9CP05091D

Jinzhe Zeng (jinzhe.zeng@rutgers.edu), Tong Zhu (tzhu@lps.ecnu.edu.cn)

==================
Features
==================
* Processing of MD trajectory containing atomic coordinates or bond orders
* Hidden Markov Model (HMM) based noise filtering
* Isomers identifying accoarding to SMILES
* Generation of reaction network for visualization using force-directed
  algorithm
* Parallel computing

==================
Simple example
==================
ReacNetGenerator can process any kind of trajectory files containing
atomic coordinates, e.g. a LAMMPS dump file prepared by running “dump 1
all custom 100 dump.reaxc id type x y z” in LAMMPS:
$ reacnetgenerator --type dump -i dump.reaxc -a C H O --nohmm
where C, H, and O are atomic names in the input file. Analysis report
will be generated automatically.

Also, ReacNetGenerator can process files containing bond information,
e.g. LAMMPS bond file:
$ reacnetgenerator --type bond -i bonds.reaxc -a C H O --nohmm

You can running the following script for help:
$ reacnetgenerator -h
"""


import gc
import itertools
import os
import time
from enum import Enum
from typing import Any, List, Tuple, Union

import numpy as np

from . import __date__, __version__
from ._detect import _Detect
from ._download import DownloadData
from ._draw import _DrawNetwork
from ._hmmfilter import _HMMFilter
from ._logging import logger
from ._matrix import _GenerateMatrix
from ._mergeiso import _mergeISO
from ._path import _CollectPaths
from ._reachtml import _HTMLResult
from .utils import must_be_list


class ReacNetGenerator:
    """Use ReacNetGenerator for trajectory analysis.

    Parameters
    ----------
    inputfiletype: str
        The type of the input file. The following type is allowed:
            - dump: LAMMPS dump file, which can be outputed by
              using `dump 1 all custom 100 dump.reaxc id type x y z`. See
              https://lammps.sandia.gov/doc/dump.html for details.
            - bond: LAMMPS ReaxFF bond file. See https://lammps.sandia.gov/doc/fix_reaxc_bonds.html
              for details.
            - xyz: XYZ file, which can also be outputed by LAMMPS using dump.
    inputfilename: str or list of strs
        The filename(s) of the input file, which can be either relative path or absolute path. If
        it is a list, the files will be read in order.
    atomname: tuple of strs
        The list of the atomic names in the input file, such as `('C', 'H', 'O')`. It should match
        the order of that in the input file.
    runHMM: bool, optional, default: True
        Process trajectory with Hidden Markov Model (HMM) or not. If the user find too many species
        are filtered, they can turn off this option.
    miso: int, optional, default: 0
        Merge the isomers and the highest frequency is used as the representative. 0, off
        two available levels:
        1, merge the isomers with same atoms and same bond-network but different bond levels;
        2, merge the isomers with same atoms with different bond-network.
    pbc: bool, optional, default: True
        Use periodic boundary conditions (PBC) or not.
    cell: (3,3) array_like or (3,) array_like or (9,) array_like, optional, default: None
        The cell (box size) of the system. If None (default), the cell will be read from the input
        file. If the input file doesn't have cell information, this parameter will be necessary.
    nproc: int, optional, default: None
        The number of processors used for analysis. If None (default), the program will try to use
        all processors.
    selectatoms: str, optional, default: None
        Select an element from the atomic names, such as `C`, and only show species with this
        element in the reaction network. If None (default), the network will show all elements.
    split: int, optional, default: None
        Split number for the time axis. For example, if set to 10, the whole trajectroy will
        be divided into 10 parts and reactions of each part will be shown.
    a: (2,2) array_like, optional, default: [[0.999, 0.001], [0.001, 0.009]]
        Transition matrix A of HMM parameters. It is recommended for users to choose their own
        parameters. See the paper for details.
    b: (2,2) array_like, optional, default: [[0.6, 0.4], [0.4, 0.6]]
        Emission matrix B of HMM parameters. It is recommended for users to choose their own
        parameters. See the paper for details.

    Examples
    --------
    >>> from reacnetgenerator import ReacNetGenerator
    >>> rng=ReacNetGenerator(inputfiletype="dump", inputfilename="dump.ch4", atomname=['C', 'H', 'O'])
    >>> rng.runanddraw()
    """

    urls: dict
    moleculetempfilename: str
    moleculetemp2filename: str
    originfilename: str
    hmmfilename: str
    resultfilename: str

    def __init__(self, **kwargs: Any) -> None:
        """Init ReacNetGenerator."""
        logger.info(doc_run)
        logger.info(f"Version: {__version__}  Creation date: {__date__}")

        try:
            nproc = len(os.sched_getaffinity(0))
        except AttributeError:
            # macos and windows
            nproc = os.cpu_count()

        # process kwargs
        necessary_key = ["inputfiletype", "inputfilename", "atomname"]
        default_value = {
            "node_color": [78 / 256, 196 / 256, 238 / 256],
            "p": [0.5, 0.5],
            "a": [[0.999, 0.001], [0.001, 0.999]],
            "b": [[0.6, 0.4], [0.4, 0.6]],
            "speciesfilter": [],
            "start_color": [0, 0, 1],
            "end_color": [1, 0, 0],
            "nproc": nproc,
            "pos": {},
            "pbc": True,
            "split": 1,
            "n_searchspecies": 2,
            "node_size": 200,
            "font_size": 6,
            "widthcoefficient": 1,
            "maxspecies": 20,
            "stepinterval": 1,
            "nolabel": False,
            "printfiltersignal": False,
            "showid": True,
            "runHMM": True,
            "SMILES": True,
            "miso": 0,
            "getoriginfile": False,
            "needprintspecies": True,
            "urls": [],
            "matrix_size": 100,
        }
        none_key = ["selectatoms", "species", "pos", "k", "speciescenter", "cell"]
        accept_keys = [
            "atomtype",
            "step",
            "hmmit",
            "timestep",
            "steplinenum",
            "N",
            "temp1it",
            "originfilename",
            "hmmfilename",
            "moleculetempfilename",
            "moleculetemp2filename",
            "allmoleculeroute",
            "splitmoleculeroute",
        ]
        nparray_key = [
            "atomname",
            "p",
            "a",
            "b",
            "node_color",
            "start_color",
            "end_color",
        ]
        file_key = {
            "moleculefilename": "moname",
            "atomroutefilename": "route",
            "reactionfilename": "reaction",
            "tablefilename": "table",
            "imagefilename": "svg",
            "speciesfilename": "species",
            "resultfilename": "html",
            "jsonfilename": "json",
            "reactionabcdfilename": "reactionabcd",
        }
        assert set(necessary_key).issubset(set(kwargs)), (
            "Must give neccessary key: %s" % ", ".join(necessary_key)
        )
        assert set(kwargs).issubset(
            set(necessary_key) | set(default_value) | set(none_key) | set(file_key)
        ), "Unsupported key"
        kwargs["inputfilename"] = must_be_list(kwargs["inputfilename"])
        for kk in default_value:
            kwargs.setdefault(kk, default_value[kk])
            if kwargs[kk] is None:
                kwargs[kk] = default_value[kk]
        for kk in itertools.chain(none_key, accept_keys):
            kwargs.setdefault(kk, None)
        for kk in file_key:
            kwargs.setdefault(kk, f"{kwargs['inputfilename'][0]}.{file_key[kk]}")
        for kk in nparray_key:
            kwargs[kk] = np.array(kwargs[kk])
        if not kwargs["runHMM"]:
            kwargs["getoriginfile"] = True
        if kwargs["selectatoms"] is None:
            kwargs["selectatoms"] = kwargs["atomname"]
        self.__dict__.update(kwargs)
        if self.cell is not None:
            if len(self.cell) == 9:
                self.cell = np.array(self.cell).reshape((3, 3))
            elif len(self.cell) == 3:
                self.cell = np.diag(self.cell)
            else:
                raise RuntimeError(
                    "cell must be (3,3) array_like or (3,) array_like or (9,) array_like"
                )

    def runanddraw(
        self, run: bool = True, draw: bool = True, report: bool = True
    ) -> None:
        """Analyze the trajectory from MD simulation.

        Parameters
        ----------
        run : bool, optional, default: True
            Process the trajectory or not, including DOWNLOAD, DETECT, HMM, PATH, and MATRIX steps.
        draw : bool, optional, default: True
            Draw the reaction network or not, i.e. NETWORK step.
        report : bool, optional, default: True
            Generate the analysis report, i.e. NETWORK step.
        """
        processthing = []
        if run:
            if self.urls:
                processthing.append(self.Status.DOWNLOAD)
            processthing.extend(
                (
                    self.Status.DETECT,
                    self.Status.MISO,
                    self.Status.HMM,
                    self.Status.PATH,
                    self.Status.MATRIX,
                )
            )
        if draw:
            processthing.append(self.Status.NETWORK)
        if report:
            processthing.append(self.Status.REPORT)
        self._process(processthing)

    def run(self) -> None:
        """Process MD trajectory, including DOWNLOAD, DETECT, HMM, PATH, and MATRIX steps."""
        processthing = []
        if self.urls:
            processthing.append(self.Status.DOWNLOAD)
        processthing.extend(
            (
                self.Status.DETECT,
                self.Status.MISO,
                self.Status.HMM,
                self.Status.PATH,
                self.Status.MATRIX,
            )
        )
        self._process(processthing)

    def draw(self) -> None:
        """Draw the reaction network, i.e. NETWORK step."""
        self._process((self.Status.NETWORK,))

    def report(self) -> None:
        """Generate the analysis report, i.e. REPORT step."""
        self._process((self.Status.REPORT,))

    class Status(Enum):
        """ReacNetGenerator status.

        The ReacNetGenerator consists of several modules and algorithms to
        process the information from the given trajectory, including:

        - DOWNLOAD: Download trajectory from urls
        - DETECT: Read bond information and detect molecules
        - HMM: HMM filter
        - MISO: Merge isomers
        - PATH: Indentify isomers and collect reaction paths
        - MATRIX: Reaction matrix generation
        - NETWORK: Draw reaction network
        - REPORT: Generate analysis report
        """

        INIT = "Init"
        DETECT = "Read bond information and detect molecules"
        MISO = "Merge isomers"
        HMM = "HMM filter"
        PATH = "Indentify isomers and collect reaction paths"
        MATRIX = "Reaction matrix generation"
        NETWORK = "Draw reaction network"
        REPORT = "Generate analysis report"
        DOWNLOAD = "Download trajectory"

        def __str__(self):
            """Return describtion of the status."""
            return self.value

    def _process(self, steps: Union[List[Status], Tuple[Status, ...]]) -> None:
        """Process steps in order.

        Parameters
        ----------
        steps : tuple of ReacNetGenerator.Status
            The process that needs to be processed.
        """
        timearray = [time.perf_counter()]
        for i, runstep in enumerate(steps, 1):
            if runstep == self.Status.DETECT:
                _Detect.gettype(self).detect()
            elif runstep == self.Status.MISO:
                _mergeISO(self).merge()
            elif runstep == self.Status.HMM:
                _HMMFilter(self).filter()
            elif runstep == self.Status.PATH:
                _CollectPaths.getstype(self).collect()
            elif runstep == self.Status.MATRIX:
                _GenerateMatrix(self).generate()
            elif runstep == self.Status.NETWORK:
                _DrawNetwork(self).draw()
            elif runstep == self.Status.REPORT:
                _HTMLResult(self).report()
            elif runstep == self.Status.DOWNLOAD:
                DownloadData(self).download_files()
            # garbage collect
            gc.collect()
            timearray.append(time.perf_counter())
            logger.info(
                f"Step {i}: Done! Time consumed (s): {timearray[-1]-timearray[-2]:.3f} ({runstep})"
            )

        # delete tempfile
        for tempfilename in (
            self.moleculetempfilename,
            self.moleculetemp2filename,
            self.originfilename,
            self.hmmfilename,
        ):
            if tempfilename is not None:
                try:
                    os.remove(tempfilename)
                except OSError:
                    pass
        # Summary
        logger.info("====== Summary ======")
        logger.info(f"Total time(s): {timearray[-1]-timearray[0]:.3f} s")


# Please import ReacNetGenerator class from reacnetgenerator instead of reacnetgenerator.reacnetgen
__all__ = []
