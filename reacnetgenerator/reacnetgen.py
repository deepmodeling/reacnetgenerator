#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# cython: language_level=3
# cython: linetrace=True
"""ReacNetGenerator.

Automatic generator of reaction network for reactive molecular dynamics
simulation.

Plase cite: J. Zeng, L. Cao, J.Z.H. Zhang, C.-H. Chin, T Zhu: ReacNetGen: an
Automatic Reaction Network Generator for Reactive Molecular Dynamic
Simulations, doi: 10.26434/chemrxiv.7421534

Author: Jinzhe Zeng, Liqun Cao, John ZH Zhang, Chih-Hao Chin, Tong Zhu

Email: tzhu@lps.ecnu.edu.cn, jzzeng@stu.ecnu.edu.cn

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
Process a LAMMPS bond file named bonds.reaxc. (See
http://lammps.sandia.gov/doc/fix_reax_bonds.html for details)
$ reacnetgenerator -i bonds.reaxc -a C H O
where C, H, and O are atomic names in the input file.

A LAMMPS dump file is also supported. You can prepare it by running "dump 1
all custom 100 dump.reaxc id type x y z" in LAMMPS. (See
https://lammps.sandia.gov/doc/dump.html for more details)
$ reacnetgenerator --dump -i dump.reaxc -a C H O

You can running the following script for help:
$ reacnetgenerator -h
"""


import gc
import logging
import os
import time
import itertools
from enum import Enum
from multiprocessing import cpu_count

import numpy as np

from . import __version__, __date__, __update__
from ._detect import InputFileType, _Detect
from ._draw import _DrawNetwork
from ._hmmfilter import _HMMFilter
from ._matrix import _GenerateMatrix
from ._path import _CollectPaths
from ._reachtml import _HTMLResult


class ReacNetGenerator:
    """Use ReacNetGenerator for trajectory analysis."""

    def __init__(self, **kwargs):
        """Init ReacNetGenerator."""
        logging.info(__doc__)
        logging.info(
            f"Version: {__version__}  Creation date: {__date__}  Update date: {__update__}")
        
        # process kwargs
        necessary_key = ['inputfiletype', 'inputfilename', 'atomname']
        default_value = {"node_color": [78/256, 196/256, 238/256], "p": [0.5, 0.5],
                            "a": [[0.999, 0.001], [0.001, 0.999]], "b": [[0.6, 0.4], [0.4, 0.6]],
                            "speciesfilter": [], "start_color": [0, 0, 1], "end_color": [1, 0, 0],
                            "nproc": cpu_count(), "pos": {}, "pbc": True, "split": 1, "n_searchspecies": 2,
                            "node_size": 200, "font_size": 6, "widthcoefficient": 1, "maxspecies": 20, "stepinterval": 1,
                            "nolabel": False,  "printfiltersignal": False, "showid": True, "runHMM": True, "SMILES": True,
                            "getoriginfile": False, "needprintspecies": True
                            }
        none_key = ['selectatoms', 'species', 'pos', 'k', 'speciescenter']
        accept_keys = ['atomtype', 'step', 'hmmit', 'timestep', 'steplinenum', 'N',
            'temp1it', 'originfilename', 'hmmfilename', 'moleculetempfilename', 'moleculetemp2filename',
            'allmoleculeroute', 'splitmoleculeroute']
        nparray_key = ["atomname", "p", "a", "b",
                        "node_color", "start_color", "end_color"]
        file_key = {"moleculefilename": "moname", "atomroutefilename": "route", "reactionfilename": "reaction",
                    "tablefilename": "table", "imagefilename": "svg", "speciesfilename": "species", "resultfilename": "html",
                    "jsonfilename": "json", "reactionabcdfilename": "reactionabcd"}
        assert set(necessary_key).issubset(
            set(kwargs)), "Must give neccessary key: %s" % ", ".join(necessary_key)
        assert set(kwargs).issubset(
            set(necessary_key) | set(default_value) | set(none_key) | set(file_key)), "Unsupported key"
        if kwargs["inputfiletype"] == "lammpsbondfile":
            kwargs["inputfiletype"] = InputFileType.LAMMPSBOND
        elif kwargs["inputfiletype"] == "lammpsdumpfile":
            kwargs["inputfiletype"] = InputFileType.LAMMPSDUMP
        else:
            raise RuntimeError("Unsupported file format!")
        if isinstance(kwargs["inputfilename"], str):
            kwargs["inputfilename"] = [kwargs["inputfilename"]]
        for kk in default_value:
            kwargs.setdefault(kk, default_value[kk])
        for kk in itertools.chain(none_key, accept_keys):
            kwargs.setdefault(kk, None)
        for kk in file_key:
            kwargs.setdefault(
                kk, f"{kwargs['inputfilename'][0]}.{file_key[kk]}")
        for kk in nparray_key:
            kwargs[kk] = np.array(kwargs[kk])
        if not kwargs["runHMM"]:
            kwargs["getoriginfile"] = True
        if kwargs["selectatoms"] is None:
            kwargs["selectatoms"] = kwargs["atomname"]
        self.__dict__.update(kwargs)

    def runanddraw(self, run=True, draw=True, report=True):
        """Analyze the trajectory from MD simulation."""
        processthing = []
        if run:
            processthing.extend((
                self.Status.DETECT,
                self.Status.HMM,
                self.Status.PATH,
                self.Status.MATRIX,
            ))
        if draw:
            processthing.append(self.Status.NETWORK)
        if report:
            processthing.append(self.Status.REPORT)
        self._process(processthing)

    def run(self):
        """Process MD trajectory."""
        self._process((
            self.Status.DETECT,
            self.Status.HMM,
            self.Status.PATH,
            self.Status.MATRIX,
        ))

    def draw(self):
        """Draw the reaction network."""
        self._process((self.Status.NETWORK,))

    def report(self):
        """Generate the analysis report."""
        self._process((self.Status.REPORT,))

    class Status(Enum):
        """ReacNetGen status.

        The ReacNetGen consists of several modules and algorithms to
        process the information from the given trajectory.
        """

        INIT = "Init"
        DETECT = "Read bond information and Detect molecules"
        HMM = "HMM filter"
        PATH = "Indentify isomers and collect reaction paths"
        MATRIX = "Reaction matrix generation"
        NETWORK = "Draw reaction network"
        REPORT = "Generate analysis report"

        def __str__(self):
            """Return describtion of the status."""
            return self.value

    def _process(self, steps):
        timearray = [time.time()]
        for i, runstep in enumerate(steps, 1):
            if runstep == self.Status.DETECT:
                _Detect.gettype(self.inputfiletype)(self).detect()
            elif runstep == self.Status.HMM:
                _HMMFilter(self).filter()
            elif runstep == self.Status.PATH:
                _CollectPaths.getstype(self.SMILES)(self).collect()
            elif runstep == self.Status.MATRIX:
                _GenerateMatrix(self).generate()
            elif runstep == self.Status.NETWORK:
                _DrawNetwork(self).draw()
            elif runstep == self.Status.REPORT:
                _HTMLResult(self).report()
            # garbage collect
            gc.collect()
            timearray.append(time.time())
            logging.info(
                f"Step {i}: Done! Time consumed (s): {timearray[-1]-timearray[-2]:.3f} ({runstep})")

        # delete tempfile
        for tempfilename in (
                self.moleculetempfilename, self.moleculetemp2filename, self.originfilename, self.hmmfilename):
            if tempfilename is not None:
                try:
                    os.remove(tempfilename)
                except OSError:
                    pass
        # Summary
        logging.info("====== Summary ======")
        logging.info(f"Total time(s): {timearray[-1]-timearray[0]:.3f} s")
