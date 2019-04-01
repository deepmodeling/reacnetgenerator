#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
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


import argparse
import gc
import logging
import os
import time
import pickle
from enum import Enum
from multiprocessing import cpu_count

import numpy as np
import lz4.frame
import pybase64

from . import __version__, __date__, __update__
from ._detect import InputFileType, _Detect
from ._draw import _DrawNetwork
from ._hmmfilter import _HMMFilter
from ._matrix import _GenerateMatrix
from ._path import _CollectPaths
from ._reachtml import _HTMLResult


class ReacNetGenerator:
    """Use ReacNetGenerator for trajectory analysis."""

    def __init__(
            self, inputfiletype='lammpsbondfile', inputfilename='bonds.reaxc',
            atomname=None, selectatoms=None, originfilename=None,
            hmmfilename=None, atomfilename=None, moleculefilename=None,
            atomroutefilename=None, reactionfilename=None, tablefilename=None,
            imagefilename=None, speciesfilename=None, resultfilename=None,
            stepinterval=1, p=None, a=None, b=None, runHMM=True, SMILES=True,
            getoriginfile=False, species=None, node_size=200, font_size=6,
            widthcoefficient=1, maxspecies=20, nolabel=False,
            needprintspecies=True, speciesfilter=None, node_color=None,
            pos=None, printfiltersignal=False, showid=True, k=None,
            start_color=None, end_color=None, nproc=None, speciescenter=None,
            n_searchspecies=2):
        """Init ReacNetGenerator."""
        print(__doc__)
        print(
            f"Version: {__version__}  Creation date: {__date__}  Update date: {__update__}")
        self.inputfiletype = InputFileType.LAMMPSBOND if inputfiletype == "lammpsbondfile" else InputFileType.LAMMPSDUMP
        self.inputfilename = inputfilename
        self.atomname = np.array(self._setparam(atomname, ["C", "H", "O"]))
        self.selectatoms = self._setparam(selectatoms, self.atomname)
        self.moleculefilename = self._setfilename(moleculefilename, "moname")
        self.atomroutefilename = self._setfilename(atomroutefilename, "route")
        self.reactionfilename = self._setfilename(reactionfilename, "reaction")
        self.tablefilename = self._setfilename(tablefilename, "table")
        self.imagefilename = self._setfilename(imagefilename, "svg")
        self.speciesfilename = self._setfilename(speciesfilename, "species")
        self.resultfilename = self._setfilename(resultfilename, "html")
        self.stepinterval = stepinterval
        self.p = np.array(self._setparam(p, [0.5, 0.5]))
        self.a = np.array(self._setparam(a, [[0.999, 0.001], [0.001, 0.999]]))
        self.b = np.array(self._setparam(b, [[0.6, 0.4], [0.4, 0.6]]))
        self.runHMM = runHMM
        self.SMILES = SMILES
        self.getoriginfile = getoriginfile if self.runHMM else True
        self.species = np.array(self._setparam(species, []))
        self.needprintspecies = needprintspecies
        self.node_size = node_size
        self.font_size = font_size
        self.widthcoefficient = widthcoefficient
        self.maxspecies = maxspecies
        self.nolabel = nolabel
        self.speciesfilter = self._setparam(speciesfilter, [])
        self.node_color = np.array(self._setparam(
            node_color, [78/256, 196/256, 238/256]))
        self.pos = self._setparam(pos, {})
        self.printfiltersignal = printfiltersignal
        self.showid = showid
        self.k = k
        self.start_color = np.array(self._setparam(start_color, [0, 0, 1]))
        self.end_color = np.array(self._setparam(end_color, [1, 0, 0]))
        self.nproc = self._setparam(nproc, cpu_count())
        self.speciescenter = speciescenter
        self.n_searchspecies = n_searchspecies
        # define attribute
        self.atomtype = None
        self.step = None
        self.hmmit = None
        self.timestep = None
        self.steplinenum = None
        self.N = None
        self.temp1it = None
        self.originfilename = None
        self.hmmfilename = None
        self.moleculetempfilename = None
        self.moleculetemp2filename = None
        self.allmoleculeroute = None

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

    @classmethod
    def produce(cls, semaphore, plist, parameter):
        """Prevent large memory usage due to slow IO."""
        for item in plist:
            semaphore.acquire()
            yield item, parameter

    @classmethod
    def compress(cls, x, isbytes=False):
        """Compress the line.

        This function reduces IO overhead to speed up the program.
        """
        if isbytes:
            return pybase64.b64encode(lz4.frame.compress(x, compression_level=0))+b'\n'
        return pybase64.b64encode(lz4.frame.compress(x.encode(), compression_level=-1))+b'\n'

    @classmethod
    def decompress(cls, x, isbytes=False):
        """Decompress the line."""
        if isbytes:
            return lz4.frame.decompress(pybase64.b64decode(x.strip(), validate=True))
        return lz4.frame.decompress(pybase64.b64decode(x.strip(), validate=True)).decode()

    @classmethod
    def _setparam(cls, x, default):
        return x if x is not None else default

    def _setfilename(self, name, suffix):
        return self._setparam(name, f"{self.inputfilename}.{suffix}")

    @classmethod
    def listtobytes(cls, x):
        return cls.compress(pickle.dumps(x), isbytes=True)

    @classmethod
    def bytestolist(cls, x):
        return pickle.loads(cls.decompress(x, isbytes=True))


def _commandline():
    parser = argparse.ArgumentParser(
        description=f'ReacNetGenerator {__version__}')
    parser.add_argument(
        '-i', '--inputfilename',
        help='Input trajectory file, e.g. bonds.reaxc', required=True)
    parser.add_argument('-a', '--atomname',
                        help='Atomic names in the trajectory, e.g. C H O',
                        nargs='*', required=True)
    parser.add_argument(
        '--nohmm', help='Process trajectory without Hidden Markov Model (HMM)',
        action="store_true")
    parser.add_argument(
        '--dump', help='Process the LAMMPS dump file', action="store_true")
    parser.add_argument(
        '-n', '-np', '--nproc', help='Number of processes', type=int)
    parser.add_argument(
        '-s', '--selectatoms',
        help='Select atoms in the reaction network, e.g. C', nargs='*')
    parser.add_argument(
        '--stepinterval', help='Step interval', type=int, default=1)
    args = parser.parse_args()
    ReacNetGenerator(
        inputfilename=args.inputfilename, atomname=args.atomname,
        runHMM=not args.nohmm,
        inputfiletype=('lammpsdumpfile' if args.dump else 'lammpsbondfile'),
        nproc=args.nproc, selectatoms=args.selectatoms,
        stepinterval=args.stepinterval
        ).runanddraw()
