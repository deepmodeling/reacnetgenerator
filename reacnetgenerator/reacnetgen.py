#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
==================
ReacNetGenerator
==================
Automatic generator of reaction network for reactive molecular dynamics simulation.

Plase cite: J. Zeng, L. Cao, J.Z.H. Zhang, C.-H. Chin, T Zhu: ReacNetGen: an Automatic Reaction Network Generator for Reactive Molecular Dynamic Simulations, to be submitted

Author: Jinzhe Zeng, Liqun Cao, John ZH Zhang, Chih-Hao Chin, Tong Zhu
Email: tzhu@lps.ecnu.edu.cn, jzzeng@stu.ecnu.edu.cn

==================
Features
==================
* Processing of MD trajectory containing atomic coordinates or bond orders
* Hidden Markov Model (HMM) based noise filtering
* Isomers identifying accoarding to SMILES
* Generation of reaction network for visualization using force-directed algorithm
* Parallel computing

==================
Simple example
==================
Process a LAMMPS bond file named bonds.reaxc. (See http://lammps.sandia.gov/doc/fix_reax_bonds.html for details)
$ reacnetgenerator -i bonds.reaxc -a C H O
where C, H, and O are atomic names in the input file.

A LAMMPS dump file is also supported. You can prepare it by running "dump 1 all custom 100 dump.reaxc id type x y z" in LAMMPS. (See https://lammps.sandia.gov/doc/dump.html for more details)
$ reacnetgenerator --dump -i dump.reaxc -a C H O

You can running the following script for help:
$ reacnetgenerator -h
"""

__version__ = '1.2.20'
__date__ = '2018-03-11'
__update__ = '2019-01-13'
__author__ = 'Jinzhe Zeng'
__email__ = 'jzzeng@stu.ecnu.edu.cn'
__credits__ = ['Jinzhe Zeng', 'Tong Zhu',
               'Liqun Cao', 'Chih-Hao Chin', 'John ZH Zhang']
__copyright__ = 'Copyright 2018, East China Normal University'


import argparse
import base64
import gc
import itertools
import logging
import os
import tempfile
import time
import zlib
from enum import Enum
from multiprocessing import Semaphore, cpu_count

from pkg_resources import DistributionNotFound, get_distribution

import numpy as np

from ._detect import _Detect, InputFileType
from ._hmmfilter import _HMMFilter
from ._path import _CollectPaths
from ._matrix import _GenerateMatrix
from ._draw import _DrawNetwork
from ._reachtml import _HTMLResult


try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass

class ReacNetGenerator:
    ''' Use ReacNetGenerator for trajectory analysis'''

    def __init__(self, inputfiletype='lammpsbondfile', inputfilename='bonds.reaxc', atomname=None, selectatoms=None, originfilename=None, hmmfilename=None, atomfilename=None, moleculefilename=None, atomroutefilename=None, reactionfilename=None, tablefilename=None, moleculestructurefilename=None, imagefilename=None, speciesfilename=None, resultfilename=None, stepinterval=1, p=None, a=None, b=None, runHMM=True, SMILES=True, getoriginfile=False, species=None, node_size=200, font_size=6, widthcoefficient=1,  maxspecies=20, nolabel=False, needprintspecies=True, speciesfilter=None, node_color=None, pos=None, printfiltersignal=False, showid=True, k=None, start_color=None, end_color=None, nproc=None, speciescenter=None, n_searchspecies=2):
        ''' Init ReacNetGenerator '''
        print(__doc__)
        print(
            f"Version: {__version__}  Creation date: {__date__}  Update date: {__update__}")
        self.inputfiletype = InputFileType.LAMMPSBOND if inputfiletype=="lammpsbondfile" else InputFileType.LAMMPSDUMP
        self.inputfilename = inputfilename
        self.atomname = self._setparam(atomname, ["C", "H", "O"])
        self.selectatoms = self._setparam(selectatoms, self.atomname)
        self.originfilename = self._setfilename(originfilename, "origin")
        self.hmmfilename = self._setfilename(hmmfilename, "hmm")
        self.atomfilename = self._setfilename(atomfilename, "atom")
        self.moleculefilename = self._setfilename(moleculefilename, "moname")
        self.atomroutefilename = self._setfilename(atomroutefilename, "route")
        self.reactionfilename = self._setfilename(reactionfilename, "reaction")
        self.tablefilename = self._setfilename(tablefilename, "table")
        self.moleculestructurefilename = self._setfilename(
            moleculestructurefilename, "structure")
        self.imagefilename = self._setfilename(imagefilename, "svg")
        self.speciesfilename = self._setfilename(speciesfilename, "species")
        self.resultfilename = self._setfilename(resultfilename, "html")
        self.stepinterval = stepinterval
        self.p = self._setparam(np.array(p), np.array([0.5, 0.5]))
        self.a = self._setparam(np.array(a), np.array(
            [[0.999, 0.001], [0.001, 0.999]]))
        self.b = self._setparam(
            np.array(b), np.array([[0.6, 0.4], [0.4, 0.6]]))
        self.runHMM = runHMM
        self.SMILES = SMILES
        self.getoriginfile = getoriginfile if self.runHMM else True
        self.species = self._setparam(species, {})
        self.needprintspecies = needprintspecies
        self.node_size = node_size
        self.font_size = font_size
        self.widthcoefficient = widthcoefficient
        self.maxspecies = maxspecies
        self.nolabel = nolabel
        self.speciesfilter = self._setparam(speciesfilter, [])
        self.node_color = self._setparam(
            np.array(node_color), np.array([78/256, 196/256, 238/256]))
        self.pos = self._setparam(pos, {})
        self.printfiltersignal = printfiltersignal
        self.showid = showid
        self.k = k
        self.start_color = self._setparam(
            np.array(start_color), np.array([0, 0, 1]))
        self.end_color = self._setparam(
            np.array(end_color), np.array([1, 0, 0]))
        self.nproc = self._setparam(nproc, cpu_count())
        self.speciescenter = speciescenter
        self.n_searchspecies = n_searchspecies
        # define attribute
        self._atomtype = None
        self._step = None
        self._hmmit = None
        self._timestep = None
        self._steplinenum = None
        self._N = None
        self._temp1it = None
        self.moleculetempfilename = None
        self.moleculetemp2filename = None
        self.allmoleculeroute = None

    def runanddraw(self, run=True, draw=True, report=True):
        ''' Analyze the trajectory from MD simulation '''
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
        """ Processing of MD trajectory """
        self._process((
            self.Status.DETECT,
            self.Status.HMM,
            self.Status.PATH,
            self.Status.MATRIX,
        ))

    def draw(self):
        """ Draw the reaction network """
        self._process((self.Status.NETWORK))

    def report(self):
        """ Generate the analysis report """
        self._process((self.Status.REPORT))

    class Status(Enum):
        INIT = "Init"
        DETECT = "Read bond information and Detect molecules"
        HMM = "HMM filter"
        PATH = "Indentify isomers and collect reaction paths"
        MATRIX = "Reaction matrix generation"
        NETWORK = "Draw reaction network"
        REPORT = "Generate analysis report"

        def __str__(self):
            return self.value

    def _process(self, steps):
        timearray = [time.time()]
        for i, runstep in enumerate(steps,1):
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
        for tempfilename in (self.moleculetempfilename, self.moleculetemp2filename):
            if tempfilename is not None:
                try:
                    os.remove(tempfilename)
                except OSError:
                    pass
        # Summary
        logging.info("====== Summary ======")
        logging.info(f"Total time(s): {timearray[-1]-timearray[0]:.3f} s")

    @classmethod
    def _produce(cls, semaphore, plist, parameter):
        for item in plist:
            semaphore.acquire()
            yield item, parameter

    @classmethod
    def _compress(cls, x):
        return base64.a85encode(zlib.compress(x.encode()))+b'\n'

    @classmethod
    def _decompress(cls, x):
        return zlib.decompress(base64.a85decode(x.strip())).decode()

    @classmethod
    def _setparam(cls, x, default):
        return x if x else default

    def _setfilename(self, name, suffix):
        return self._setparam(name, f"{self.inputfilename}.{suffix}")

def _commandline():
    parser = argparse.ArgumentParser(
        description=f'ReacNetGenerator {__version__}')
    parser.add_argument('-i', '--inputfilename',
                        help='Input trajectory file, e.g. bonds.reaxc', required=True)
    parser.add_argument(
        '-a', '--atomname', help='Atomic names in the trajectory, e.g. C H O', nargs='*', required=True)
    parser.add_argument(
        '--nohmm', help='Process trajectory without Hidden Markov Model (HMM)', action="store_true")
    parser.add_argument(
        '--dump', help='Process the LAMMPS dump file', action="store_true")
    parser.add_argument(
        '-n', '-np', '--nproc', help='Number of processes', type=int)
    parser.add_argument('-s', '--selectatoms',
                        help='Select atoms in the reaction network, e.g. C', nargs='*')
    args = parser.parse_args()
    ReacNetGenerator(inputfilename=args.inputfilename, atomname=args.atomname, runHMM=not args.nohmm,
                     inputfiletype=('lammpsdumpfile' if args.dump else 'lammpsbondfile'), nproc=args.nproc, selectatoms=args.selectatoms).runanddraw()
