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
import math
import os
import tempfile
import time
import zlib
from enum import Enum
from collections import Counter
from multiprocessing import Pool, Semaphore, cpu_count

from pkg_resources import DistributionNotFound, get_distribution

from tqdm import tqdm
from rdkit import Chem
import numpy as np
import networkx.algorithms.isomorphism as iso
import networkx as nx

from ._detect import _Detect, InputFileType
from ._hmmfilter import _HMMFilter
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
        self._timearray = []
        self._status = self.Status.INIT
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
        self._printtime()
        for runstep in steps:
            self._status = runstep
            if runstep == self.Status.DETECT:
                _Detect.gettype(self.inputfiletype)(self).detect()
            elif runstep == self.Status.HMM:
                _HMMFilter(self).filter()
            elif runstep == self.Status.PATH:
                if self.SMILES:
                    self._printmoleculeSMILESname()
                else:
                    self._printmoleculename()
                atomeach = self._getatomeach()
                allmoleculeroute = self._printatomroute(atomeach)
            elif runstep == self.Status.MATRIX:
                allroute = self._getallroute(allmoleculeroute)
                self._printtable(allroute)
                if self.needprintspecies:
                    self._printspecies()
            elif runstep == self.Status.NETWORK:
                _DrawNetwork(self).draw()
            elif runstep == self.Status.REPORT:
                _HTMLResult(self).report()
                logging.info(
                    f"Report is generated. Please see {self.resultfilename} for more details.")
            # garbage collect
            gc.collect()
            self._printtime()

        # delete tempfile
        for tempfilename in (self.moleculetempfilename, self.moleculetemp2filename):
            if tempfilename is not None:
                try:
                    os.remove(tempfilename)
                except OSError:
                    pass
        # Summary
        self._summary()

    def _printtime(self):
        self._timearray.append(time.time())
        if self._status != self.Status.INIT:
            logging.info(
                f"Step {len(self._timearray)-1}: Done! Time consumed (s): {self._timearray[-1]-self._timearray[-2]:.3f} ({self._status})")

    def _summary(self):
        logging.info("====== Summary ======")
        logging.info(
            f"Total time(s): {self._timearray[-1]-self._timearray[0]:.3f} s")

    @classmethod
    def _produce(cls, semaphore, plist, parameter):
        for item in plist:
            semaphore.acquire()
            yield item, parameter

    def _getatomroute(self, item):
        (i, (atomeachi, atomtypei)), _ = item
        routestrarr = []
        moleculeroute = []
        right = -1
        for j, atomeachij in enumerate(atomeachi.tolist()):
            if atomeachij > 0 and atomeachij != right:
                routestrarr.append(
                    f"{self._mname[atomeachij-1]} ({atomeachij} step { self._timestep[j]})")
                left, right = right, atomeachij
                if self.atomname[atomtypei-1] in self.selectatoms:
                    if left >= 0 and (left, right) not in moleculeroute:
                        moleculeroute.append((left, right))
        routestr = f"Atom {i} {self.atomname[atomtypei-1]}: " + \
            " -> ".join(routestrarr)
        return moleculeroute, routestr

    def _printatomroute(self, atomeach):
        with open(self.atomroutefilename, 'w') as f, Pool(self.nproc, maxtasksperchild=1000) as pool:
            allmoleculeroute = []
            semaphore = Semaphore(self.nproc*15)
            results = pool.imap(self._getatomroute, self._produce(
                semaphore, enumerate(zip(atomeach, self._atomtype), start=1), ()), 10)
            for route in tqdm(results, total=self._N, desc="Collect reaction paths", unit="atom"):
                moleculeroute, routestr = route
                print(routestr, file=f)
                for mroute in moleculeroute:
                    if mroute not in allmoleculeroute:
                        allmoleculeroute.append(mroute)
                semaphore.release()
        return allmoleculeroute

    @classmethod
    def _makemoleculegraph(cls, atoms, bonds):
        G = nx.Graph()
        for line in bonds:
            G.add_edge(line[0], line[1], level=line[2])
        for atom in atoms:
            atomnumber, atomtype = atom
            G.add_node(atomnumber, atom=atomtype)
        return G

    def _getstructure(self, name, atoms, bonds):
        index = {}
        for i, atom in enumerate(atoms, start=1):
            index[atom] = i
        return " ".join((name, ",".join([self.atomname[self._atomtype[x-1]-1] for x in atoms]), ";".join([",".join((str(index[x[0]]), str(index[x[1]]), str(x[2]))) for x in bonds])))

    def _readstrcture(self):
        with open(self.moleculestructurefilename) as f:
            d = {}
            for line in f:
                s = line.split()
                name = s[0]
                atoms = [x for x in s[1].split(",")]
                bonds = [tuple(int(y) for y in x.split(","))
                         for x in s[2].split(";")] if len(s) == 3 else []
                d[name] = (atoms, bonds)
        return d

    def _printmoleculename(self):
        mname = []
        d = {}
        em = iso.numerical_edge_match(['atom', 'level'], ["None", 1])
        with open(self.moleculefilename, 'w') as fm, open(self.moleculetemp2filename, 'rb') as ft, open(self.moleculestructurefilename, 'w') as fs:
            for line in ft:
                s = self._decompress(line).split()
                atoms = np.array([int(x) for x in s[0].split(",")])
                bonds = np.array([tuple(int(y) for y in x.split(","))
                                  for x in s[1].split(";")] if len(s) == 3 else [])
                typenumber = np.zeros(len(self.atomname), dtype=np.int)
                atomtypes = []
                for atomnumber in atoms:
                    typenumber[self._atomtype[atomnumber-1]-1] += 1
                    atomtypes.append(
                        (atomnumber, self._atomtype[atomnumber-1]))
                G = self._makemoleculegraph(atomtypes, bonds)
                name = "".join([self.atomname[i]+(str(typenumber[i] if typenumber[i] > 1 else ""))
                                if typenumber[i] > 0 else "" for i in range(0, len(self.atomname))])
                if name in d:
                    for j in range(len(d[name])):
                        if nx.is_isomorphic(G, d[name][j], em):
                            if j > 0:
                                name += "_"+str(j+1)
                            break
                    else:
                        d[name].append(G)
                        name += "_"+str(len(d[name]))
                        print(self._getstructure(name, atoms, bonds), file=fs)
                else:
                    d[name] = [G]
                    print(self._getstructure(name, atoms, bonds), file=fs)
                mname.append(name)
                print(name, ",".join([str(x) for x in atoms]), ";".join(
                    [",".join([str(y) for y in x]) for x in bonds]), file=fm)
        self._mname = mname

    def _calmoleculeSMILESname(self, item):
        line, _ = item
        s = self._decompress(line).split()
        atoms = [int(x) for x in s[0].split(",")]
        bonds = [tuple(int(y) for y in x.split(","))
                 for x in s[1].split(";")] if len(s) == 3 else []
        name = self._convertSMILES(atoms, bonds)
        return name, atoms, bonds

    def _printmoleculeSMILESname(self):
        mname = []
        with open(self.moleculefilename, 'w') as fm, open(self.moleculetemp2filename, 'rb') as ft, Pool(self.nproc, maxtasksperchild=1000) as pool:
            semaphore = Semaphore(self.nproc*15)
            results = pool.imap(self._calmoleculeSMILESname,
                                self._produce(semaphore, ft, ()), 10)
            for name, atoms, bonds in tqdm(results, total=self._hmmit, desc="Indentify isomers", unit="molecule"):
                mname.append(name)
                print(name, ",".join([str(x) for x in atoms]), ";".join(
                    [",".join([str(y) for y in x]) for x in bonds]), file=fm)
                semaphore.release()
        self._mname = mname

    def _convertSMILES(self, atoms, bonds):
        m = Chem.RWMol(Chem.MolFromSmiles(''))
        d = {}
        for atomnumber in atoms:
            d[atomnumber] = m.AddAtom(
                Chem.Atom(self.atomname[self._atomtype[atomnumber-1]-1]))
        for atom1, atom2, level in bonds:
            m.AddBond(d[atom1], d[atom2], Chem.BondType.DOUBLE if level == 2 else (
                Chem.BondType.TRIPLE if level == 3 else (Chem.BondType.AROMATIC if level == 12 else Chem.BondType.SINGLE)))
        name = Chem.MolToSmiles(m)
        return name

    @classmethod
    def _compress(cls, x):
        return base64.a85encode(zlib.compress(x.encode()))+b'\n'

    @classmethod
    def _decompress(cls, x):
        return zlib.decompress(base64.a85decode(x.strip())).decode()

    def _getatomeach(self):
        atomeach = np.zeros((self._N, self._step), dtype=np.int)
        with open(self.hmmfilename if self.runHMM else self.originfilename, 'rb') as fh, open(self.moleculetemp2filename, 'rb') as ft:
            for i, (linehz, linetz) in enumerate(zip(fh, ft), start=1):
                lineh = self._decompress(linehz)
                linet = self._decompress(linetz)
                s = linet.split()
                key1 = np.array([int(x) for x in s[0].split(",")])
                index = np.array(
                    [j for j in range(len(lineh)) if lineh[j] == "1"])
                if index.size:
                    atomeach[key1[:, None]-1, index] = i
        return atomeach

    def _getallroute(self, allmoleculeroute):
        allroute = Counter()
        for moleculeroute in allmoleculeroute:
            leftname = self._mname[moleculeroute[0]-1]
            rightname = self._mname[moleculeroute[1]-1]
            if leftname == rightname:
                continue
            equation = (leftname, rightname)
            allroute[equation] += 1
        return allroute

    def _searchspecies(self, originspec, sortedreactions, species):
        searchedspecies = []
        for reaction, n_reaction in sortedreactions:
            ii = 1
            if originspec == reaction[1-ii]:
                if not reaction[ii] in species:
                    searchedspecies.append(
                        (reaction[ii], (reaction, n_reaction)))
            if len(searchedspecies) >= self.n_searchspecies:
                break
        return searchedspecies

    def _printtable(self, allroute, maxsize=100):
        species = []
        table = np.zeros((maxsize, maxsize), dtype=np.int)
        reactionnumber = np.zeros((2), dtype=np.int)
        sortedreactions = sorted(
            allroute.items(), key=lambda d: d[1], reverse=True)
        # added on Nov 17, 2018
        if self.speciescenter:
            newreactions = []
            species = [self.speciescenter]
            newspecies = [self.speciescenter]
            while len(species) < maxsize and newspecies:
                newnewspecies = []
                for newspec in newspecies:
                    searchedspecies = self._searchspecies(
                        newspec, sortedreactions, species)
                    for searchedspec, searchedreaction in searchedspecies:
                        if len(species) < maxsize:
                            newnewspecies.append(searchedspec)
                            species.append(searchedspec)
                            newreactions.append(searchedreaction)
                newspecies = newnewspecies
            for reac in sortedreactions:
                if reac not in newreactions:
                    newreactions.append(reac)
            sortedreactions = newreactions

        with open(self.reactionfilename, 'w') as f:
            for reaction, n_reaction in sortedreactions:
                print(n_reaction, "->".join(reaction), file=f)
                for i, spec in enumerate(reaction):
                    if spec in species:
                        number = species.index(spec)
                    elif len(species) < 100:
                        species.append(spec)
                        number = species.index(spec)
                    else:
                        number = -1
                    reactionnumber[i] = number
                if all(reactionnumber >= 0):
                    table[reactionnumber[0]][reactionnumber[1]] = n_reaction

        with open(self.tablefilename, 'w') as f:
            print("\t"+"\t".join(species), file=f)
            for i, speciesi in enumerate(species):
                print(speciesi, end='\t', file=f)
                for j in range(len(species)):
                    print(table[i][j], end='\t', file=f)
                print(file=f)

    def _convertstructure(self, atoms, bonds):
        atomtypes = []
        for i, atom in enumerate(atoms, start=1):
            atomtypes.append((i, self.atomname.index(atom)))
        G = self._makemoleculegraph(atomtypes, bonds)
        return G

    @classmethod
    def _setparam(cls, x, default):
        return x if x else default

    def _setfilename(self, name, suffix):
        return self._setparam(name, f"{self.inputfilename}.{suffix}")

    def _printspecies(self):
        with open(self.moleculetemp2filename, 'rb') as f2, open(self.speciesfilename, 'w') as fw:
            d = [Counter() for i in range(len(self._timestep))]
            for name, line2 in zip(self._mname, f2):
                for t in [int(x) for x in self._decompress(line2).split()[-1].split(",")]:
                    d[t][name] += 1
            for t in range(len(self._timestep)):
                buff=[f"Timestep {self._timestep[t]}:"]
                buff.extend([f"{name} {num}" for name, num in d[t].items()])
                buff.append('\n')
                fw.write(' '.join(buff))

    def __enter__(self):
        return self

    def __exit__(self, Type, value, traceback):
        pass

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
