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


import itertools
import argparse
from functools import reduce
from multiprocessing import Pool, Semaphore, cpu_count
import math
import gc
import time
import zlib
import base64
from io import StringIO
from collections import Counter, defaultdict
import numpy as np
import networkx as nx
import networkx.algorithms.isomorphism as iso
import matplotlib as mpl
import matplotlib.pyplot as plt
from hmmlearn import hmm
from rdkit.Chem import Draw
from rdkit import Chem
import openbabel
from ase import Atom, Atoms
import scour.scour
from ._reachtml import _HTMLResult

plt.switch_backend('Agg')


class ReacNetGenerator(object):
    ''' Use ReacNetGenerator for trajectory analysis'''

    def __init__(self, inputfiletype='lammpsbondfile', inputfilename='bonds.reaxc', atomname=["C", "H", "O"], selectatoms=None, originfilename=None, hmmfilename=None, atomfilename=None, moleculefilename=None, atomroutefilename=None, reactionfilename=None, tablefilename=None, moleculetempfilename=None, moleculetemp2filename=None, moleculestructurefilename=None, imagefilename=None, speciesfilename=None, resultfilename=None, stepinterval=1, p=[0.5, 0.5], a=[[0.999, 0.001], [0.001, 0.999]], b=[[0.6, 0.4], [0.4, 0.6]], runHMM=True, SMILES=True, getoriginfile=False, species={}, node_size=200, font_size=6, widthcoefficient=1,  maxspecies=20, nolabel=False, needprintspecies=True, filter=[], node_color=[78/256, 196/256, 238/256], pos={}, printfiltersignal=False, showid=True, k=None, start_color=[0, 0, 1], end_color=[1, 0, 0], nproc=None, speciescenter=None, n_searchspecies=2):
        ''' Init ReacNetGenerator '''
        print(__doc__)
        print(
            f"Version: {__version__}  Creation date: {__date__}  Update date: {__update__}")
        self.inputfiletype = inputfiletype
        self.inputfilename = inputfilename
        self.atomname = atomname
        self.selectatoms = self._setparam(selectatoms, self.atomname)
        self.originfilename = self._setfilename(originfilename, "origin")
        self.hmmfilename = self._setfilename(hmmfilename, "hmm")
        self.atomfilename = self._setfilename(atomfilename, "atom")
        self.moleculefilename = self._setfilename(moleculefilename, "moname")
        self.atomroutefilename = self._setfilename(atomroutefilename, "route")
        self.reactionfilename = self._setfilename(reactionfilename, "reaction")
        self.tablefilename = self._setfilename(tablefilename, "table")
        self.moleculetempfilename = self._setfilename(
            moleculetempfilename, "temp")
        self.moleculetemp2filename = self._setfilename(
            moleculetemp2filename, "temp2")
        self.moleculestructurefilename = self._setfilename(
            moleculestructurefilename, "structure")
        self.imagefilename = self._setfilename(imagefilename, "svg")
        self.speciesfilename = self._setfilename(speciesfilename, "species")
        self.resultfilename = self._setfilename(resultfilename, "html")
        self.stepinterval = stepinterval
        self.p = np.array(p)
        self.a = np.array(a)
        self.b = np.array(b)
        self.runHMM = runHMM
        self.SMILES = SMILES
        self.getoriginfile = getoriginfile if self.runHMM else True
        self.species = species
        self.needprintspecies = needprintspecies
        self.node_size = node_size
        self.font_size = font_size
        self.widthcoefficient = widthcoefficient
        self.maxspecies = maxspecies
        self.nolabel = nolabel
        self.filter = filter
        self.node_color = np.array(node_color)
        self.pos = pos
        self.printfiltersignal = printfiltersignal
        self.showid = showid
        self.k = k
        self.start_color = np.array(start_color)
        self.end_color = np.array(end_color)
        self.nproc = self._setparam(nproc, cpu_count())
        self.loggingfreq = 1000
        self.speciescenter = speciescenter
        self.n_searchspecies = n_searchspecies
        self._timearray = []
        self._statusid = 0
        self._statusidmax = 0

    def runanddraw(self, run=True, draw=True, report=True):
        ''' Analyze the trajectory from MD simulation '''
        self._statusidmax = max(
            self._statusidmax, (6 if report else (5 if draw else 4)))
        self._printtime(0)
        if run:
            self.run()
        if draw:
            self.draw()
        if report:
            self.report()

    def run(self):
        """ Processing of MD trajectory """
        self._statusidmax = max(self._statusidmax, 4)
        self._printtime(0)
        for runstep in range(1, 5):
            if runstep == 1:
                self._readinputfile()
            elif runstep == 2:
                if self.runHMM:
                    self._initHMM()
                self._calhmm()
            elif runstep == 3:
                if self.SMILES:
                    self._printmoleculeSMILESname()
                else:
                    self._printmoleculename()
                atomeach = self._getatomeach()
                allmoleculeroute = self._printatomroute(atomeach)
            elif runstep == 4:
                allroute = self._getallroute(allmoleculeroute)
                self._printtable(allroute)
                if self.needprintspecies:
                    self._printspecies()
            # garbage collect
            gc.collect()
            self._printtime(runstep)

    def draw(self):
        """ Draw the reaction network """
        self._statusidmax = max(self._statusidmax, 5)
        self._printtime(0)
        # read table
        table, name = self._readtable()
        species, showname = self._handlespecies(name)

        G = nx.DiGraph()
        for i, tablei in enumerate(table):
            if name[i] in species and not name[i] in self.filter:
                G.add_node(showname[name[i]] if name[i]
                           in showname else name[i])
                for j, tableij in enumerate(tablei):
                    if name[j] in species and not name[j] in self.filter:
                        if tableij > 0:
                            G.add_weighted_edges_from([((showname[name[i]] if name[i] in showname else name[i]), (
                                showname[name[j]] if name[j] in showname else name[j]), tableij)])
        weights = np.array([math.log(G[u][v]['weight']+1)
                            for u, v in G.edges()])
        widths = [weight/max(weights) * self.widthcoefficient*2 if weight > max(weights)
                  * 0.7 else weight/max(weights) * self.widthcoefficient*0.5 for weight in weights]
        colors = [self.start_color + weight /
                  max(weights) * (self.end_color-self.start_color) for weight in weights]
        try:
            self.pos = (nx.spring_layout(G) if not self.pos else nx.spring_layout(G, pos=self.pos, fixed=[p for p in self.pos])) if not self.k else (
                nx.spring_layout(G, k=self.k) if not self.pos else nx.spring_layout(G, pos=self.pos, fixed=[p for p in self.pos], k=self.k))
            if self.pos:
                self._logging("The position of the species in the network is:")
                self._logging(self.pos)
            for with_labels in ([True] if not self.nolabel else [True, False]):
                nx.draw(G, pos=self.pos, width=widths, node_size=self.node_size, font_size=self.font_size,
                        with_labels=with_labels, edge_color=colors, node_color=self.node_color)
                imagefilename = "".join(
                    (("" if with_labels else "nolabel_"), self.imagefilename))
                with StringIO() as stringio, open(imagefilename, 'w') as f:
                    plt.savefig(stringio, format='svg')
                    f.write(scour.scour.scourString(stringio.getvalue()))
                plt.close()
        except Exception as e:
            self._logging("Error: cannot draw images. Details:", e)
        self._printtime(5)

    def report(self):
        """ Generate the analysis report """
        self._statusidmax = max(self._statusidmax, 6)
        self._printtime(0)
        _HTMLResult(self)._report()
        self._logging(
            f"Report is generated. Please see {self.resultfilename} for more details.")
        self._printtime(6)

    def _logging(self, *message, end='\n'):
        if message:
            localtime = time.asctime(time.localtime(time.time()))
            print(f"{localtime} ReacNetGenerator {__version__}",
                  *message, end=end)
        else:
            print(end=end)

    def _loggingprocessing(self, index):
        if index % self.loggingfreq == 0:
            self._logging(f"processing {index} ...", end='\r')

    @property
    def _status(self):
        return ["Init", "Read bond information and Detect molecules", "HMM filter", "Indentify isomers and collect reaction paths", "Reaction matrix generation", "Draw reaction network", "Generate analysis report"][self._statusid]

    def _printtime(self, statusid):
        self._statusid = statusid
        if not self._timearray or self._statusid > 0:
            self._timearray.append(time.time())
            if statusid > 0:
                self._logging(
                    f"Step {len(self._timearray)-1}: Done! Time consumed (s): {self._timearray[-1]-self._timearray[-2]:.3f} ({self._status})")
            if statusid >= self._statusidmax:
                self._logging("====== Summary ======")
                self._logging(
                    f"Total time(s): {self._timearray[-1]-self._timearray[0]:.3f} s")

    def _mo(self, i, bond, level, molecule, done, bondlist):
        # connect molecule
        molecule.append(i)
        done[i-1] = True
        for b, l in zip(bond[i-1], level[i-1]):
            if not done[b-1]:
                bondlist.append((i, b, l) if i < b else (b, i, l))
                molecule, done, bondlist = self._mo(
                    b, bond, level, molecule, done, bondlist)
        return molecule, done, bondlist

    @property
    def _readNfunc(self):
        if self.inputfiletype == "lammpsbondfile":
            return self._readlammpsbondN
        return self._readlammpscrdN

    @property
    def _readstepfunc(self):
        if self.inputfiletype == "lammpsbondfile":
            return self._readlammpsbondstep
        return self._readlammpscrdstep

    def _readinputfile(self):
        self._steplinenum = self._readNfunc()
        self._getdandtimestep()

    def _readlammpsbondN(self):
        with open(self.inputfilename) as file:
            iscompleted = False
            for index, line in enumerate(file):
                if line.startswith("#"):
                    if line.startswith("# Number of particles"):
                        if iscompleted:
                            stepbindex = index
                            break
                        else:
                            iscompleted = True
                            stepaindex = index
                        N = [int(s) for s in line.split() if s.isdigit()][0]
                        atomtype = np.zeros(N, dtype=np.int)
                else:
                    s = line.split()
                    atomtype[int(s[0])-1] = int(s[1])
        steplinenum = stepbindex-stepaindex
        self._N = N
        self._atomtype = atomtype
        return steplinenum

    def _readlammpsbondstep(self, item):
        (step, lines), _ = item
        bond = [None for x in range(self._N)]
        level = [None for x in range(self._N)]
        for line in lines:
            if line:
                if line.startswith("#"):
                    if line.startswith("# Timestep"):
                        timestep = int(line.split()[-1])
                else:
                    s = line.split()
                    bond[int(s[0])-1] = [int(x) for x in s[3:3+int(s[2])]]
                    level[int(s[0])-1] = [max(1, round(float(x)))
                                          for x in s[4+int(s[2]):4+2*int(s[2])]]
        molecules = self._connectmolecule(bond, level)
        return molecules, (step, timestep)

    def _readlammpscrdN(self):
        with open(self.inputfilename) as f:
            iscompleted = False
            for index, line in enumerate(f):
                if line.startswith("ITEM:"):
                    linecontent = 4 if line.startswith("ITEM: TIMESTEP") else (3 if line.startswith(
                        "ITEM: ATOMS") else (1 if line.startswith("ITEM: NUMBER OF ATOMS") else 2))
                else:
                    if linecontent == 1:
                        if iscompleted:
                            stepbindex = index
                            break
                        else:
                            iscompleted = True
                            stepaindex = index
                        N = int(line.split()[0])
                        atomtype = np.zeros(N, dtype=np.int)
                    elif linecontent == 3:
                        s = line.split()
                        atomtype[int(s[0])-1] = int(s[1])
        steplinenum = stepbindex-stepaindex
        self._N = N
        self._atomtype = atomtype
        return steplinenum

    def _readlammpscrdstep(self, item):
        (step, lines), _ = item
        step_atoms = []
        for line in lines:
            if line:
                if line.startswith("ITEM:"):
                    linecontent = 4 if line.startswith("ITEM: TIMESTEP") else (3 if line.startswith(
                        "ITEM: ATOMS") else (1 if line.startswith("ITEM: NUMBER OF ATOMS") else 2))
                else:
                    if linecontent == 3:
                        s = line.split()
                        step_atoms.append(
                            (int(s[0]), Atom(self.atomname[int(s[1])-1], [float(x) for x in s[2:5]])))
                    elif linecontent == 4:
                        timestep = step, int(line.split()[0])
        _, step_atoms = zip(*sorted(step_atoms, key=lambda a: a[0]))
        step_atoms = Atoms(step_atoms)
        bond, level = self._getbondfromcrd(step_atoms)
        molecules = self._connectmolecule(bond, level)
        return molecules, timestep

    def _getdandtimestep(self):
        d = defaultdict(list)
        timestep = {}
        with open(self.inputfilename) as file, Pool(self.nproc, maxtasksperchild=1000) as pool:
            semaphore = Semaphore(360)
            results = pool.imap_unordered(self._readstepfunc, self._produce(semaphore, enumerate(itertools.islice(
                itertools.zip_longest(*[file]*self._steplinenum), 0, None, self.stepinterval)), None), 10)
            for index, (molecules, (step, thetimestep)) in enumerate(results):
                self._loggingprocessing(index)
                for molecule in molecules:
                    d[molecule].append(step)
                timestep[step] = thetimestep
                semaphore.release()
        self._writemoleculetempfile(d)
        self._timestep = timestep
        self._step = len(timestep)-1

    def _connectmolecule(self, bond, level):
        molecules = []
        done = np.zeros(self._N, dtype=bool)
        for i in range(1, self._N+1):
            if not done[i-1]:
                mole, done, bondlist = self._mo(i, bond, level, [], done, [])
                moleculestr = ' '.join((",".join((str(x) for x in sorted(mole))), ";".join(
                    (",".join([str(y) for y in x]) for x in sorted(bondlist)))))
                molecules.append(self._compress(moleculestr))
        return molecules

    def _writemoleculetempfile(self, d):
        with open(self.moleculetempfilename, 'wb') as f:
            for key, value in d.items():
                f.write(self._compress(
                    ' '.join((self._decompress(key), ",".join((str(x) for x in value))))))

    def _getbondfromcrd(self, step_atoms):
        atomnumber = len(step_atoms)
        xyzstring = f"{atomnumber}\nReacNetGenerator\n"+"\n".join(
            [f'{s:2s} {x:22.15f} {y:22.15f} {z:22.15f}' for s, (x, y, z) in zip(step_atoms.get_chemical_symbols(), step_atoms.positions)])
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats('xyz', 'mol2')
        mol = openbabel.OBMol()
        conv.ReadString(mol, xyzstring)
        mol2string = conv.WriteString(mol)
        linecontent = -1
        bond = [[] for i in range(atomnumber)]
        bondlevel = [[] for i in range(atomnumber)]
        for line in mol2string.split('\n'):
            if line.startswith("@<TRIPOS>BOND"):
                linecontent = 0
            else:
                if linecontent == 0:
                    s = line.split()
                    if len(s) > 3:
                        bond[int(s[1])-1].append(int(s[2]))
                        bond[int(s[2])-1].append(int(s[1]))
                        level = 12 if s[3] == 'ar' else int(s[3])
                        bondlevel[int(s[1])-1].append(level)
                        bondlevel[int(s[2])-1].append(level)
        return bond, bondlevel

    def _initHMM(self):
        self._model = hmm.MultinomialHMM(n_components=2)
        self._model.startprob_ = self.p
        self._model.transmat_ = self.a
        self._model.emissionprob_ = self.b

    def _produce(self, semaphore, plist, parameter):
        for item in plist:
            semaphore.acquire()
            yield item, parameter

    def _getoriginandhmm(self, item):
        line_c, _ = item
        line = self._decompress(line_c)
        s = line.split()
        value = np.array([int(x)-1 for x in s[-1].split(",")])
        origin = np.zeros(self._step, dtype=np.int)
        origin[value] = 1
        if self.runHMM:
            _, hmm = self._model.decode(
                np.array([origin]).T, algorithm="viterbi")
        return origin, (np.array(hmm) if self.runHMM else np.array([])), line

    def _calhmm(self):
        with open(self.originfilename, 'wb') if self.getoriginfile or not self.runHMM else Placeholder() as fo, open(self.hmmfilename, 'wb') if self.runHMM else Placeholder() as fh, open(self.moleculetempfilename, 'rb') as ft, open(self.moleculetemp2filename, 'wb') as ft2, Pool(self.nproc, maxtasksperchild=1000) as pool:
            semaphore = Semaphore(360)
            results = pool.imap_unordered(
                self._getoriginandhmm, self._produce(semaphore, ft, ()), 10)
            for index, (originsignal, hmmsignal, mlist) in enumerate(results):
                self._loggingprocessing(index)
                if 1 in hmmsignal or self.printfiltersignal or not self.runHMM:
                    if self.getoriginfile:
                        fo.write(self._compress(
                            "".join([str(i) for i in originsignal.tolist()])))
                    if self.runHMM:
                        fh.write(self._compress(
                            "".join([str(i) for i in hmmsignal.tolist()])))
                    ft2.write(self._compress(mlist.strip()))
                semaphore.release()

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
            semaphore = Semaphore(360)
            results = pool.imap(self._getatomroute, self._produce(
                semaphore, enumerate(zip(atomeach, self._atomtype), start=1), ()), 10)
            for index, route in enumerate(results):
                self._loggingprocessing(index)
                moleculeroute, routestr = route
                print(routestr, file=f)
                for mroute in moleculeroute:
                    if not mroute in allmoleculeroute:
                        allmoleculeroute.append(mroute)
                semaphore.release()
        return allmoleculeroute

    def _makemoleculegraph(self, atoms, bonds):
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
        return name+" "+",".join([self.atomname[self._atomtype[x-1]-1] for x in atoms])+" "+";".join([str(index[x[0]])+","+str(index[x[1]])+","+str(x[2]) for x in bonds])

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
            semaphore = Semaphore(360)
            results = pool.imap(self._calmoleculeSMILESname,
                                self._produce(semaphore, ft, ()), 10)
            for index, (name, atoms, bonds) in enumerate(results):
                self._loggingprocessing(index)
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

    def _compress(self, x):
        return base64.a85encode(zlib.compress(x.encode()))+b'\n'

    def _decompress(self, x):
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
                if(len(index)) > 0:
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
                if not reac in newreactions:
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
            for i in range(len(species)):
                print(species[i], end='\t', file=f)
                for j in range(len(species)):
                    print(table[i][j], end='\t', file=f)
                print(file=f)

    def _readtable(self):
        table = []
        name = []
        with open(self.tablefilename) as file:
            for line in itertools.islice(file, 1, None):
                name.append(line.split()[0])
                table.append([int(s) for s in line.split()[1:]])
        return np.array(table), name

    def _convertstructure(self, atoms, bonds):
        atomtypes = []
        for i, atom in enumerate(atoms, start=1):
            atomtypes.append((i, self.atomname.index(atom)))
        G = self._makemoleculegraph(atomtypes, bonds)
        return G

    def _handlespecies(self, name):
        showname = {}
        if self.species == {}:
            species_out = dict([(x, {}) for x in (name if len(
                name) <= self.maxspecies else name[0:self.maxspecies])])
        else:
            species_out = {}
            b = True
            for spec in self.species.items():
                specname, value = spec
                if "structure" in value:
                    atoms, bonds = value["structure"]
                    G1 = self._convertstructure(atoms, bonds)
                    if b:
                        structures = self._readstrcture()
                        em = iso.numerical_edge_match(
                            ['atom', 'level'], ["None", 1])
                        b = False
                    i = 1
                    while (specname+"_"+str(i) if i > 1 else specname) in structures:
                        G2 = self._convertstructure(structures[(
                            specname+"_"+str(i) if i > 1 else specname)][0], structures[(specname+"_"+str(i) if i > 1 else specname)][1])
                        if nx.is_isomorphic(G1, G2, em):
                            if i > 1:
                                specname += "_"+str(i)
                            break
                        i += 1
                species_out[specname] = {}
                if "showname" in value:
                    showname[specname] = value["showname"]
        if self.showid:
            if species_out:
                self._logging()
                self._logging("Species are:")
                for n, (specname, value) in enumerate(species_out.items(), start=1):
                    showname[specname] = str(n)
                    print(n, specname)
        return species_out, showname

    def _setparam(self, x, default):
        return x if x else default

    def _setfilename(self, name, suffix):
        return self._setparam(name, self.inputfilename+"."+suffix)

    def _printspecies(self):
        with open(self.moleculetemp2filename, 'rb') as f2, open(self.speciesfilename, 'w') as fw:
            d = [Counter() for i in range(len(self._timestep))]
            for name, line2 in zip(self._mname, f2):
                for t in [int(x) for x in self._decompress(line2).split()[-1].split(",")]:
                    d[t][name] += 1
            for t in range(len(self._timestep)):
                print("Timestep", self._timestep[t], ":", end=' ', file=fw)
                for name, num in d[t].items():
                    print(name, num, end=' ', file=fw)
                print(file=fw)

    def __enter__(self):
        return self

    def __exit__(self, Type, value, traceback):
        pass


class Placeholder(object):
    def __init__(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, Type, value, traceback):
        pass


def _commandline():
    parser = argparse.ArgumentParser(description='ReacNetGenerator')
    parser.add_argument('-i', '--inputfilename',
                        help='Input trajectory file, e.g. bonds.reaxc', required=True)
    parser.add_argument(
        '-a', '--atomname', help='Atomic names in the trajectory, e.g. C H O', nargs='*', required=True)
    parser.add_argument(
        '--nohmm', help='Process trajectory without Hidden Markov Model (HMM)', action="store_true")
    parser.add_argument(
        '--dump', help='Process the LAMMPS dump file', action="store_true")
    parser.add_argument(
        '-n', '-np', '--nproc', help='Number of processes')
    parser.add_argument('-s', '--selectatoms',
                        help='Select atoms in the reaction network, e.g. C', nargs='*')
    args = parser.parse_args()
    ReacNetGenerator(inputfilename=args.inputfilename, atomname=args.atomname, runHMM=not args.nohmm,
                     inputfiletype=('lammpsdumpfile' if args.dump else 'lammpsbondfile'), nproc=args.nproc, selectatoms=args.selectatoms).runanddraw()
