"""Generate Matrix.

A reaction network cannot accommodate too many species, so only the first
species which have the most reactions are taken. A reaction matrix can be
generated.
"""

import itertools
from collections import Counter

import numpy as np
import pandas as pd


class _GenerateMatrix:
    def __init__(self, rng):
        self.rng = rng
        self.tablefilename = rng.tablefilename
        self.speciesfilename = rng.speciesfilename
        self.reactionfilename = rng.reactionfilename
        self.moleculetemp2filename = rng.moleculetemp2filename
        self.n_searchspecies = rng.n_searchspecies
        self.needprintspecies = rng.needprintspecies
        self.allmoleculeroute = rng.allmoleculeroute
        self.speciescenter = rng.speciescenter
        self._mname = rng.mname
        self._timestep = rng.timestep
        self._decompress = rng.decompress
        self._bytestolist = rng.bytestolist

    def generate(self):
        """Generate a reaction matrix and print species.

        A reaction matrix can be generated as
            R=[a_ij ], i=1,2,…,100;j=1,2,…,100
        where aij is the number of reactions from species si to sj.
        """
        allroute = self._getallroute(self.allmoleculeroute)
        self._printtable(allroute)
        if self.needprintspecies:
            self._printspecies()

    def _getallroute(self, allmoleculeroute):
        allroute = Counter()
        names = self._mname[allmoleculeroute-1]
        names = names[names[:, 0] != names[:, 1]]
        if names.size > 0:
            equations = np.unique(names, return_counts=True, axis=0)
            return zip(equations[0].tolist(), equations[1].tolist())
        return []

    def _printtable(self, allroute, maxsize=100):
        species = []
        sortedreactions = sorted(
            allroute, key=lambda d: d[1], reverse=True)
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

        table = np.zeros((maxsize, maxsize), dtype=np.int)
        reactionnumber = np.zeros((2), dtype=np.int)
        with open(self.reactionfilename, 'w') as f:
            for reaction, n_reaction in sortedreactions:
                f.write(f"{n_reaction} {'->'.join(reaction)}\n")
                for i, spec in enumerate(reaction):
                    if spec in species:
                        number = species.index(spec)
                    elif len(species) < maxsize:
                        species.append(spec)
                        number = species.index(spec)
                    else:
                        number = -1
                    reactionnumber[i] = number
                if all(reactionnumber >= 0):
                    table[tuple(reactionnumber)] = n_reaction

        df = pd.DataFrame(table[:len(species), :len(
            species)], index=species, columns=species)
        df.to_csv(self.tablefilename, sep=' ')

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

    def _printspecies(self):
        with open(self.moleculetemp2filename, 'rb') as ft, open(self.speciesfilename, 'w') as fw:
            d = [Counter() for i in range(len(self._timestep))]
            for name, line in zip(self._mname, itertools.zip_longest(*[ft] * 3)):
                for t in self._bytestolist(line[-1], nparray=True).tolist():
                    d[t][name] += 1
            for t, ts in enumerate(self._timestep):
                buff = [f"Timestep {ts}:"]
                buff.extend(
                    map(lambda item: 'f" {item[0]} {item[1]}"', d[t].items()))
                buff.append('\n')
                if len(buff) > 200:
                    fw.write(''.join(buff))
            if buff:
                fw.write(''.join(buff))
