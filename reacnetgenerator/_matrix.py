# SPDX-License-Identifier: LGPL-3.0-or-later
# cython: language_level=3
# cython: linetrace=True
"""Generate Matrix.

A reaction network cannot accommodate too many species, so only the first
species which have the most reactions are taken. A reaction matrix can be
generated.
"""

import itertools
import operator
from collections import Counter
from typing import List

import numpy as np
import pandas as pd

from .utils import SharedRNGData, WriteBuffer, bytestolist, read_compressed_block


class _GenerateMatrix(SharedRNGData):
    tablefilename: str
    speciesfilename: str
    reactionfilename: str
    moleculetemp2filename: str
    n_searchspecies: int
    needprintspecies: bool
    allmoleculeroute: np.ndarray
    speciescenter: str
    matrix_size: int
    mname: np.ndarray
    timestep: np.ndarray
    splitmoleculeroute: List[np.ndarray]

    def __init__(self, rng):
        SharedRNGData.__init__(
            self,
            rng,
            [
                "tablefilename",
                "speciesfilename",
                "reactionfilename",
                "moleculetemp2filename",
                "n_searchspecies",
                "needprintspecies",
                "allmoleculeroute",
                "speciescenter",
                "mname",
                "timestep",
                "splitmoleculeroute",
                "matrix_size",
            ],
            [],
        )

    def generate(self):
        """Generate a reaction matrix and print species.

        A reaction matrix can be generated as
            R=[a_ij ], i=1,2,…,100;j=1,2,…,100
        where aij is the number of reactions from species si to sj.
        """
        self._printtable(self._getallroute(self.allmoleculeroute))
        if self.splitmoleculeroute is not None:
            for i, smr in enumerate(self.splitmoleculeroute):
                self._printtable(self._getallroute(smr), timeaxis=i)
        if self.needprintspecies:
            self._printspecies()

    def _getallroute(self, allmoleculeroute):
        names = self.mname[allmoleculeroute - 1]
        names = names[names[:, 0] != names[:, 1]]
        if names.size > 0:
            equations = np.unique(names, return_counts=True, axis=0)
            return zip(equations[0].tolist(), equations[1].tolist())
        return []

    def _printtable(self, allroute, timeaxis=None):
        maxsize = self.matrix_size
        species = []
        sortedreactions = sorted(allroute, key=operator.itemgetter(1, 0), reverse=True)
        # added on Nov 17, 2018
        if self.speciescenter:
            newreactions = []
            species = [self.speciescenter]
            newspecies = [self.speciescenter]
            while len(species) < maxsize and newspecies:
                newnewspecies = []
                for newspec in newspecies:
                    searchedspecies = self._searchspecies(
                        newspec, sortedreactions, species
                    )
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

        table = np.zeros((maxsize, maxsize), dtype=int)
        reactionnumber = np.zeros((2), dtype=int)
        with open(
            self.reactionfilename
            if timeaxis is None
            else f"{self.reactionfilename}.{timeaxis}",
            "w",
        ) as f:
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

        df = pd.DataFrame(
            table[: len(species), : len(species)], index=species, columns=species
        )
        df.to_csv(
            self.tablefilename
            if timeaxis is None
            else f"{self.tablefilename}.{timeaxis}",
            sep=" ",
        )

    def _searchspecies(self, originspec, sortedreactions, species):
        searchedspecies = []
        for reaction, n_reaction in sortedreactions:
            ii = 1
            if originspec == reaction[1 - ii]:
                if reaction[ii] not in species:
                    searchedspecies.append((reaction[ii], (reaction, n_reaction)))
            if len(searchedspecies) >= self.n_searchspecies:
                break
        return searchedspecies

    def _printspecies(self):
        with open(self.moleculetemp2filename, "rb") as ft, WriteBuffer(
            open(self.speciesfilename, "w")
        ) as fw:
            d = [Counter() for i in range(len(self.timestep))]
            for name, line in zip(
                self.mname, itertools.zip_longest(*[read_compressed_block(ft)] * 4)
            ):
                for t in bytestolist(line[-1]).tolist():
                    d[t][name] += 1
            for t in range(len(self.timestep)):
                fw.append(f"Timestep {self.timestep[t]}:")
                fw.extend(f" {item[0]} {item[1]}" for item in d[t].items())
                fw.append("\n")
