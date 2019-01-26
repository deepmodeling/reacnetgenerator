''' Draw reaction network '''

import logging
import math
from io import StringIO
import itertools

import networkx.algorithms.isomorphism as iso
import networkx as nx
import numpy as np
import scour.scour
import matplotlib.pyplot as plt
import pandas as pd

from ._path import _CollectMolPaths


class _DrawNetwork:
    def __init__(self, rng):
        self.atomname = rng.atomname
        self.tablefilename = rng.tablefilename
        self.imagefilename = rng.imagefilename
        self.moleculestructurefilename = rng.moleculestructurefilename
        self.maxspecies = rng.maxspecies
        self.species = rng.species
        self.speciesfilter = rng.speciesfilter
        self.start_color = rng.start_color
        self.end_color = rng.end_color
        self.node_size = rng.node_size
        self.node_color = rng.node_color
        self.font_size = rng.font_size
        self.widthcoefficient = rng.widthcoefficient
        self.k = rng.k
        self.pos = rng.pos
        self.nolabel = rng.nolabel
        self.showid = rng.showid

    def draw(self):
        ''' Draw the network '''
        table, name = self._readtable()
        species, showname = self._handlespecies(name)

        G = nx.DiGraph()
        for i, tablei in enumerate(table):
            if name[i] in species and not name[i] in self.speciesfilter:
                G.add_node(showname[name[i]] if name[i]
                           in showname else name[i])
                for j, tableij in enumerate(tablei):
                    if name[j] in species and not name[j] in self.speciesfilter:
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
                logging.info("The position of the species in the network is:")
                logging.info(self.pos)
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
            logging.error(f"Error: cannot draw images. Details: {e}")

    def _readtable(self):
        df = pd.read_csv(self.tablefilename, sep=' ', index_col=0, header=0)
        return df.values, df.index

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
                    while (f"{specname}_{i}" if i > 1 else specname) in structures:
                        G2 = self._convertstructure(structures[(
                            f"{specname}_{i}" if i > 1 else specname)][0], structures[(f"{specname}_{i}" if i > 1 else specname)][1])
                        if nx.is_isomorphic(G1, G2, em):
                            if i > 1:
                                specname += f"_{i}"
                            break
                        i += 1
                species_out[specname] = {}
                if "showname" in value:
                    showname[specname] = value["showname"]
        if self.showid:
            if species_out:
                print()
                logging.info("Species are:")
                for n, (specname, value) in enumerate(species_out.items(), start=1):
                    showname[specname] = str(n)
                    print(n, specname)
        return species_out, showname

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

    def _convertstructure(self, atoms, bonds):
        atomtypes = []
        for i, atom in enumerate(atoms, start=1):
            atomtypes.append((i, self.atomname.index(atom)))
        G = _CollectMolPaths._makemoleculegraph(atomtypes, bonds)
        return G
