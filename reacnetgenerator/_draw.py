# cython: language_level=3
# cython: linetrace=True
"""Draw reaction network.

With the reaction matrix, a reaction network can be drawn. Here the NetworkX
package is used to make a graph which indicates reactions between species.
Fruchterman-Reingold force-directed algorithm is used to make layout of nodes
relate to reaction quantity and different colors and widths of lines are drawn
depending on reaction quantity. The distance of two species in the network,
the color and thickness of the line between them are determined by the number
of their reactions, making the reaction network more intuitive.

Reference:
[1] Hagberg, A.; Swart, P.; Daniel, S. C. Exploring network structure,
dynamics, and function using NetworkX; Los Alamos National Lab. (LANL), Los
Alamos, NM (United States): 2008.
[2] Fruchterman, T. M.; Reingold, E. M. Graph drawing by force-directed
placement. Software: Practice and experince. 1991, 21(11),1129-1164.
"""

import logging
import math
import traceback
from io import StringIO
from itertools import permutations

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scour.scour

from .utils import SCOUROPTIONS, SharedRNGData


class _DrawNetwork(SharedRNGData):
    def __init__(self, rng):
        SharedRNGData.__init__(self, rng, ['atomname', 'tablefilename', 'imagefilename', 'maxspecies', 'species', 'speciesfilter',
                                  'start_color', 'end_color', 'node_size', 'node_color', 'font_size', 'widthcoefficient', 'k', 'pos', 'nolabel',
                                  'showid', 'split'], [])

    def draw(self):
        """Draw the network."""
        self._draw()
        if self.split > 1:
            for st in range(self.split):
                self._draw(timeaxis=st)

    def _draw(self, timeaxis=None):
        table, name = self._readtable(
            self.tablefilename if timeaxis is None else f"{self.tablefilename}.{timeaxis}")
        species, showname = self._handlespecies(name)

        G = nx.DiGraph()
        idx = []
        for i, _ in enumerate(table):
            if name[i] in species:
                G.add_node(showname[name[i]])
                idx.append(i)
        for i, j in permutations(idx, 2):
            if table[i][j] > 0:
                G.add_weighted_edges_from(
                    [(showname[name[i]], showname[name[j]], table[i][j])])
        weights = np.array([math.log(G[u][v]['weight']+1)
                            for u, v in G.edges()])
        widths = weights / np.max(weights) * self.widthcoefficient * np.array(
            [0.5, 2])[(weights > np.max(weights) * 0.7)+0] if weights.size else np.zeros(0)
        colors = self.start_color + weights[:, np.newaxis] / np.max(weights) * (
            self.end_color - self.start_color) if weights.size else np.zeros(0)
        try:
            pos = nx.spring_layout(G,
                                   pos=self.pos if self.pos else None,
                                   fixed=list(self.pos) if self.pos else None,
                                   k=self.k)
            if pos:
                logging.info("The position of the species in the network is:")
                logging.info(pos)
            for with_labels in ([True] if not self.nolabel else [True, False]):
                nx.draw(
                    G, pos=pos, width=widths, node_size=self.node_size,
                    font_size=self.font_size, with_labels=with_labels,
                    edge_color=colors, node_color=[self.node_color]*len(pos))
                imagefilename = "".join(
                    (("" if with_labels else "nolabel_"), self.imagefilename))
                with StringIO() as stringio, open(imagefilename if timeaxis is None else f"{imagefilename}.{timeaxis}", 'w') as f:
                    plt.savefig(stringio, format='svg')
                    f.write(scour.scour.scourString(
                        stringio.getvalue(), SCOUROPTIONS))
                plt.close()
        except Exception as e:
            logging.error(f"Error: cannot draw images. Details: {e}")
            traceback.print_tb(e.__traceback__)

    def _readtable(self, filename):
        df = pd.read_csv(filename, sep=' ', index_col=0, header=0)
        return df.values, df.index

    def _handlespecies(self, name):
        species = self.species if self.species else name[:min(
            len(name), self.maxspecies)]
        # filter
        species = [spec for spec in species if spec not in self.speciesfilter]

        if self.showid:
            showname = dict([(v, u) for u, v in enumerate(species, start=1)])
        else:
            showname = dict([(u, u) for u in species])
        if species:
            logging.info("Species are:")
            for specname, n in showname.items():
                logging.info("{} {}".format(n, specname))
        return species, showname
