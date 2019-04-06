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

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scour.scour


class _DrawNetwork:
    def __init__(self, rng):
        self.atomname = rng.atomname
        self.tablefilename = rng.tablefilename
        self.imagefilename = rng.imagefilename
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
        self._split = rng.split

    def draw(self):
        """Draw the network."""
        self._draw()
        if self._split > 1:
            for st in range(self._split):
                self._draw(timeaxis=st)

    def _draw(self, timeaxis=None):
        table, name = self._readtable(
            self.tablefilename if timeaxis is None else f"{self.tablefilename}.{timeaxis}")
        species, showname = self._handlespecies(name)

        G = nx.DiGraph()
        for i, tablei in enumerate(table):
            if name[i] in species and not name[i] in self.speciesfilter:
                G.add_node(showname[name[i]] if name[i]
                           in showname else name[i])
                for j, tableij in enumerate(tablei):
                    if name[j] in species and not name[j] in self.speciesfilter:
                        if tableij > 0:
                            G.add_weighted_edges_from(
                                [((showname[name[i]]
                                   if name[i] in showname else name[i]),
                                  (showname[name[j]]
                                   if name[j] in showname else name[j]),
                                  tableij)])
        weights = np.array([math.log(G[u][v]['weight']+1)
                            for u, v in G.edges()])
        widths = weights / np.max(weights) * self.widthcoefficient * \
            np.array([0.5, 2])[(weights > np.max(weights) * 0.7)+0]
        colors = self.start_color + weights[:, np.newaxis] / \
            np.max(weights) * (self.end_color - self.start_color)
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
                    f.write(scour.scour.scourString(stringio.getvalue()))
                plt.close()
        except Exception as e:
            logging.error(f"Error: cannot draw images. Details: {e}")
            traceback.print_tb(e.__traceback__)

    def _readtable(self, filename):
        df = pd.read_csv(filename, sep=' ', index_col=0, header=0)
        return df.values, df.index

    def _handlespecies(self, name):
        showname = {}
        species = self.species if self.species.size > 0 else name[:min(
            len(name), self.maxspecies)]

        if self.showid:
            if species.size > 0:
                print()
                logging.info("Species are:")
                showname = dict([(v, u)
                                 for u, v in enumerate(species, start=1)])
                for specname, n in showname.items():
                    print(n, specname)
        return species, showname
