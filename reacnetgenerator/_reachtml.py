# SPDX-License-Identifier: LGPL-3.0-or-later
# cython: linetrace=True
# cython: language_level=3
"""Generate a web page to show the analysis report.

For the convenience analyze of the trajectories, we put all the results
generated by ReacNetGenerator (Network, Species, and Reactions) in an interactive
web page. By default, 20 species with the most reactions are taken to draw the
network. However, by clicking on a given species, one can check the special
network which starts from it.
"""

import json
import os
import re
from collections import defaultdict
from typing import Dict, Union

import numpy as np
import scour.scour

from ._logging import logger
from .utils import SCOUROPTIONS, SharedRNGData


class _HTMLResult(SharedRNGData):
    reactionfilename: str
    resultfilename: str
    imagefilename: str
    reactionabcdfilename: str
    jsonfilename: str
    nproc: int
    split: int
    atomname: np.ndarray

    def __init__(self, rng):
        SharedRNGData.__init__(
            self,
            rng,
            [
                "reactionfilename",
                "resultfilename",
                "imagefilename",
                "reactionabcdfilename",
                "jsonfilename",
                "nproc",
                "split",
                "atomname",
            ],
            [],
            ["_specs", "_reaction", "_reactionsabcd", "_svgfiles"],
        )
        self._templatedict: Dict[str, Union[str, int, list, dict]] = {
            "speciesshownum": 30,
            "reactionsshownum": 20,
        }
        self._linkreac = defaultdict(list)
        # define instance
        self._svgspecs = set()

    def report(self):
        """Generate a web page to show the result."""
        self._readdata()
        self._generateresult()
        logger.info(
            f"Report is generated. Please see {self.resultfilename} for more details."
        )

    def _re(self, smi):
        """If you use RDkit to convert a methyl radical to SMILES, you will get something
        like [H]C([H])[H]. However, OpenBabel will consider it as a methane molecule. So,
        you have to use [H][C]([H])[H], if you need to process some radicals.

        Examples
        --------
        >>> self._re('C')
        [C]
        >>> self._re('[C]')
        [C]
        >>> self._re('[CH]')
        [CH]
        >>> self._re('Na')
        [Na]
        >>> self._re('[H]c(Cl)C([H])Cl')
        [H][c]([Cl])[C]([H])[Cl]
        """
        if "_unknownSMILES" in smi:
            # not SMILES
            return smi
        Satom = sorted(self.atomname, key=lambda i: len(i), reverse=True)
        elements = "|".join(
            [
                ((an.upper() + "|" + an.lower()) if len(an) == 1 else an)
                for an in Satom
                if an != "H"
            ]
        )
        smi = re.sub(r"(?<!\[)(" + elements + r")(?!H)", r"[\1]", smi)
        return smi.replace("[HH]", "[H]")

    def _handlereaction(self, line):
        sx = line.split()
        left, right = sx[1].split("->")
        left = [self._re(spec) for spec in left.split("+")]
        right = [self._re(spec) for spec in right.split("+")]
        num = int(sx[0])
        return left, right, num

    def _readreaction(self, timeaxis=None, linknum=6):
        reaction = []
        with open(
            self.reactionfilename
            if timeaxis is None
            else f"{self.reactionfilename}.{timeaxis}"
        ) as f:
            for i, line in enumerate(f, 1):
                left, right, num = self._handlereaction(line)
                reaction.append({"i": i, "l": left, "r": right, "n": num})
                for start, end in [(left[0], right[0]), (right[0], left[0])]:
                    if (
                        timeaxis is None
                        and len(self._linkreac[start]) < linknum
                        and end not in self._linkreac[start]
                    ):
                        self._linkreac[start].append(end)
        return reaction

    def _readreactionabcd(self):
        reactionsabcd = []
        if os.path.isfile(self.reactionabcdfilename):
            with open(self.reactionabcdfilename) as f:
                for i, line in enumerate(f, 1):
                    left, right, num = self._handlereaction(line)
                    reactionsabcd.append({"i": i, "l": left, "r": right, "n": num})
                    for spec in left + right:
                        self._svgspecs.add(spec)
        return reactionsabcd

    def _readspecies(self, reaction, timeaxis=None):
        specs = []
        for reac in reaction:
            for spec in (reac["l"][0], reac["r"][0]):
                if spec not in specs:
                    specs.append(spec)
                    if timeaxis is None:
                        self._svgspecs.add(spec)
        # return list of dict
        return [{"s": spec, "i": i} for i, spec in enumerate(specs, 1)]

    def _readdata(self):
        self._reaction = [self._readreaction()]
        self._specs = [self._readspecies(self._reaction[0])]
        self._reactionsabcd = self._readreactionabcd()
        if self.split > 1:
            for i in range(self.split):
                reaction = self._readreaction(timeaxis=i)
                self._reaction.append(reaction)
                self._specs.append(self._readspecies(reaction, timeaxis=i))

    def _generateresult(self):
        network = [self._generatenetwork()]
        if self.split > 1:
            for i in range(self.split):
                network.append(self._generatenetwork(timeaxis=i))
        self._templatedict["network"] = network
        self._templatedict["species"] = self._specs
        self._templatedict["reactions"] = self._reaction
        self._templatedict["reactionsabcd"] = self._reactionsabcd
        self._templatedict["linkreac"] = self._linkreac
        rngdata = json.dumps(self._templatedict, separators=(",", ":"))
        with open(
            os.path.join(os.path.dirname(__file__), "static", "webpack", "bundle.html")
        ) as f:
            template = f.read()
        webpage = template.replace("PUTREACNETGENERATORDATAHERE", rngdata)
        with open(self.jsonfilename, "w") as f:
            f.write(rngdata)
        with open(self.resultfilename, "w", encoding="utf-8") as f:
            f.write(webpage)

    def _generatenetwork(self, timeaxis=None):
        with open(
            self.imagefilename
            if timeaxis is None
            else f"{self.imagefilename}.{timeaxis}"
        ) as f:
            svgdata = f.read().strip()
            svgdata = re.sub(r'width="\d+(\.\d+)?pt"', 'width="100%"', svgdata, count=2)
            svgdata = re.sub(r'height="\d+(\.\d+)?pt"', "", svgdata, count=2)
            svgdata = re.sub(
                r"""<(\?xml|\!DOCTYPE|\!\-\-)("[^"]*"|'[^']*'|[^'">])*>""", "", svgdata
            )
            svgdata = svgdata.replace(
                r"""<style type="text/css">*{""",
                r"""<style type="text/css">#network svg *{""",
            )
        return scour.scour.scourString(svgdata, SCOUROPTIONS)
