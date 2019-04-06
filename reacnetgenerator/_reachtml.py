# -*- coding: utf-8 -*-
"""Generate a web page to show the analysis report.

For the convenience analyze of the trajectories, we put all the results
generated by ReacNetGen (Network, Species, and Reactions) in an interactive
web page. By default, 20 species with the most reactions are taken to draw the
network. However, by clicking on a given species, one can check the special
network which starts from it.
"""

import json
import logging
import re
from collections import defaultdict
from multiprocessing import Pool
import pkg_resources

import htmlmin
import openbabel
import scour.scour
from jinja2 import Template


class _HTMLResult:
    def __init__(self, rng):
        self._reactionfile = rng.reactionfilename
        self._resultfile = rng.resultfilename
        self._imagefile = rng.imagefilename
        self._nproc = rng.nproc
        self._split = rng.split
        self._templatedict = {
            "speciesshownum": 30,
            "reactionsshownum": 20,
        }
        self._linkreac = defaultdict(list)
        # define instance
        self._specs = None
        self._reaction = None
        self._svgfiles = {}

    def report(self):
        """Generate a web page to show the result."""
        self._readdata()
        self._generateresult()
        logging.info(
            f"Report is generated. Please see {self._resultfile} for more details.")

    @classmethod
    def _re(cls, smi):
        return smi.replace(
            "O", "[O]").replace(
            "C", "[C]").replace(
            "[HH]", "[H]")

    def _readreaction(self, timeaxis=None):
        reaction = []
        with open(self._reactionfile if timeaxis is None else f"{self._reactionfile}.{timeaxis}") as f:
            for line in f:
                sx = line.split()
                s = sx[1].split("->")
                left, right, num = self._re(s[0]), self._re(s[1]), int(sx[0])
                reaction.append((left, right, num))
                if timeaxis is None and len(self._linkreac[left]) < 5:
                    self._linkreac[left].append(right)
        return reaction

    @classmethod
    def _convertsvg(cls, smiles):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("smi", "svg")
        obConversion.AddOption('x')
        mol = openbabel.OBMol()
        obConversion.ReadString(mol, smiles)
        svgdata = obConversion.WriteString(mol)
        svgdata = scour.scour.scourString(svgdata)
        svgdata = re.sub(r"\d+(\.\d+)?px", "100%", svgdata, count=2)
        svgdata = re.sub(
            r"""<rect("[^"]*"|'[^']*'|[^'">])*>""", '', svgdata)
        svgdata = re.sub(
            r"""<\?xml("[^"]*"|'[^']*'|[^'">])*>""", '', svgdata)
        svgdata = re.sub(r"""<title>.*?<\/title>""", '', svgdata)
        return smiles, svgdata

    def _readspecies(self, reaction, timeaxis=None):
        specs = []
        for reac in reaction:
            for spec in reac[:2]:
                if spec not in specs:
                    specs.append(spec)
        if timeaxis is None:
            with Pool(self._nproc) as pool:
                results = pool.imap_unordered(self._convertsvg, specs)
                for spec, svgfile in results:
                    self._svgfiles[spec] = svgfile
            pool.join()
            pool.close()
        return specs

    def _readdata(self):
        self._reaction = [self._readreaction()]
        self._specs = [self._readspecies(self._reaction[0])]
        if self._split > 1:
            for i in range(self._split):
                reaction = self._readreaction(timeaxis=i)
                self._reaction.append(reaction)
                self._specs.append(self._readspecies(reaction, timeaxis=i))

    def _generateresult(self):
        self._templatedict["network_time"] = [self._generatenetwork()]
        if self._split > 1:
            for i in range(self._split):
                self._templatedict["network_time"].append(self._generatenetwork(timeaxis=i))
        self._generatesvg()
        self._templatedict["speciestime"] = self._specs
        self._templatedict["reactionstime"] = self._reaction
        self._templatedict["javascript"] = pkg_resources.resource_string(
            __name__, 'static/webpack/bundle.js').decode()
        self._templatedict["linkreac"] = json.dumps(
            self._linkreac, separators=(',', ':'))
        template = Template(pkg_resources.resource_string(
            __name__, 'static/template.html').decode())
        webpage = template.render(**self._templatedict)
        with open(self._resultfile, 'w', encoding="utf-8") as f:
            f.write(htmlmin.minify(webpage))

    def _generatenetwork(self, timeaxis=None):
        with open(self._imagefile if timeaxis is None else f"{self._imagefile}.{timeaxis}") as f:
            svgdata = f.read().strip()
            svgdata = re.sub(r"\d+(\.\d+)?pt", "100%", svgdata, count=2)
            svgdata = re.sub(
                r"""<(\?xml|\!DOCTYPE|\!\-\-)("[^"]*"|'[^']*'|[^'">])*>""", '',
                svgdata)
        return svgdata

    def _generatesvg(self):
        self._templatedict["speciessvg"] = list(
            [{"name": spec, "svg": self._svgfiles[spec]}
             for spec in self._specs[0]])
