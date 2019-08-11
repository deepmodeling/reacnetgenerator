# -*- coding: utf-8 -*-
# cython: language_level=3
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
import os
from collections import defaultdict
from multiprocessing import Pool
import pkg_resources

import htmlmin
import openbabel
import scour.scour
from jinja2 import Template
from tqdm import tqdm


class _HTMLResult:
    def __init__(self, rng):
        self._reactionfile = rng.reactionfilename
        self._resultfile = rng.resultfilename
        self._imagefile = rng.imagefilename
        self._reactionabcdfilename = rng.reactionabcdfilename
        self._nproc = rng.nproc
        self._split = rng.split
        self._atomname = rng.atomname
        self._templatedict = {
            "speciesshownum": 30,
            "reactionsshownum": 20,
        }
        self._linkreac = defaultdict(list)
        # define instance
        self._specs = None
        self._reaction = None
        self._reactionsabcd = None
        self._svgfiles = {}
        self.scouroptions = rng.SCOUROPTIONS

    def report(self):
        """Generate a web page to show the result."""
        self._readdata()
        self._generateresult()
        logging.info(
            f"Report is generated. Please see {self._resultfile} for more details.")

    def _re(self, smi):
        for an in self._atomname:
            if an != 'H':
                smi = smi.replace(an.upper(), f"[{an.upper()}]").replace(an.lower(), f"[{an.lower()}]")
        return smi.replace("[HH]", "[H]")

    def _handlereaction(self, line):
        sx = line.split()
        left, right = sx[1].split("->")
        left = list([self._re(spec) for spec in left.split("+")])
        right = list([self._re(spec) for spec in right.split("+")])
        num = int(sx[0])
        return left, right, num

    def _readreaction(self, timeaxis=None, linknum=6):
        reaction = []
        with open(self._reactionfile if timeaxis is None else f"{self._reactionfile}.{timeaxis}") as f:
            for i, line in enumerate(f, 1):
                left, right, num = self._handlereaction(line)
                reaction.append({"i":i, "l":left, "r":right, "n":num})
                for start, end in [(left[0], right[0]), (right[0], left[0])]:
                    if timeaxis is None and len(self._linkreac[start]) < linknum:
                        self._linkreac[start].append(end)
        return reaction

    def _readreactionabcd(self):
        reactionsabcd = []
        if os.path.isfile(self._reactionabcdfilename):
            append_spec = set()
            with open(self._reactionabcdfilename) as f:
                for i, line in enumerate(f, 1):
                    left, right, num = self._handlereaction(line)
                    reactionsabcd.append({"i":i, "l": left, "r": right, "n": num})
                    for spec in left + right:
                        if spec not in self._svgfiles:
                            append_spec.add(spec)
            if append_spec:
                with Pool(self._nproc) as pool:
                    results = pool.imap_unordered(self._convertsvg, tqdm(append_spec))
                    for spec, svgfile in results:
                        self._svgfiles[spec] = svgfile
                pool.join()
                pool.close()
        return reactionsabcd

    def _convertsvg(self, smiles):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("smi", "svg")
        obConversion.AddOption('x')
        mol = openbabel.OBMol()
        obConversion.ReadString(mol, smiles)
        svgdata = obConversion.WriteString(mol)
        svgdata = scour.scour.scourString(svgdata, self.scouroptions)
        svgdata = re.sub(r"\d+(\.\d+)?px", "100%", svgdata, count=2)
        svgdata = re.sub(
            r"""<rect("[^"]*"|'[^']*'|[^'">])*>""", '', svgdata)
        svgdata = re.sub(
            r"""<\?xml("[^"]*"|'[^']*'|[^'">])*>""", '', svgdata)
        return smiles, svgdata

    def _readspecies(self, reaction, timeaxis=None):
        specs = []
        for reac in reaction:
            for spec in (reac['l'][0], reac['r'][0]):
                if spec not in specs:
                    specs.append(spec)
        if timeaxis is None:
            with Pool(self._nproc) as pool:
                results = pool.imap_unordered(self._convertsvg, tqdm(specs))
                for spec, svgfile in results:
                    self._svgfiles[spec] = svgfile
            pool.join()
            pool.close()
        # return list of dict
        return list([{"s": spec, "i":i} for i, spec in enumerate(specs, 1)])

    def _readdata(self):
        self._reaction = [self._readreaction()]
        self._specs = [self._readspecies(self._reaction[0])]
        self._reactionsabcd = self._readreactionabcd()
        if self._split > 1:
            for i in range(self._split):
                reaction = self._readreaction(timeaxis=i)
                self._reaction.append(reaction)
                self._specs.append(self._readspecies(reaction, timeaxis=i))

    def _generateresult(self):
        network = [self._generatenetwork()]
        if self._split > 1:
            for i in range(self._split):
                network.append(self._generatenetwork(timeaxis=i))
        self._templatedict["network"] = json.dumps(network, separators=(',', ':'))
        self._templatedict["speciessvg"] = json.dumps(self._svgfiles, separators=(',', ':'))
        self._templatedict["species"] = json.dumps(self._specs, separators=(',', ':'))
        self._templatedict["reactions"] = json.dumps(self._reaction, separators=(',', ':'))
        self._templatedict["reactionsabcd"] = json.dumps(self._reactionsabcd, separators=(',', ':'))
        self._templatedict["linkreac"] = json.dumps(self._linkreac, separators=(',', ':'))
        template = Template(pkg_resources.resource_string(
            __name__, 'static/webpack/bundle.html').decode())
        webpage = template.render(**self._templatedict)
        with open(self._resultfile, 'w', encoding="utf-8") as f:
            f.write(htmlmin.minify(webpage))

    def _generatenetwork(self, timeaxis=None):
        with open(self._imagefile if timeaxis is None else f"{self._imagefile}.{timeaxis}") as f:
            svgdata = f.read().strip()
            svgdata = re.sub(r'width="\d+(\.\d+)?pt"', 'width="100%"', svgdata, count=2)
            svgdata = re.sub(r'height="\d+(\.\d+)?pt"', '', svgdata, count=2)
            svgdata = re.sub(
                r"""<(\?xml|\!DOCTYPE|\!\-\-)("[^"]*"|'[^']*'|[^'">])*>""", '',
                svgdata)
            svgdata = svgdata.replace(r"""<style type="text/css">*{""",r"""<style type="text/css">#network svg *{""")
        return htmlmin.minify(svgdata)
