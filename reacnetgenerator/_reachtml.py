# -*- coding: utf-8 -*-
# Generate a web page to show the analysis report

import re
from collections.abc import Mapping
from multiprocessing import Pool
import htmlmin
import openbabel
from ._htmlstatic import _static_js, _static_css, _static_img
from ._htmltemplate import _html


class _HTMLResult(object):
    def __init__(self, ReacNetGenerator):
        self._reactionfile = ReacNetGenerator.reactionfilename
        self._resultfile = ReacNetGenerator.resultfilename
        self._imagefile = ReacNetGenerator.imagefilename
        self._nproc = ReacNetGenerator.nproc
        self._script = []
        self._result = []
        self._svgfiles = {}

    def _report(self):
        self._readdata()
        self._generateresult()

    def _re(self, smi):
        return smi.replace("O", "[O]").replace("C", "[C]").replace("[HH]", "[H]")

    def _readreaction(self):
        reaction = []
        with open(self._reactionfile) as f:
            for line in f:
                sx = line.split()
                s = sx[1].split("->")
                reaction.append((self._re(s[0]), self._re(s[1]), int(sx[0])))
        return reaction

    def _convertsvg(self, smiles):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("smi", "svg")
        mol = openbabel.OBMol()
        obConversion.ReadString(mol, smiles)
        svgfile = obConversion.WriteString(mol)
        return smiles, svgfile

    def _readspecies(self):
        specs = []
        for reac in self._reaction:
            for spec in reac[:2]:
                if not spec in specs:
                    specs.append(spec)
        with Pool(self._nproc) as pool:
            results = pool.imap(self._convertsvg, specs)
            for spec, svgfile in results:
                self._svgfiles[spec] = svgfile
        return specs

    def _readdata(self):
        self._reaction = self._readreaction()
        self._specs = self._readspecies()

    def _generateresult(self):
        self._result.append(_html['page-top'].format("".join([
            _static_css["bootstrap.min.css"],
            _static_css["creative.min.css"],
            _static_css["magnific-popup.min.css"],
            _html['bk-css'] % (_static_img["fire.jpg"])])))
        self._generatenetwork()
        self._generatesvg()
        self._generatespecies()
        self._generatereaction()
        self._result.append(_html['page-bottom'].format("".join([
            _static_js["jquery.min.js"],
            _static_js["bootstrap.bundle.min.js"],
            _static_js["jquery.easing.min.js"],
            _static_js["scrollreveal.min.js"],
            _static_js["jquery.magnific-popup.min.js"],
            _static_js["creative.min.js"],
            "".join(self._script)])))
        with open(self._resultfile, 'w', encoding="utf-8") as f:
            print(htmlmin.minify("".join(self._result)), file=f)

    def _generatenetwork(self):
        with open(self._imagefile) as f:
            svgdata = f.read().strip()
            svgdata = re.sub(r"\d+(\.\d+)?pt", "100%", svgdata, count=2)
        self._result.append(_html['network'] % svgdata)

    def _generatesvg(self):
        buff = [_html['speciessvg-top']]
        for spec in self._specs:
            svgdata = self._svgfiles[spec]
            svgdata = re.sub(r"\d+(\.\d+)?px", "100%", svgdata, count=2)
            buff.append(_html['speciessvg-each'] % (spec, svgdata))
        buff.append(_html['speciessvg-bottom'])
        self._result.extend(buff)

    def _generatespecies(self, line=10, shownum=30):
        buff = [_html['species-top']]
        for i, spec in enumerate(self._specs):
            buff.append(_html['species-each'] % (((_html['tr-top'] if i < shownum else _html['tr-specnone'])
                                                  if i % line == 0 else ""), spec, str(i+1), (_html['tr-bottom'] if i % line == line-1 else "")))
        buff.append(_html['species-bottom'])
        if len(self._specs) <= shownum:
            self._script .append(_html['script-hidespec'])
        self._result.extend(buff)

    def _generatereaction(self, line=4, reacnum=True, shownum=20):
        buff = [_html['reactions-top']]
        for i, reac in enumerate(self._reaction):
            buff.append(_html['reactions-each'] % (((_html['tr-top'] if i < shownum else _html['tr-reacnone']) if i % line == 0 else ""), str(
                i+1), reac[0], (str(reac[2]) if reacnum else ""), _html['narrowurl'], reac[1], (_html['tr-bottom'] if i % line == line-1 else "")))
        buff.append(_html['reactions-bottom'])
        if len(self._reaction) <= shownum:
            self._script.append(_html['script-hidereac'])
        self._result.extend(buff)
