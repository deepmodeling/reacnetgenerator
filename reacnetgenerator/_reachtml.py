# -*- coding: utf-8 -*-
# Generate a web page to show the analysis report

import re
from multiprocessing import Pool
import htmlmin
import openbabel
from ._htmlstatic import _static_js, _static_css, _static_img


class _HTMLResult(object):
    def __init__(self, ReacNetGenerator):
        self._reactionfile = ReacNetGenerator.reactionfilename
        self._resultfile = ReacNetGenerator.resultfilename
        self._imagefile = ReacNetGenerator.imagefilename
        self._nproc = ReacNetGenerator.nproc
        self._script = ""
        self._svgfiles = {}

    def _report(self):
        self._readdata()
        self._generateresult()

    def _re(self, smi):
        return smi.replace("O", "[O]").replace("C", "[C]").replace("[HH]", "[H]")

    def _readreaction(self):
        reaction = []
        with open(_self.reactionfile) as f:
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
        svgfile=obConversion.WriteString(mol)
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
        self._result = self._html['page-top']
        self._generatenetwork()
        self._generatesvg()
        self._generatespecies()
        self._generatereaction()
        self._result += self._html['page-bottom']
        with open(self._resultfile, 'w', encoding="utf-8") as f:
            print(htmlmin.minify(self._result), file=f)

    def _generatenetwork(self):
        with open(self._imagefile) as f:
            svgdata = f.read().strip()
            svgdata = re.sub(r"\d+(\.\d+)?pt", "100%", svgdata, count=2)
        self._result += self._html['network'] % svgdata

    def _generatesvg(self):
        buff = self._html['speciessvg-top']
        for spec in self._specs:
            svgdata = self._svgfiles[spec]
            svgdata = re.sub(r"\d+(\.\d+)?px", "100%", svgdata, count=2)
            buff += self._html['speciessvg-each'] % (spec, svgdata)
        buff += self._html['speciessvg-bottom']
        self._result += buff

    def _generatespecies(self, line=10, shownum=30):
        buff = self._html['species-top']
        for i, spec in enumerate(self._specs):
            buff += self._html['species-each'] % ((("<tr>" if i < shownum else "<tr class='specnone'>") if i %
                                                  line == 0 else ""), spec, str(i+1), ("</tr>" if i % line == line-1 else ""))
        buff += self._html['species-bottom']
        if len(self._specs) <= shownum:
            self._script += self._html['script-hidespec']
        self._result += buff

    def _generatereaction(self, line=4, reacnum=True, shownum=20):
        buff = self._html['reactions-top']
        for i, reac in enumerate(self._reaction):
            buff += self._html['reactions-each'] % ((("<tr>" if i < shownum else "<tr class='reacnone'>") if i % line == 0 else ""), str(
                i+1), reac[0], (str(reac[2]) if reacnum else ""), self._html['narrowurl'], reac[1], ("</tr>" if i % line == line-1 else ""))
        buff += self._html['reactions-bottom']
        if len(self._reaction) <= shownum:
            self._script += self._html['script-hidereac']
        self._result += buff

    @property
    def _html(self):
        return {
            'page-top': """
                <html><head>
                <meta charset="utf-8">
                <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
                <title>Analysis Report</title>
                </head><body id="page-top">
                <style>{}</style>
                <svg class="d-none">
                    <defs>
                        <path id="narrow" style="text-align:start;line-height:100%;-inkscape-font-specification:Arial Unicode MS" d="M24.35 7.613c-3.035 1.11-5.407 2.908-7.113 5.395h-1.299c.585-1.743 1.567-3.39 2.945-4.938H.649V6.278h18.234c-1.378-1.548-2.36-3.2-2.945-4.956h1.299c1.706 2.487 4.078 4.285 7.114 5.395v.896z" font-size="39.506" font-weight="400" font-family="Arial Unicode MS"/>
                    </defs>
                </svg>
                <nav class="navbar navbar-expand-lg navbar-dark fixed-top" id="mainNav">
                    <div class="container">
                        <a class="navbar-brand js-scroll-trigger" href="#page-top">ReacNetGenerator</a>
                        <button class="navbar-toggler navbar-toggler-right" type="button" data-toggle="collapse" data-target="#navbarResponsive" aria-controls="navbarResponsive" aria-expanded="false" aria-label="Toggle navigation">
                            <span class="navbar-toggler-icon"></span>
                        </button>
                        <div class="collapse navbar-collapse" id="navbarResponsive">
                            <ul class="navbar-nav ml-auto">
                                <li class="nav-item"><a class="nav-link js-scroll-trigger" href="#network">Network</a></li>
                                <li class="nav-item"><a class="nav-link js-scroll-trigger" href="#species">Species</a></li>
                                <li class="nav-item"><a class="nav-link js-scroll-trigger" href="#reactions">Reactions</a></li>
                            </ul>
                        </div>
                    </div>
                </nav>
                <header class="masthead text-center text-white d-flex bg-dark imgbg" id="info">
                    <div class="container my-auto">
                        <div class="row">
                            <div class="col-lg-10 mx-auto">
                                <h1 class="text-uppercase">Analysis Report</h1>
                            </div>
                            <div class="col-lg-10 mx-auto">
                                <p><strong>Please cite:</strong> J. Zeng, L. Cao, J.Z.H. Zhang, C.-H. Chen, T. Zhu: ReacNetGen: an Automatic Reaction Network Generator for Reactive Molecular Dynamic Simulations, to be submitted</p>
                            </div>
                        </div>
                        <div class="row">
                            <div class="col-lg-3 mx-auto text-center">
                                <div class="service-box mt-5 mx-auto">
                                    <a class="btn btn-primary btn-xl js-scroll-trigger" href="#network">Network</a>
                                </div>
                            </div>
                            <div class="col-lg-3 mx-auto text-center">
                                <div class="service-box mt-5 mx-auto">
                                    <a class="btn btn-primary btn-xl js-scroll-trigger" href="#species">Species</a>
                                </div>
                            </div>
                            <div class="col-lg-3 mx-auto text-center">
                                <div class="service-box mt-5 mx-auto">
                                    <a class="btn btn-primary btn-xl js-scroll-trigger" href="#reactions">Reactions</a>
                                </div>
                            </div>
                        </div>
                    </div>
                </header>
                """.format("".join([
                _static_css["bootstrap.min.css"],
                _static_css["creative.min.css"],
                _static_css["magnific-popup.min.css"],
                """.spec{height:100px;width:100px;}.content-table{font-size:14;margin-left:auto;margin-right:auto;text-align:center;}.reacnum{color:blue;}.specnone,.reacnone{display:none;}#info{background-image:url(%s)}""" % (
                    _static_img["fire.jpg"])
            ])),
            'page-bottom': """
                <section class="bg-dark text-white" id="foot">
                  <div class="container">
                    <div class="row">
                      <div class="col-lg-8 mx-auto text-center">
                        <p>Generated by <a href="https://njzjz.github.io/reacnetgenerator">ReacNetGenerator</a></p>
                        <p><strong>Citation:</strong> J. Zeng, L. Cao, J.Z.H. Zhang, C.-H. Chin, T. Zhu: ReacNetGen: an Automatic Reaction Network Generator for Reactive Molecular Dynamic Simulations, to be submitted</p>
                        <p><strong>Author:</strong> <a href="https://cv.njzjz.win/">Jinzhe Zeng</a>, Liqun Cao, <a href="https://research.shanghai.nyu.edu/centers-and-institutes/chemistry/people/john-zenghui-zhang">John ZH Zhang</a>, Chih-Hao Chin, <a href="http://computchem.cn/people/">Tong Zhu</a></p>
                        <p><strong>Email:</strong> <a href="mailto:tzhu@lps.ecnu.edu.cn">tzhu@lps.ecnu.edu.cn</a>, <a href="mailto:jzzeng@stu.ecnu.edu.cn">jzzeng@stu.ecnu.edu.cn</a></p>
                      </div>
                    </div>
                  </div>
                </section>
                <script>{}</script>
                </body></html>
                """.format("".join([
                _static_js["jquery.min.js"],
                _static_js["bootstrap.bundle.min.js"],
                _static_js["jquery.easing.min.js"],
                _static_js["scrollreveal.min.js"],
                _static_js["jquery.magnific-popup.min.js"],
                _static_js["creative.min.js"],
                self._script
            ])),
            'network': """<section id="network" class='bg-white mx-auto text-center'>
                    <div class="container my-auto">
                        <div class="row">
                            <div class="mx-auto w-100 text-center">
                                <h2>Reaction Network</h2>
                                %s
                            </div>
                        </div>
                    </div>
                </section>
                """,
            'species-top': """<section id="species" class='bg-white mx-auto'>
                <div class="container">
                    <div class="row">
                        <div class="mx-auto text-center">
                            <h2 class='text-center'>Species</h2>
                            <table class='content-table'>
                """,
            'species-each': """%s<td><svg class='spec'><use xlink:href="#%s"></use></svg><br/>%s</td>%s
                """,
            'species-bottom': """</table>
                        </div>
                        <div class="mx-auto text-center mt-5">
                            <a class="btn btn-primary btn-xl js-scroll-trigger" href="javascript:$('.specnone').show();$('#showmorespec').hide();" id="showmorespec">Show all</a>
                        </div>
                    </div>
                </section>
                """,
            'reactions-top': """<section id="reactions" class='bg-white mx-auto'>
                <div class="container">
                    <div class="row">
                        <div class="mx-auto text-center">
                            <h2 class='text-center'>Reactions</h2>
                            <h4 class='text-center'>(sorted by frequency)</h4>
                            <table class='content-table'>
                """,
            'reactions-each': """%s<td>%s</td><td><svg class='spec'><use xlink:href="#%s"></use></svg></td>
                <td class="reacnum">%s<br/>%s</td>
                <td><svg class='spec'><use xlink:href="#%s"></use></svg></td>
                %s
                """,
            'reactions-bottom': """</table>
                        </div>
                        <div class="mx-auto text-center mt-5">
                            <a class="btn btn-primary btn-xl js-scroll-trigger" href="javascript:$('.reacnone').show();$('#showmorereac').hide();" id="showmorereac">Show all</a>
                        </div>
                    </div>
                </section>
                """,
            'narrowurl': '''<svg width="25" height="14.33" version="1"><use xlink:href="#narrow"/></svg>''',
            'speciessvg-top': '''<svg class="d-none"><defs>
            ''',
            'speciessvg-each': '''<svg id="%s">%s</svg>
            ''',
            'speciessvg-bottom': '''</defs></svg>
            ''',
            'script-hidereac': '''$('#showmorereac').hide();''',
            'script-hidespec': '''$('#showmorespec').hide();''',
        }
