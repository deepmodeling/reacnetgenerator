# -*- coding: utf-8 -*-
import os
from multiprocessing.pool import ThreadPool
from multiprocessing import cpu_count
import htmlmin
import re
from ReacNetGenerator.htmlstatic import static_js,static_css,static_img

class HTMLResult(object):
    def __init__(self,reactionfile,resultfile,imagefile,specfolder="species",n_thread=None):
        self.reactionfile=reactionfile
        self.resultfile=resultfile
        self.specfolder=specfolder
        self.imagefile=imagefile
        self.n_thread=n_thread if n_thread else cpu_count()

    def re(self,smi):
        return smi.replace("O","[O]").replace("C","[C]").replace("[HH]","[H]")

    def readreaction(self):
        reaction=[]
        with open(self.reactionfile) as f:
            for line in f:
                sx=line.split()
                s=sx[1].split("->")
                reaction.append((self.re(s[0]),self.re(s[1]),int(sx[0])))
        return reaction

    def svgfilename(self,smiles,url=False):
        if url:
            smiles=smiles.replace("#","%23")
        return os.path.join(self.specfolder,smiles+'.svg')

    def convertsvg(self,smiles,folder="species"):
        filename=self.svgfilename(smiles)
        command=False
        if not os.path.exists(filename):
            if not os.path.exists(folder):
                os.makedirs(folder)
            command="obabel -:'%s' -osvg -O '%s'"%(smiles,filename)
        return command

    def readspecies(self):
        specs=[]
        commands=[]
        for reac in self.reaction:
            for spec in reac[:2]:
                if not spec in specs:
                    specs.append(spec)
                    command=self.convertsvg(spec)
                    if command:
                        commands.append(command)
        with ThreadPool(self.n_thread) as pool:
            results=pool.imap(os.system,commands)
            for result in results:
                pass
        return specs

    def report(self):
        self.readdata()
        self.generateresult()

    def readdata(self):
        self.reaction=self.readreaction()
        self.specs=self.readspecies()
 
    def generateresult(self):
        self.result=self.html['page-top']
        self.generatenetwork()
        self.generatesvg()
        self.generatespecies()
        self.generatereaction()
        self.result+=self.html['page-bottom']
        with open(self.resultfile,'w') as f:
            print(htmlmin.minify(self.result),file=f)

    def generatenetwork(self):
        with open(self.imagefile) as f:
            svgdata=f.read().strip()
            svgdata=re.sub(r"\d+(\.\d+)?pt","100%",svgdata,count=2)
        self.result+=self.html['network']%svgdata
    
    def generatesvg(self):
        buff=self.html['speciessvg-top']
        for spec in self.specs:
            with open(self.svgfilename(spec)) as f:
                svgdata=f.read().strip()
                svgdata=re.sub(r"\d+(\.\d+)?px","100%",svgdata,count=2)
            buff+=self.html['speciessvg-each']%(spec,svgdata)
        buff+=self.html['speciessvg-bottom']
        self.result+=buff
    
    def generatespecies(self,line=10,shownum=30):
        buff=self.html['species-top']
        for i,spec in enumerate(self.specs):
            buff+=self.html['species-each']%((("<tr>" if i<shownum else "<tr class='specnone'>") if i%line==0 else ""),spec,str(i+1),("</tr>" if i%line==line-1 else ""))
        buff+=self.html['species-bottom']
        self.result+=buff

    def generatereaction(self,line=4,reacnum=True,shownum=20):
        buff=self.html['reactions-top']
        for i,reac in enumerate(self.reaction):
            buff+=self.html['reactions-each']%((("<tr>" if i<shownum else "<tr class='reacnone'>") if i%line==0 else ""),str(i+1),reac[0],(str(reac[2]) if reacnum else ""),self.html['narrowurl'],reac[1],("</tr>" if i%line==line-1 else ""))
        buff+=self.html['reactions-bottom']
        self.result+=buff

    @property
    def html(self):
        return {
            'page-top':"""
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
                                <p><strong>Please cite:</strong> J. Zeng, L. Cao, J.Z.H. Zhang, T Zhu: ReacNetGen: an Automatic Reaction Network Generator for Reactive Molecular Dynamic Simulations, to be submitted</p>
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
                    static_css["bootstrap.min.css"],
                    static_css["creative.min.css"],
                    static_css["magnific-popup.min.css"],
                    """.spec{height:100px;width:100px;}.content-table{font-size:14;margin-left:auto;margin-right:auto;text-align:center;}.reacnum{color:blue;}.specnone,.reacnone{display:none;}#info{background-image:url(%s)}"""%(static_img["fire.jpg"])
                ])),
            'page-bottom':"""
            	<section class="bg-dark text-white" id="foot">
                  <div class="container">
                    <div class="row">
                      <div class="col-lg-8 mx-auto text-center">
                        <p>Generated by <a href="https://njzjz.github.io/ReacNetGenerator">ReacNetGenerator</a></p>
                        <p><strong>Citation:</strong> J. Zeng, L. Cao, J.Z.H. Zhang, T. Zhu: ReacNetGen: an Automatic Reaction Network Generator for Reactive Molecular Dynamic Simulations, to be submitted</p>
                        <p><strong>Author:</strong> <a href="https://cv.njzjz.win/">Jinzhe Zeng</a>, Liqun Cao, <a href="https://research.shanghai.nyu.edu/centers-and-institutes/chemistry/people/john-zenghui-zhang">John ZH Zhang</a>, <a href="http://computchem.cn/people/">Tong Zhu</a></p>
                        <p><strong>Email:</strong> <a href="mailto:tzhu@lps.ecnu.edu.cn">tzhu@lps.ecnu.edu.cn</a>, <a href="mailto:jzzeng@stu.ecnu.edu.cn">jzzeng@stu.ecnu.edu.cn</a></p>
                      </div>
                    </div>
                  </div>
                </section>
                <script>{}</script>
                </body></html>
                """.format("".join([
                    static_js["jquery.min.js"],
                    static_js["bootstrap.bundle.min.js"],
                    static_js["jquery.easing.min.js"],
                    static_js["scrollreveal.min.js"],
                    static_js["jquery.magnific-popup.min.js"],
                    static_js["creative.min.js"],
                ])),
            'network':"""<section id="network" class='bg-white mx-auto text-center'>
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
            'species-top':"""<section id="species" class='bg-white mx-auto'>
                <div class="container">
                    <div class="row">
                        <div class="mx-auto text-center">
                            <h2 class='text-center'>Species</h2>
                            <table class='content-table'>
                """,
            'species-each':"""%s<td><svg class='spec'><use xlink:href="#%s"></use></svg><br/>%s</td>%s
                """,
            'species-bottom':"""</table>
                        </div>
                        <div class="mx-auto text-center mt-5">
                            <a class="btn btn-primary btn-xl js-scroll-trigger" href="javascript:$('.specnone').show();$('#showmorespec').hide();" id="showmorespec">Show all</a>
                        </div>
                    </div>
                </section>
                """,
            'reactions-top':"""<section id="reactions" class='bg-white mx-auto'>
                <div class="container">
                    <div class="row">
                        <div class="mx-auto text-center">
                            <h2 class='text-center'>Reactions</h2>
                            <h4 class='text-center'>(sorted by frequency)</h4>
                            <table class='content-table'>
                """,
            'reactions-each':"""%s<td>%s</td><td><svg class='spec'><use xlink:href="#%s"></use></svg></td>
                <td class="reacnum">%s<br/>%s</td>
                <td><svg class='spec'><use xlink:href="#%s"></use></svg></td>
                %s
                """,
            'reactions-bottom':"""</table>
                        </div>
                        <div class="mx-auto text-center mt-5">
                            <a class="btn btn-primary btn-xl js-scroll-trigger" href="javascript:$('.reacnone').show();$('#showmorereac').hide();" id="showmorereac">Show all</a>
                        </div>
                    </div>
                </section>
                """,
            'narrowurl':'''<svg width="25" height="14.33" version="1"><use xlink:href="#narrow"/></svg>''',
            'speciessvg-top':'''<svg class="d-none"><defs>
            ''',
            'speciessvg-each':'''<svg id="%s">%s</svg>
            ''',
            'speciessvg-bottom':'''</defs></svg>
            '''
        }
