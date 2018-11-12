import os
from multiprocessing.pool import ThreadPool
from multiprocessing import cpu_count

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
        self.generatespecies()
        self.generatereaction()
        self.result+=self.html['page-bottom']
        with open(self.resultfile,'w') as f:
            print(self.result,file=f)

    def generatenetwork(self):
        self.result+=self.html['network']%self.imagefile

    def generatespecies(self,line=10,shownum=30):
        buff=self.html['species-top']
        for i,spec in enumerate(self.specs):
            buff+=self.html['species-each']%((("<tr>" if i<shownum else "<tr class='specnone'>") if i%line==0 else ""),self.svgfilename(spec,url=True),str(i+1),("</tr>" if i%line==line-1 else ""))
        buff+=self.html['species-bottom']
        self.result+=buff

    def generatereaction(self,line=4,reacnum=True,shownum=20):
        buff=self.html['reactions-top']
        for i,reac in enumerate(self.reaction):
            buff+=self.html['reactions-each']%((("<tr>" if i<shownum else "<tr class='reacnone'>") if i%line==0 else ""),str(i+1),self.svgfilename(reac[0],url=True),(str(reac[2]) if reacnum else ""),self.html['narrowurl'],self.svgfilename(reac[1],url=True),("</tr>" if i%line==line-1 else ""))
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
                <link href="https://lib.baomitu.com/twitter-bootstrap/4.1.3/css/bootstrap.min.css" rel="stylesheet">
                <link href="https://lib.baomitu.com/font-awesome/4.7.0/css/font-awesome.min.css" rel="stylesheet" type="text/css">
                <link href='https://fonts.googleapis.com/css?family=Open+Sans:300italic,400italic,600italic,700italic,800italic,400,300,600,700,800' rel='stylesheet' type='text/css'>
                <link href='https://fonts.googleapis.com/css?family=Merriweather:400,300,300italic,400italic,700,700italic,900,900italic' rel='stylesheet' type='text/css'>
                <link href="https://lib.baomitu.com/startbootstrap-creative/5.0.1/css/creative.min.css" rel="stylesheet">
                <link href="https://lib.baomitu.com/magnific-popup.js/1.1.0/magnific-popup.min.css" rel="stylesheet">
                </head><body id="page-top">
                <style>
                    .spec{height:100px;width:100px;}
                    .content-table{font-size:14;margin-left:auto;margin-right:auto;text-align:center;}
                    .reacnum{color:blue;}
                    .specnone,.reacnone{display:none;}
                    #info{background-image:url(https://ww1.sinaimg.cn/large/005YhI8igy1fx5kf3u1v9j30xc0kmn12);}
                </style>
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
                <img src="https://ww1.sinaimg.cn/large/005YhI8igy1fx5kf3u1v9j30xc0kmn12" class="d-none">
                <header class="masthead text-center text-white d-flex bg-dark imgbg" id="info">
                    <div class="container my-auto">
                        <div class="row">
                            <div class="col-lg-10 mx-auto">
                                <h1 class="text-uppercase">Analysis Report</h1>
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
                """,
            'page-bottom':"""
            	<section class="bg-dark text-white" id="foot">
                  <div class="container">
                    <div class="row">
                      <div class="col-lg-8 mx-auto text-center">
                        <p >Generated by <a href="https://github.com/njzjz/ReacNetGenerator">ReacNetGenerator</a></p>
                        <p>Author: <a href="https://cv.njzjz.win/">Jinzhe Zeng</a></p>
                        <p>Email: <a href="mailto:jzzeng@stu.ecnu.edu.cn">jzzeng@stu.ecnu.edu.cn</a></p>
                      </div>
                    </div>
                  </div>
                </section>
                <script src="https://lib.baomitu.com/jquery/3.2.1/jquery.min.js"></script>
                <script src="https://lib.baomitu.com/twitter-bootstrap/4.1.3/js/bootstrap.bundle.min.js"></script>
                <script src="https://lib.baomitu.com/jquery-easing/1.4.1/jquery.easing.min.js"></script>
                <script src="https://lib.baomitu.com/scrollReveal.js/4.0.0-beta.31/scrollreveal.min.js"></script>
                <script src="https://lib.baomitu.com/magnific-popup.js/1.1.0/jquery.magnific-popup.min.js"></script>
                <script src="https://lib.baomitu.com/startbootstrap-creative/5.0.1/js/creative.min.js"></script>
                </body></html>
                """,
            'network':"""<section id="network" class='bg-white mx-auto text-center'>
                    <div class="container my-auto">
                        <div class="row">
                            <div class="mx-auto w-100 text-center">
                                <h2>Reaction Network</h2>
                                <img src='%s' class="w-100">
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
            'species-each':"""%s<td><img src='%s' class='spec'><br/>%s</td>%s
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
            'reactions-each':"""%s<td>%s</td><td><img src='%s' class='spec'></td>
                <td class="reacnum">%s<br/><img src='%s'></td>
                <td><img src='%s' class='spec'></td>
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
            'narrowurl':'https://upload.wikimedia.org/wikipedia/commons/8/8d/U%2B2192.svg',
            
        }
