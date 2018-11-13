#!/usr/bin/env python3
# -*- coding: UTF-8 -*-  
###################################
## Reaction Network Generator(ReacNetGenerator)
## An automatic generator of reaction network for reactive molecular dynamics simulation.
###################################
_version='1.2.11'
_date='2018/11/13'
_author='Jinzhe Zeng'
_email='jzzeng@stu.ecnu.edu.cn'
#########     Features    #########
## * Processing of MD trajectory containing atomic coordinates or bond orders
## * Hidden Markov Model (HMM) based noise filtering
## * Isomers identifying accoarding to SMILES
## * Generation of reaction network for visualization using force-directed algorithm
## * Parallel computing
#########  Simple example #########
## Process a LAMMPS bond file named bonds.reaxc. (See http://lammps.sandia.gov/doc/fix_reax_bonds.html for details)
## >>> from ReacNetGenerator import ReacNetGenerator
## >>> ReacNetGenerator(inputfiletype="lammpsbondfile",inputfilename="bonds.reaxc",atomname=["C","H","O"]).runanddraw()
###################################
#########       Script    #########
###################################
######### import #########
import time
import os
import gc
import math
from multiprocessing import Pool, Semaphore, cpu_count
import numpy as np
from functools import reduce
import itertools
try:
    from hmmlearn import hmm
    hmmlearn_installed=True
except ImportError as e:
    hmmlearn_installed=False
try:
    import networkx as nx
    import networkx.algorithms.isomorphism as iso
    networkx_installed=True
except ImportError as e:
    networkx_installed=False
try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib_installed=True
except ImportError as e:
    matplotlib_installed=False
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    rdkit_installed=True
except ImportError as e:
    rdkit_installed=False
from ReacNetGenerator.reachtml import HTMLResult
######## class ########
class ReacNetGenerator(object):
    def __init__(self,inputfiletype="lammpsbondfile",inputfilename="bonds.reaxc",atomname=["C","H","O"],selectatoms=None,originfilename=None,hmmfilename=None,atomfilename=None,moleculefilename=None,atomroutefilename=None,reactionfilename=None,tablefilename=None,moleculetempfilename=None,moleculetemp2filename=None,moleculestructurefilename=None,imagefilename=None,speciesfilename=None,resultfilename=None,stepinterval=1,p=[0.5,0.5],a=[[0.999,0.001],[0.001,0.999]],b=[[0.6, 0.4],[0.4, 0.6]],runHMM=True,SMILES=True,getoriginfile=False,species={},node_size=200,font_size=6,widthcoefficient=1,show=False,maxspecies=20,n_color=256,drawmolecule=False,nolabel=False,needprintspecies=True,filter=[],node_color=[78/256,196/256,238/256],pos={},printfiltersignal=False,showid=True,k=None,start_color=[0,0,1],end_color=[1,0,0],nproc=None):
        self.version=_version
        self.author=_author
        self.email=_email
        self.date=_date
        self.logging("======= ReacNetGenerator",self.version," ======")
        self.logging("Author:",self.author)
        self.logging("Email:",self.email)
        self.logging()
        self.inputfiletype=inputfiletype
        self.inputfilename=inputfilename
        self.atomname=atomname
        self.selectatoms=selectatoms if selectatoms else self.atomname
        self.originfilename=originfilename if originfilename else inputfilename+".origin"
        self.hmmfilename=(hmmfilename if hmmfilename else inputfilename+".hmm" )if runHMM else self.originfilename
        self.atomfilename=atomfilename if atomfilename else inputfilename+".atom"
        self.moleculefilename=moleculefilename if moleculefilename else inputfilename+".moname"
        self.atomroutefilename=atomroutefilename if atomroutefilename else inputfilename+".route"
        self.reactionfilename=reactionfilename if reactionfilename else inputfilename+".reaction"
        self.tablefilename=tablefilename if tablefilename else inputfilename+".table"
        self.moleculetempfilename=moleculetempfilename if moleculetempfilename else inputfilename+".temp"
        self.moleculetemp2filename=moleculetemp2filename if moleculetemp2filename else inputfilename+".temp2"
        self.moleculestructurefilename=moleculestructurefilename if moleculestructurefilename else inputfilename+".structure"
        self.imagefilename=imagefilename if imagefilename else inputfilename+".svg"
        self.speciesfilename=speciesfilename if speciesfilename else inputfilename+".species"
        self.resultfilename=resultfilename if resultfilename else inputfilename+".html"
        self.stepinterval=stepinterval
        self.p=np.array(p)
        self.a=np.array(a)
        self.b=np.array(b)
        self.runHMM=runHMM
        self.SMILES=SMILES
        self.getoriginfile=getoriginfile
        self.species=species
        self.needprintspecies=needprintspecies
        self.node_size=node_size
        self.font_size=font_size
        self.widthcoefficient=widthcoefficient
        self.show=show
        self.maxspecies=maxspecies
        self.n_color=n_color
        self.drawmolecule=drawmolecule
        self.nolabel=nolabel
        self.filter=filter
        self.node_color=np.array(node_color)
        self.pos=pos
        self.printfiltersignal=printfiltersignal
        self.showid=showid
        self.k=k
        self.start_color=np.array(start_color)
        self.end_color=np.array(end_color)
        self.nproc=nproc if nproc else cpu_count()
        
    #### run and draw ####
    def runanddraw(self,run=True,draw=True,report=True):
        if run:
            self.run()
        if draw:
            self.draw()
        if report:
            self.report()

    ######## run ########
    def run(self):
        if self.SMILES:
            if not rdkit_installed:
                if networkx_installed:
                    self.SMILES=False
                    self.logging("RDkit is not installed. SMILES cannot be used to identify isomers. VF2 algorithm is used.")
                else:
                    self.logging("RDKit is not installed. If you have installed Anaconda, try \"conda install rdkit -c rdkit\".")
                    raise InstallError("RDkit")
        else:
            if not networkx_installed:
                if rdkit_installed:
                    self.SMILES=True
                    self.logging("networkx is not installed. VF2 cannot be used to identify isomers. SMILES algorithm is used.")
                else:
                    self.logging("networkx is not installed. Try \"pip install networkx\".")
                    raise InstallError("networkx")

        if self.runHMM and not hmmlearn_installed:
            self.runHMM=False
            self.logging("Since you have not installed hmmlearn, HMM cannot be used to filter noise.")

        ######start#####
        self.logging("======Run ReacNetGenerator:======")
        timearray=self.printtime([])
        for runstep in range(1,5):
            ######## step 1 ##### 
            if(runstep==1):
                self.step1()
            ######## step 2 ##### 
            elif(runstep==2):
                self.step2()
            ######## step 3 ##### 
            elif(runstep==3):
                allmoleculeroute=self.step3()
            ######## step 4 ##### 
            elif(runstep==4):
                self.step4(allmoleculeroute)
            #garbage collect
            gc.collect()
            timearray=self.printtime(timearray)
        ####### end #######
        self.logging()
        self.logging("Time consumed:")
        for i in range(1,len(timearray)):
            self.logging("Step %d consumed: %.3f s"%(i,timearray[i]-timearray[i-1]))
        self.logging("Total time: %.3f s"%(timearray[-1]-timearray[0]))
        self.logging()

    #####draw#####    
    def draw(self):
        #start
        if not matplotlib_installed:
            self.logging("You must install matplotlib if you want to draw the reaction network. Try \"pip install matplotlib\" to install it.")
            raise InstallError("matplotlib")
        elif not networkx_installed:
            self.logging("You must install networkx if you want to draw the reaction network. Try \"pip install networkx\" to install it.")
            raise InstallError("networkx")
        self.logging("======Draw the image:======")
        timearray=self.printtime([])
        #read table
        table,name=self.readtable()
        species,showname=self.handlespecies(name)

        #make color
        colorsRGB=[(self.start_color+i*(self.end_color-self.start_color) / self.n_color) for i in np.arange(self.n_color)]

        G = nx.DiGraph()
        for i in range(len(table)):
            if name[i] in species and not name[i] in self.filter:
                G.add_node(showname[name[i]] if name[i] in showname else name[i])
                for j in range(len(table)):
                    if name[j] in species and not name[j] in self.filter:
                        if table[i][j]>0:
                            G.add_weighted_edges_from([((showname[name[i]] if name[i] in showname else name[i]),(showname[name[j]] if name[j] in showname else name[j]),table[i][j])])
        weights = np.array([math.log(G[u][v]['weight']+1) for u,v in G.edges()])
        #widths=[weight/max(weights) *self.widthcoefficient for weight in weights]
        widths=[weight/max(weights) *self.widthcoefficient*2 if weight>max(weights)*0.7 else weight/max(weights) *self.widthcoefficient*0.5 for weight in weights]
        colors=[colorsRGB[math.floor(weight/max(weights)*(self.n_color-1))] for weight in weights]
        try:
            pos = (nx.spring_layout(G) if not self.pos else nx.spring_layout(G,pos=self.pos,fixed=[p for p in self.pos])) if not self.k else (nx.spring_layout(G,k=self.k) if not self.pos else nx.spring_layout(G,pos=self.pos,fixed=[p for p in self.pos],k=self.k))
            self.logging()
            self.logging("The position of the species in the network is:")
            self.logging(pos)
            self.logging()

            for with_labels in ([True] if not self.nolabel else [True,False]):
                nx.draw(G,pos = pos,width=widths,node_size=self.node_size,font_size=self.font_size,with_labels=with_labels,edge_color=colors,node_color=self.node_color)

                plt.savefig(self.imagefilename if with_labels else "nolabel_"+self.imagefilename)

                if self.show:
                    plt.show()

                plt.close()

        except Exception as e:
            self.logging("Error: cannot draw images. Details:",e)

        timearray=self.printtime(timearray)
        ####### end #######
        self.logging()
        self.logging("Time consumed:")
        for i in range(1,len(timearray)):
            self.logging("Step %d consumed: %.3f s"%(i,timearray[i]-timearray[i-1]))
        self.logging("Total time: %.3f s"%(timearray[-1]-timearray[0]))
        self.logging()
        return pos

    ######## report #####
    def report(self):
        self.logging("======Report:======")
        timearray=self.printtime([])
        HTMLResult(reactionfile=self.reactionfilename,resultfile=self.resultfilename,imagefile=self.imagefilename,n_thread=self.nproc).report()
        timearray=self.printtime(timearray)
        ####### end #######
        self.logging()
        self.logging("Time consumed:")
        for i in range(1,len(timearray)):
            self.logging("Step %d consumed: %.3f s"%(i,timearray[i]-timearray[i-1]))
        self.logging("Total time: %.3f s"%(timearray[-1]-timearray[0]))
        self.logging()
        self.logging("Please view %s for more details."%self.reactionfilename)
        
    ######## steps ######
    def step1(self):
        if self.inputfiletype=="lammpsbondfile":
            readNfunc=self.readlammpsbondN
            readstepfunc=self.readlammpsbondstep
        elif self.inputfiletype=="lammpscrdfile" or self.inputfiletype=="lammpsdumpfile":
            readNfunc=self.readlammpscrdN
            readstepfunc=self.readlammpscrdstep
        self.readinputfile(readNfunc,readstepfunc)

    def step2(self):
        if self.runHMM:
            self.initHMM()
        self.calhmm()

    def step3(self):
        if self.SMILES:
            self.printmoleculeSMILESname()
        else:
            self.printmoleculename()
        atomeach=self.getatomeach()
        allmoleculeroute=self.printatomroute(atomeach)
        return allmoleculeroute

    def step4(self,allmoleculeroute):
        allroute=self.getallroute(allmoleculeroute)
        self.printtable(allroute)
        if self.needprintspecies:
            self.printspecies()

    ####### functions #######
    def logging(self,*message):
        if message:
            localtime = time.asctime( time.localtime(time.time()) )
            print(localtime,'ReacNetGenerator',self.version,*message)
        else:
            print()

    def printtime(self,timearray):
        timearray.append(time.time())
        if len(timearray)>1:
            self.logging("Step %d has been completed. Time consumed: %f s"%(len(timearray)-1,timearray[-1]-timearray[-2]))
        return timearray

    def union_dict(self,x,y):
        for k, v in y.items():
            if k in x.keys():
                x[k] += v
            else:
                x[k] = v
        return x

    def mo(self,i,bond,level,molecule,done,bondlist): #connect molecule
        molecule.append(i)
        done[i]=True
        for j in range(len(bond[i])):
            b=bond[i][j]
            l=level[i][j]
            bo=(i,b,l) if i<b else (b,i,l)
            if not bo in bondlist:
                bondlist.append(bo)
            if not done[b]:
                molecule,done,bondlist=self.mo(b,bond,level,molecule,done,bondlist)
        return molecule,done,bondlist

    def readinputfile(self,readNfunc,readstepfunc):
        steplinenum=readNfunc()
        self.getdandtimestep(readstepfunc,steplinenum)

    def readlammpsbondN(self):
        with open(self.inputfilename) as file:
            iscompleted=False
            for index,line in enumerate(file):
                if line.startswith("#"):
                    if line.startswith("# Number of particles"):
                        if iscompleted:
                            stepbindex=index
                            break
                        else:
                            iscompleted=True
                            stepaindex=index
                        N=[int(s) for s in line.split() if s.isdigit()][0]
                        atomtype=np.zeros(N+1,dtype=np.int)
                else:
                    s=line.split()
                    atomtype[int(s[0])]=int(s[1])
        steplinenum=stepbindex-stepaindex
        self.N=N
        self.atomtype=atomtype
        return steplinenum

    def readlammpsbondstep(self,item):
        (step,lines),_=item
        bond=[[] for x in range(self.N+1)]
        level=[[] for x in range(self.N+1)]
        for line in lines:
            if line:
                if line.startswith("#"):
                    if line.startswith("# Timestep"):
                        timestep=step,[int(s) for s in line.split() if s.isdigit()][0]
                else:
                    s=line.split()
                    for i in range(int(s[2])):
                        bond[int(s[0])].append(int(s[i+3]))
                        bondlevel=round(float(s[i+4+int(s[2])]))
                        if bondlevel==0:
                            bondlevel=1
                        level[int(s[0])].append(bondlevel)
        d=self.connectmolecule({},step,bond,level)
        return d,timestep

    def readlammpscrdN(self):
        with open(self.inputfilename) as f:
            iscompleted=False
            for index,line in enumerate(f):
                if line.startswith("ITEM:"):
                    if line.startswith("ITEM: TIMESTEP"):
                        linecontent=4
                    elif line.startswith("ITEM: ATOMS"):
                        linecontent=3
                    elif line.startswith("ITEM: NUMBER OF ATOMS"):
                        linecontent=1
                    elif line.startswith("ITEM: BOX BOUNDS"):
                        linecontent=2
                else:
                    if linecontent==1:
                        if iscompleted:
                            stepbindex=index
                            break
                        else:
                            iscompleted=True
                            stepaindex=index
                        N=int(line.split()[0])
                        atomtype=np.zeros(N+1,dtype=np.int)
                    elif linecontent==3:
                        s=line.split()
                        atomtype[int(s[0])]=int(s[1])
        steplinenum=stepbindex-stepaindex
        self.N=N
        self.atomtype=atomtype
        return steplinenum

    def readlammpscrdstep(self,item):
        (step,lines),_=item
        atomtype=np.zeros((self.N),dtype=np.int)
        atomcrd=np.zeros((self.N,3))
        for line in lines:
            if line:
                if line.startswith("ITEM:"):
                    if line.startswith("ITEM: TIMESTEP"):
                        linecontent=4
                    elif line.startswith("ITEM: ATOMS"):
                        linecontent=3
                    elif line.startswith("ITEM: NUMBER OF ATOMS"):
                        linecontent=1
                    elif line.startswith("ITEM: BOX BOUNDS"):
                        linecontent=2
                else:
                    if linecontent==3:
                        s=line.split()
                        atomtype[int(s[0])-1]=int(s[1])
                        atomcrd[int(s[0])-1]=float(s[2]),float(s[3]),float(s[4])
                    elif linecontent==4:
                        timestep=step,int(line.split()[0])
        bond,level=self.getbondfromcrd(atomtype,atomcrd,step)
        d=self.connectmolecule({},step,bond,level)
        return d,timestep

    def getdandtimestep(self,readfunc,steplinenum):
        d={}
        timestep={}
        with open(self.inputfilename) as file,Pool(self.nproc,maxtasksperchild=100) as pool:
            semaphore = Semaphore(360)
            results=pool.imap_unordered(readfunc,self.produce(semaphore,enumerate(itertools.islice(itertools.zip_longest(*[file]*steplinenum),0,None,self.stepinterval)),None),10)
            for dstep,timesteptuple in results:
                d=self.union_dict(d,dstep)
                step,thetimestep=timesteptuple
                timestep[step]=thetimestep
                semaphore.release()
        self.writemoleculetempfile(d)
        self.timestep=timestep
        self.step=len(timestep)-1

    def connectmolecule(self,d,step,bond,level):
        #init
        done=np.zeros(self.N+1,dtype=bool)
        #connect molecule
        for i in range(1,self.N+1):
            if not done[i]:
                mole,done,bondlist=self.mo(i,bond,level,[],done,[])
                mole.sort()
                bondlist.sort()
                if (tuple(mole),tuple(bondlist)) in d:
                    d[(tuple(mole),tuple(bondlist))].append(step)
                else:
                    d[(tuple(mole),tuple(bondlist))]=[step]
        return d

    def writemoleculetempfile(self,d):
        with open(self.moleculetempfilename,'w') as f:
            for item in d.items():
                key,value=item
                print(",".join([str(x) for x in key[0]]),";".join([",".join([str(y) for y in x]) for x in key[1]]),",".join([str(x) for x in value]),file=f)

    def getbondfromcrd(self,atomtype,atomcrd,step,filename="crd"):
        xyzfilename=filename+"_"+str(step)+".xyz"
        mol2filename=filename+"_"+str(step)+".mol2"
        self.convertxyz(atomtype,atomcrd,xyzfilename)
        os.system("obabel -ixyz "+xyzfilename+" -omol2 -O "+mol2filename+" >/dev/null")
        bond,bondlevel=self.getbondfrommol2(len(atomcrd),mol2filename)
        return bond,bondlevel

    def convertxyz(self,atomtype,atomcrd,xyzfilename):
        with open(xyzfilename,'w') as f:
            print(len(atomcrd),file=f)
            print("by ReacNetGenerator",file=f)
            for type,(x,y,z) in zip(atomtype,atomcrd):
                print(self.atomname[type-1],x,y,z,file=f)

    def getbondfrommol2(self,atomnumber,mol2filename):
        linecontent=-1
        bond=[[] for i in range(atomnumber+1)]
        bondlevel=[[] for i in range(atomnumber+1)]
        with open(mol2filename) as f:
            for line in f:
                if line.startswith("@<TRIPOS>BOND"):
                    linecontent=0
                else:
                    if linecontent==0:
                        s=line.split()
                        bond[int(s[1])].append(int(s[2]))
                        bond[int(s[2])].append(int(s[1]))
                        level=12 if s[3]=='ar' else int(s[3])
                        bondlevel[int(s[1])].append(level)
                        bondlevel[int(s[2])].append(level)
        return bond,bondlevel

    def initHMM(self):
        self.model = hmm.MultinomialHMM(n_components=2)
        self.model.startprob_= self.p
        self.model.transmat_= self.a
        self.model.emissionprob_= self.b

    def produce(self,semaphore, list,parameter):
        for item in list:
            # Reduce Semaphore by 1 or wait if 0
            semaphore.acquire()
            # Now deliver an item to the caller (pool)
            yield item,parameter

    def getoriginandhmm(self,item):
        line,_=item
        list=line.split()
        value=np.array([int(x)-1 for x in list[-1].split(",")])
        origin=np.zeros(self.step,dtype=np.int)
        origin[value]=1
        if self.runHMM:
            logprob,hmm=self.model.decode(np.array([origin]).T,algorithm="viterbi")
            return origin,np.array(hmm),line
        else:
            return origin,line

    def calhmm(self):
        with open(self.originfilename, 'w') if self.getoriginfile or not self.runHMM else Placeholder() as fo,open(self.hmmfilename, 'w') if self.runHMM else Placeholder() as fh,open(self.moleculetempfilename) as ft,open(self.moleculetemp2filename,'w') as ft2,Pool(self.nproc,maxtasksperchild=100) as pool:
            semaphore = Semaphore(360)
            results=pool.imap_unordered(self.getoriginandhmm,self.produce(semaphore,ft,()),10)
            if self.runHMM:
                for originsignal,hmmsignal,mlist in results:
                    if 1 in hmmsignal or self.printfiltersignal:
                        if self.getoriginfile:
                            print("".join([str(i) for i in originsignal]), file=fo)
                        print("".join([str(i) for i in hmmsignal]), file=fh)
                        print(mlist,end='',file=ft2)
                    semaphore.release()
            else:
                for originsignal,mlist in results:
                    print("".join([str(i) for i in originsignal]), file=fo)
                    print(mlist,end='',file=ft2)
                    semaphore.release()

    def getatomroute(self,item):
        (i,(atomeachi,atomtypei)),_=item
        route=[]
        routestrarr=[]
        moleculeroute=[]
        molecule=-1
        right=-1
        for j in range(0,self.step):
            if atomeachi[j]>0 and atomeachi[j]!=molecule:
                routestrarr.append("%s (%d step %d)"%(self.mname[atomeachi[j]-1],atomeachi[j],self.timestep[j]))
                left=right
                molecule=atomeachi[j]
                right=molecule
                if self.atomname[atomtypei-1] in self.selectatoms:
                    if left>=0 and not (left,right) in moleculeroute:
                        moleculeroute.append((left,right))
        routestr="Atom %d %s: "%(i,self.atomname[atomtypei-1])+" -> ".join(routestrarr)
        return moleculeroute,routestr

    def printatomroute(self,atomeach):
        with open(self.atomroutefilename, 'w') as f,Pool(self.nproc,maxtasksperchild=100) as pool:
            allmoleculeroute=[]
            semaphore = Semaphore(360)
            results=pool.imap(self.getatomroute,self.produce(semaphore,enumerate(zip(atomeach[1:],self.atomtype[1:]),start=1),()),10)
            for route in results:
                moleculeroute,routestr=route
                print(routestr, file=f)
                for mroute in moleculeroute:
                    if not mroute in allmoleculeroute:
                        allmoleculeroute.append(mroute)
                semaphore.release()
        return allmoleculeroute

    def makemoleculegraph(self,atoms,bonds):
        G=nx.Graph()
        for line in bonds:
            G.add_edge(line[0],line[1],level=line[2])
        for atom in atoms:
            atomnumber,atomtype=atom
            G.add_node(atomnumber, atom=atomtype)
        return G

    def getstructure(self,name,atoms,bonds):
        index={}
        for i,atom in enumerate(atoms,start=1):
            index[atom]=i
        return name+" "+",".join([self.atomname[self.atomtype[x]-1] for x in atoms])+" "+";".join([str(index[x[0]])+","+str(index[x[1]])+","+str(x[2]) for x in bonds])

    def readstrcture(self):
        with open(self.moleculestructurefilename) as f:
            d={}
            for line in f:
                list=line.split()
                name=list[0]
                atoms=[x for x in list[1].split(",")]
                bonds=[tuple(int(y) for y in x.split(",")) for x in list[2].split(";")] if len(list)==3 else []
                d[name]=(atoms,bonds)
        return d

    def printmoleculename(self):
        mname=[]
        d={}
        em = iso.numerical_edge_match(['atom','level'], ["None",1])
        with open(self.moleculefilename, 'w') as fm,open(self.moleculetemp2filename) as ft,open(self.moleculestructurefilename,'w') as fs:
            for line in ft:
                list=line.split()
                atoms=np.array([int(x) for x in list[0].split(",")])
                bonds=np.array([tuple(int(y) for y in x.split(",")) for x in list[1].split(";")] if len(list)==3 else [])
                typenumber=np.zeros(len(self.atomname),dtype=np.int)
                atomtypes=[]
                for atomnumber in atoms:
                    typenumber[self.atomtype[atomnumber]-1]+=1
                    atomtypes.append((atomnumber,self.atomtype[atomnumber]))
                G=self.makemoleculegraph(atomtypes,bonds)
                name="".join([self.atomname[i]+(str(typenumber[i] if typenumber[i]>1 else "")) if typenumber[i]>0 else "" for i in range(0,len(self.atomname))])
                if name in d:
                    for j in range(len(d[name])):
                        if nx.is_isomorphic(G,d[name][j],em):
                            if j>0:
                                name+="_"+str(j+1)
                            break
                    else:
                        d[name].append(G)
                        name+="_"+str(len(d[name]))
                        print(self.getstructure(name,atoms,bonds),file=fs)
                else:
                    d[name]=[G]
                    print(self.getstructure(name,atoms,bonds),file=fs)
                mname.append(name)
                print(name,",".join([str(x) for x in atoms]),";".join([",".join([str(y) for y in x]) for x in bonds]), file=fm)
        self.mname=mname

    def calmoleculeSMILESname(self,item):
        line,_=item
        list=line.split()
        atoms=np.array([int(x) for x in list[0].split(",")])
        bonds=np.array([tuple(int(y) for y in x.split(",")) for x in list[1].split(";")] if len(list)==3 else [])
        type={}
        for atomnumber in atoms:
            type[atomnumber]=self.atomname[self.atomtype[atomnumber]-1]
        name=self.convertSMILES(atoms,bonds,type)
        return name,atoms,bonds

    def printmoleculeSMILESname(self):
        mname=[]
        with open(self.moleculefilename, 'w') as fm,open(self.moleculetemp2filename) as ft,Pool(self.nproc,maxtasksperchild=100) as pool:
            semaphore = Semaphore(360)
            results=pool.imap(self.calmoleculeSMILESname,self.produce(semaphore,ft,()),10)
            for result in results:
                name,atoms,bonds=result
                mname.append(name)
                print(name,",".join([str(x) for x in atoms]),";".join([",".join([str(y) for y in x]) for x in bonds]),file=fm)
                semaphore.release()
        self.mname=mname

    def convertSMILES(self,atoms,bonds,type):
        m = Chem.RWMol(Chem.MolFromSmiles(''))
        d={}
        for atomnumber in atoms:
            d[atomnumber]=m.AddAtom(Chem.Atom(type[atomnumber]))
        for bond in bonds:
            atom1,atom2,level=bond
            m.AddBond(d[atom1],d[atom2], Chem.BondType.DOUBLE if level==2 else (Chem.BondType.TRIPLE if level==3 else (Chem.BondType.AROMATIC if level==12 else Chem.BondType.SINGLE)))
        name=Chem.MolToSmiles(m)
        return name

    def getatomeach(self):
        atomeach=np.zeros((self.N+1,self.step),dtype=np.int)
        with open(self.hmmfilename) as fh,open(self.moleculetemp2filename) as ft:
            for i,(lineh,linet) in enumerate(zip(fh,ft),start=1):
                list=linet.split()
                key1=np.array([int(x) for x in list[0].split(",")])
                index=np.array([j for j in range(len(lineh)) if lineh[j]=="1"])
                if(len(index))>0:
                    atomeach[key1[:,None],index]=i
        with open(self.atomfilename, 'w') as f:
            for atom in atomeach[1:]:
                print(atom, file=f)
        return atomeach

    def getallroute(self,allmoleculeroute):
        allroute={}
        for moleculeroute in allmoleculeroute:
            leftname=self.mname[moleculeroute[0]-1]
            rightname=self.mname[moleculeroute[1]-1]
            if leftname==rightname:
                continue
            equation=leftname+"->"+rightname
            if equation in allroute:
                allroute[equation]+=1
            else:
                allroute[equation]=1
        return allroute

    def printtable(self,allroute):
        species=[]
        table=np.zeros((100,100),dtype=np.int)
        reactionnumber=np.zeros((2),dtype=np.int)
        with open(self.reactionfilename,'w') as f:
            for k, v in sorted(allroute.items(), key=lambda d: d[1] ,reverse=True):
                print(v,k,file=f)
                reaction=k.split("->")
                for i,spec in enumerate(reaction):
                    if spec in species:
                        number=species.index(spec)
                    elif len(species)<100:
                        species.append(spec)
                        number=species.index(spec)
                    else:
                        number=-1
                    reactionnumber[i]=number
                if all(reactionnumber>=0):
                    table[reactionnumber[0]][reactionnumber[1]]=v
        with open(self.tablefilename,'w') as f:
            print("\t"+"\t".join(species),file=f)
            for i in range(len(species)):
                print(species[i],end='\t',file=f)
                for j in range(len(species)):
                    print(table[i][j],end='\t',file=f)
                print(file=f)

    def readtable(self):
        table=[]
        name=[]
        with open(self.tablefilename) as file:
            for line in itertools.islice(file, 1, None):
                name.append(line.split()[0])
                table.append([int(s) for s in line.split()[1:]])
        return table,name

    def convertstructure(self,atoms,bonds):
        types={}
        atomtypes=[]
        for atom in enumerate(atoms,start=1):
            atomtypes.append((i,self.atomname.index(atom)))
        G=self.makemoleculegraph(atomtypes,bonds)
        return G

    def handlespecies(self,name):
        showname={}
        n=0
        if self.species=={}:
            species_out=dict([(x,{}) for x in (name if len(name)<=self.maxspecies else name[0:self.maxspecies])])
        else:
            species_out={}
            b=True
            for spec in self.species.items():
                specname,value=spec
                if "structure" in value:
                    atoms,bonds=value["structure"]
                    G1=convertstructure(atoms,bonds)
                    if b:
                        structures=readstrcture()
                        em = iso.numerical_edge_match(['atom','level'], ["None",1])
                        b=False
                    i=1
                    while (specname+"_"+str(i) if i>1 else specname) in structures:
                        G2=convertstructure(structures[(specname+"_"+str(i) if i>1 else specname)][0],structures[(specname+"_"+str(i) if i>1 else specname)][1])
                        if nx.is_isomorphic(G1,G2,em):
                            if i>1:
                                specname+="_"+str(i)
                            break
                        i+=1
                species_out[specname]={}
                if "showname" in value:
                    showname[specname]=value["showname"]
        if self.showid:
            print()
            self.logging("Species are:")
            for specname,value in species_out.items():
                n+=1
                showname[specname]=str(n)
                print(n,specname)
        return species_out,showname

    def printspecies(self):
        with open(self.moleculetemp2filename) as f2,open(self.speciesfilename,'w') as fw:
            d=[{} for i in range(len(self.timestep))]
            for name,line2 in zip(self.mname,f2):
                for t in [int(x) for x in line2.split()[-1].split(",")]:
                    if name in d[t]:
                        d[t][name]+=1
                    else:
                        d[t][name]=1
            for t in range(len(self.timestep)):
                print("Timestep",self.timestep[t],":",end=' ',file=fw)
                for name,num in d[t].items():
                    print(name,num,end=' ',file=fw)
                print(file=fw)


    def __enter__(self):return self
    def __exit__(self,Type, value, traceback):pass

class InstallError(Exception):
    def __init__(self,ErrorInfo):
        super().__init__(self)
        self.errorinfo=ErrorInfo
    def __str__(self):
        return self.errorinfo+" is not installed."

class Placeholder(object):
    def __init__(self):pass
    def __enter__(self):return self
    def __exit__(self,Type, value, traceback):pass
