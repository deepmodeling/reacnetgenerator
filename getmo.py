# -*- coding: UTF-8 -*-  
# updated at 2018/4/20 16:00
#########  Usage #########
## import getmo
## getmo.run()
## getmo.draw()
##########################
######### import #########
import time
import gc
import math
from multiprocessing import Pool, Semaphore
import numpy as np
from hmmlearn import hmm
from functools import reduce
import networkx as nx
import networkx.algorithms.isomorphism as iso
from itertools import islice
import matplotlib.pyplot as plt
######## function ########
def printtime(timearray):
    timearray.append(time.time())
    if len(timearray)>1:
       print("Step ",len(timearray)-1," has been completed. Time consumed: ",round(timearray[-1]-timearray[-2],3),"s")
    return timearray

def readatomnumber(bondfilename):
    with open(bondfilename) as file: #read Atom number
        for line in file:
            if line.startswith("# Number of particles"):
                N=[int(s) for s in line.split() if s.isdigit()][0]
                break
    return N
    
def mo(i,bond,level,molecule,done,bondlist): #connect molecule
    molecule.append(i)
    done[i]=True
    for j in range(len(bond[i])):
        b=bond[i][j]
        l=level[i][j]
        bo=(i,b,l) if i<b else (b,i,l)
        if not bo in bondlist:
            bondlist.append(bo)
        if not done[b]:
            molecule,done,bondlist=mo(b,bond,level,molecule,done,bondlist)
    return molecule,done,bondlist

def openbondfile(bondfilename,moleculetempfilename,N,stepinterval):
    step=0
    realstep=0
    iscontent=False
    isfirststep=True
    timestep=[]
    atomtype=[0 for x in range(N+1)]
    d={}
    with open(bondfilename) as file:
        for line in file:
            if line.startswith("#"):
                if iscontent and realstep%stepinterval==0:
                    #a step finished
                    iscontent=False
                    if isfirststep:
                        isfirststep=False
                    #init
                    done=[False for x in range(N+1)]
                    #connect molecule
                    for i in range(1,N+1):
                        if not done[i]:
                            mole,done,bondlist=mo(i,bond,level,[],done,[])
                            mole.sort()
                            bondlist.sort()
                            if (tuple(mole),tuple(bondlist)) in d:
                                d[(tuple(mole),tuple(bondlist))].append(step)
                            else:
                                d[(tuple(mole),tuple(bondlist))]=[step]
                if line.startswith("# Timestep"):
                    #read step
                    realstep+=1
                    if realstep%stepinterval==0:
                        step+=1
                        timestep.append([int(s) for s in line.split() if s.isdigit()][0])
            else:
                if not realstep%stepinterval==0:
                    continue
                if not iscontent:
                    #new step start
                    iscontent=True
                    #init
                    bond=[[] for x in range(N+1)]
                    level=[[] for x in range(N+1)]
                #read
                s=line.split()
                for i in range(0,int(s[2])):
                    bond[int(s[0])].append(int(s[i+3]))
                    bondlevel=round(float(s[i+4+int(s[2])]))
                    if bondlevel==0:
                        bondlevel=1
                    level[int(s[0])].append(bondlevel)
                #read atom type
                if isfirststep:
                    atomtype[int(s[0])]=int(s[1])
    with open(moleculetempfilename,'w') as f:
        for item in d.items():
            key,value=item
            print(",".join([str(x) for x in key[0]]),";".join([",".join([str(y) for y in x]) for x in key[1]]),",".join([str(x) for x in value]),file=f)
    return atomtype,step,timestep

def initHMM(states,observations,p,a,b):
    n_states = len(states)
    n_observations = len(observations)
    model = hmm.MultinomialHMM(n_components=n_states)
    model.startprob_= np.array(p)
    model.transmat_= np.array(a)
    model.emissionprob_= np.array(b)
    return model

def gethmm(ori,model,states):
    o = np.array([ori]).T
    logprob, h = model.decode(o, algorithm="viterbi")
    hmmlist=list(map(lambda x: states[x], h))
    return hmmlist
            
def produce(semaphore, list,parameter):
    for item in list:
        # Reduce Semaphore by 1 or wait if 0
        semaphore.acquire()
        # Now deliver an item to the caller (pool)
        yield item,parameter

def getoriginandhmm(item):
    line,parameter=item
    step,model,states=parameter
    list=line.split()
    value=[int(x) for x in list[-1].split(",")]
    origin= [1 if i in value else 0 for i in range(1, step+1)]
    hmm=gethmm(origin,model,states)
    return origin,hmm,line

def calhmm(originfilename,hmmfilename,moleculetempfilename,moleculetemp2filename,model,states,step,getoriginfile):
    with open(originfilename, 'w') as fo,open(hmmfilename, 'w') as fh,open(moleculetempfilename) as ft,open(moleculetemp2filename,'w') as ft2,Pool(maxtasksperchild=100) as pool:
        semaphore = Semaphore(360)
        results=pool.imap_unordered(getoriginandhmm,produce(semaphore,ft,(step,model,states)),10)
        for originsignal,hmmsignal,mlist in results:
            if 1 in hmmsignal:
                if getoriginfile:
                    print("".join([str(i) for i in originsignal]), file=fo)
                print("".join([str(i) for i in hmmsignal]), file=fh)
                print(mlist,end='',file=ft2)
            semaphore.release()
     
def getorigin(item):
    line,parameter=item
    step,=parameter
    list=line.split()
    value=[int(x) for x in list[-1].split(",")]
    origin= [1 if i in value else 0 for i in range(1, step+1)]
    return origin,line     
         
def noHMM(originfilename,moleculetempfilename,moleculetemp2filename,step):
    with open(originfilename, 'w') as fh,open(moleculetempfilename) as ft,open(moleculetemp2filename,'w') as ft2,Pool(maxtasksperchild=100) as pool:
        semaphore = Semaphore(360)
        results=pool.imap_unordered(getorigin,produce(semaphore,ft,(step,)),10)
        for originsignal,mlist in results:
            print("".join([str(i) for i in originsignal]), file=fh)
            print(mlist,end='',file=ft2)
            semaphore.release()

def getatomroute(item):
    itemi,parameter=item
    i,atomeachi,atomtypei=itemi
    step,atomname,mname,timestep=parameter
    route=[]
    routestrarr=[]
    moleculeroute=[]
    molecule=-1
    right=-1
    for j in range(0,step):
        if atomeachi[j]>0 and atomeachi[j]!=molecule:
            routestrarr.append(mname[atomeachi[j]-1] + " ("+ str(atomeachi[j])+" step "+str(timestep[j])+")")
            left=right
            molecule=atomeachi[j]
            right=molecule
            if left>=0 and not (left,right) in moleculeroute:
                moleculeroute.append((left,right))
    routestr="Atom "+str(i)+" "+atomname[atomtypei-1]+": "+" -> ".join(routestrarr)
    return moleculeroute,routestr

def printatomroute(atomroutefilename,N,step,atomeach,atomtype,atomname,mname,timestep):
    with open(atomroutefilename, 'w') as f,Pool(maxtasksperchild=100) as pool:
        allmoleculeroute=[]
        semaphore = Semaphore(360)
        results=pool.imap(getatomroute,produce(semaphore,zip(range(1,N+1),atomeach[1:],atomtype[1:]),(step,atomname,mname,timestep)),10)
        for route in results:
            moleculeroute,routestr=route
            print(routestr, file=f)
            for mroute in moleculeroute:
                if not mroute in allmoleculeroute:
                    allmoleculeroute.append(mroute)
            semaphore.release() 
    return allmoleculeroute

def makemoleculegraph(atoms,bonds):
    G=nx.Graph()
    for line in bonds:
        G.add_edge(line[0],line[1],level=line[2])
    for atom in atoms:
        atomnumber,atomtype=atom
        G.add_node(atomnumber, atom=atomtype)
    return G

def getstructure(name,atoms,bonds,atomtype,atomname):
    index={}
    i=1
    for atom in atoms:
        index[atom]=i
        i+=1
    return name+" "+",".join([atomname[atomtype[x]-1] for x in atoms])+" "+";".join([str(index[x[0]])+","+str(index[x[1]])+","+str(x[2]) for x in bonds])

def readstrcture(moleculestructurefilename):
    with open(moleculestructurefilename) as f:
        d={}
        for line in f:
            list=line.split()
            name=list[0]
            atoms=[x for x in list[1].split(",")]
            bonds=[tuple(int(y) for y in x.split(",")) for x in list[2].split(";")] if len(list)==3 else []
            d[name]=(atoms,bonds)
    return d

def printmoleculename(moleculefilename,moleculetempfilename,moleculestructurefilename,atomname,atomtype):
    mname=[]
    d={}
    em = iso.numerical_edge_match(['atom','level'], ["None",1])
    with open(moleculefilename, 'w') as fm,open(moleculetempfilename) as ft,open(moleculestructurefilename,'w') as fs:
        for line in ft:
            list=line.split()
            atoms=[int(x) for x in list[0].split(",")]
            bonds=[tuple(int(y) for y in x.split(",")) for x in list[1].split(";")] if len(list)==3 else []
            typenumber=[0 for i in range(len(atomname))]
            atomtypes=[]
            for atomnumber in atoms:
                typenumber[atomtype[atomnumber]-1]+=1
                atomtypes.append((atomnumber,atomtype[atomnumber]))
            G=makemoleculegraph(atomtypes,bonds)
            name="".join([atomname[i]+(str(typenumber[i] if typenumber[i]>1 else "")) if typenumber[i]>0 else "" for i in range(0,len(atomname))])                
            if name in d:
                for j in range(len(d[name])):
                    if nx.is_isomorphic(G,d[name][j],em):
                        if j>0:
                            name+="_"+str(j+1)
                        break
                else:
                    d[name].append(G)
                    name+="_"+str(len(d[name]))
                    print(getstructure(name,atoms,bonds,atomtype,atomname),file=fs)
            else:
                d[name]=[G]
                print(getstructure(name,atoms,bonds,atomtype,atomname),file=fs)
            mname.append(name)
            print(name,atoms,bonds, file=fm)
    return mname

def getatomeach(hmmfilename,moleculetemp2filename,atomfilename,N,step):
    atomeach=[[0 for j in range(0,step)] for i in range(0,N+1)]
    with open(hmmfilename) as fh,open(moleculetemp2filename) as ft:
        i=0
        for lineh,linet in zip(fh,ft):
            list=linet.split()
            key1=[int(x) for x in list[0].split(",")]
            for j in range(0,step):#j is step
                if lineh[j]=="1":
                    for a in key1:
                        atomeach[a][j]=i+1
            i+=1
    with open(atomfilename, 'w') as f:
        for atom in atomeach[1:]:
            print(atom, file=f)
    return atomeach

def getallroute(reactionfilename,allmoleculeroute,mname):
    allroute={}
    for moleculeroute in allmoleculeroute:
        leftname=mname[moleculeroute[0]-1]
        rightname=mname[moleculeroute[1]-1]
        if leftname==rightname:
            continue
        equation=leftname+"->"+rightname
        if equation in allroute:
            allroute[equation]+=1
        else:
            allroute[equation]=1            
    return allroute

def printtable(tablefilename,reactionfilename,allroute):
    species=[]
    table=[[0 for x in range(100)] for x in range(100)]
    with open(reactionfilename,'w') as f:
        for k, v in sorted(allroute.items(), key=lambda d: d[1] ,reverse=True):
            print(v,k,file=f)
            l=k.split("->")
            left=l[0]
            right=l[1]
            if left in species:
                leftnumber=species.index(left)    
            else:
                if len(species)<100:
                    species.append(left)
                    leftnumber=species.index(left)
                else:
                    leftnumber=-1
            if right in species:
                rightnumber=species.index(right)    
            else:
                if len(species)<100:
                    species.append(right)
                    rightnumber=species.index(right)
                else:
                    rightnumber=-1
            if leftnumber>=0 and rightnumber>=0:
                table[leftnumber][rightnumber]=v
    with open(tablefilename,'w') as f:
        print("\t"+"\t".join(species),file=f)
        for i in range(len(species)):
            print(species[i],end='\t',file=f)
            for j in range(len(species)):
                print(table[i][j],end='\t',file=f)
            print(file=f)

def readtable(tablefilename):
    table=[]
    name=[]
    with open(tablefilename) as file:
        for line in islice(file, 1, None):
            name.append(line.split()[0])
            table.append([int(s) for s in line.split()[1:]])
    return table,name
    
def convertstructure(atoms,bonds,atomname):
    types={}
    i=0
    atomtypes=[]
    for atom in atoms:
        i+=1
        atomtypes.append((i,atomname.index(atom)))
    G=makemoleculegraph(atomtypes,bonds)
    return G
    
def handlespecies(species,name,maxspecies,atomname,moleculestructurefilename):
    showname={}
    if species=={}:
        species_out=dict([(x,{}) for x in (name if len(name)<=maxspecies else name[0:(maxspecies-1)])])
    else:
        species_out={}
        b=True
        for spec in species.items():
            specname,value=spec
            if "structure" in value:
                atoms,bonds=value["structure"]
                G1=convertstructure(atoms,bonds,atomname)
                if b:
                    structures=readstrcture(moleculestructurefilename)
                    em = iso.numerical_edge_match(['atom','level'], ["None",1])
                    b=False
                i=1
                while (specname+"_"+str(i) if i>1 else specname) in structures:
                    G2=convertstructure(structures[(specname+"_"+str(i) if i>1 else specname)][0],structures[(specname+"_"+str(i) if i>1 else specname)][1],atomname)
                    if nx.is_isomorphic(G1,G2,em):
                        if i>1:
                            specname+="_"+str(i)
                        break
                    i+=1
            species_out[specname]={}
            if "showname" in value:
                showname[specname]=value["showname"]
    return species_out,showname
######## steps ######
def step1(bondfilename,moleculetempfilename,stepinterval):
    N=readatomnumber(bondfilename)
    atomtype,step,timestep=openbondfile(bondfilename,moleculetempfilename,N,stepinterval)
    return N,atomtype,step,timestep
    
def step2(states,observations,p,a,b,originfilename,hmmfilename,moleculetempfilename,moleculetemp2filename,step,runHMM,getoriginfile):
    if runHMM:
        model=initHMM(states,observations,p,a,b)
        calhmm(originfilename,hmmfilename,moleculetempfilename,moleculetemp2filename,model,states,step,getoriginfile)
    else:
        noHMM(originfilename,moleculetempfilename,moleculetemp2filename,step)

def step3(atomname,atomtype,N,step,timestep,moleculefilename,hmmfilename,atomfilename,moleculetemp2filename,atomroutefilename,moleculestructurefilename):
    mname=printmoleculename(moleculefilename,moleculetemp2filename,moleculestructurefilename,atomname,atomtype)
    atomeach=getatomeach(hmmfilename,moleculetemp2filename,atomfilename,N,step)
    allmoleculeroute=printatomroute(atomroutefilename,N,step,atomeach,atomtype,atomname,mname,timestep)
    return allmoleculeroute,mname

def step4(allmoleculeroute,mname,reactionfilename,tablefilename):
    allroute=getallroute(reactionfilename,allmoleculeroute,mname)
    printtable(tablefilename,reactionfilename,allroute)
######## run ########
def run(bondfilename="bonds.reaxc",atomname=["C","H","O"],originfilename="originsignal.txt",hmmfilename="hmmsignal.txt",atomfilename="atom.txt",moleculefilename="moleculename.txt",atomroutefilename="atomroute.txt",reactionfilename="reaction.txt",tablefilename="table.txt",moleculetempfilename="moleculetemp.txt",moleculetemp2filename="moleculetemp2.txt",moleculestructurefilename="moleculestructure.txt",stepinterval=1,states=[0,1],observations=[0,1],p=[0.5,0.5],a=[[0.999,0.001],[0.001,0.999]],b=[[0.6, 0.4],[0.4, 0.6]],runHMM=True,getoriginfile=False):
    ######### Parameter above ############
    ######start#####
    print("Run HMM calculation:")
    timearray=printtime([])
    for runstep in range(1,5):
        ######## step 1 ##### 
        if(runstep==1):
            N,atomtype,step,timestep=step1(bondfilename,moleculetempfilename,stepinterval)
        ######## step 2 ##### 
        elif(runstep==2):
            step2(states,observations,p,a,b,originfilename,hmmfilename,moleculetempfilename,moleculetemp2filename,step,runHMM,getoriginfile)
        ######## step 3 ##### 
        elif(runstep==3):
            allmoleculeroute,mname=step3(atomname,atomtype,N,step,timestep,moleculefilename,hmmfilename if runHMM else originfilename,atomfilename,moleculetemp2filename,atomroutefilename,moleculestructurefilename)
        ######## step 4 ##### 
        elif(runstep==4):
            step4(allmoleculeroute,mname,reactionfilename,tablefilename)
        #garbage collect
        gc.collect()
        timearray=printtime(timearray)
    ####### end #######
    print()
    print("Time consumed:")
    for i in range(1,len(timearray)):
        print("Step ",i," consumed: ",round(timearray[i]-timearray[i-1],3),"s")
    print("Total time:",round(timearray[-1]-timearray[0],3),"s")
    print()

#####draw#####    
def draw(tablefilename="table.txt",imagefilename="image.svg",moleculestructurefilename="moleculestructure.txt",species={},node_size=200,font_size=6,widthcoefficient=3,show=False,maxspecies=20,n_color=256,atomname=["C","H","O"]):
    #start
    print("Draw the image:")
    timearray=printtime([])
    #read table
    table,name=readtable(tablefilename)
    species,showname=handlespecies(species,name,maxspecies,atomname,moleculestructurefilename)

    #make color
    start = np.array([0, 1, 0])
    end = np.array([0, 0, 1])
    colorsRGB=[(start + i*(end-start) / n_color) for i in range(n_color)]
          
    G = nx.DiGraph()
    for i in range(len(table)):
        if (name[i] in species):
            G.add_node(showname[name[i]] if name[i] in showname else name[i])
        for j in range(len(table)):
            if (name[i] in species and name[j] in species):
                if(table[i][j]>0):
                    G.add_weighted_edges_from([((showname[name[i]] if name[i] in showname else name[i]),(showname[name[j]] if name[j] in showname else name[j]),table[i][j])])
    weights = [G[u][v]['weight'] for u,v in G.edges()]
    widths=[weight/max(weights) *widthcoefficient for weight in weights]
    colors=[colorsRGB[math.floor(width/max(widths)*(n_color-1))] for width in widths]
    try:
        nx.draw(G,pos = nx.spring_layout(G),width=widths,node_size=node_size,font_size=font_size,with_labels=True,edge_color=colors)            
        plt.savefig(imagefilename)
        if show:
            plt.show()
    except:
        print("Error: cannot draw images")

    timearray=printtime(timearray)
    ####### end #######
    print()
    print("Time consumed:")
    for i in range(1,len(timearray)):
        print("Step ",i," consumed: ",round(timearray[i]-timearray[i-1],3),"s")
    print("Total time:",round(timearray[-1]-timearray[0],3),"s")
    print()
#### run and draw ####
def runanddraw(bondfilename="bonds.reaxc",atomname=["C","H","O"],originfilename="originsignal.txt",hmmfilename="hmmsignal.txt",atomfilename="atom.txt",moleculefilename="moleculename.txt",atomroutefilename="atomroute.txt",reactionfilename="reaction.txt",tablefilename="table.txt",moleculetempfilename="moleculetemp.txt",moleculetemp2filename="moleculetemp2.txt",moleculestructurefilename="moleculestructure.txt",imagefilename="image.svg",stepinterval=1,states=[0,1],observations=[0,1],p=[0.5,0.5],a=[[0.999,0.001],[0.001,0.999]],b=[[0.6, 0.4],[0.4, 0.6]],runHMM=True,getoriginfile=False,species={},node_size=200,font_size=6,widthcoefficient=3,show=False,maxspecies=20,n_color=256):
    run(bondfilename,atomname,originfilename,hmmfilename,atomfilename,moleculefilename,atomroutefilename,reactionfilename,tablefilename,moleculetempfilename,moleculetemp2filename,moleculestructurefilename,stepinterval,states,observations,p,a,b,runHMM,getoriginfile)
    draw(tablefilename,imagefilename,moleculestructurefilename,species,node_size,font_size,widthcoefficient,show,maxspecies,n_color,atomname)
    
##### main #####
if __name__ == '__main__':
    runanddraw()
