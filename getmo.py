# -*- coding: UTF-8 -*-  
#updated at 2018/4/6 23:00

######### Parameter ############
bondfilename="bonds.reaxc"#filename
atomname=["C","H","O"]#atom type
originfilename="originsignal.txt"
hmmfilename="hmmsignal.txt"
atomfilename="atom.txt"
moleculefilename="moleculename.txt"
atomroutefilename="atomroute.txt"
reactionfilename="reaction.txt"
tablefilename="table.txt"
moleculetempfilename="moleculetemp.txt"
moleculetemp2filename="moleculetemp2.txt"
stepinterval=1
######## HMM Parameter ########
states = [0,1]
observations = [0,1]
p=[0.5,0.5]
a=[[0.999,0.001],
    [0.001,0.999]]
b=[[0.6, 0.4],
    [0.4, 0.6]]
######### import #########
import time
import gc
from multiprocessing import Pool, Semaphore
import numpy as np
from hmmlearn import hmm
from functools import reduce
import networkx as nx
import networkx.algorithms.isomorphism as iso
######## function ########
def printtime(timearray):
    timearray.append(time.time())
    if len(timearray)>1:
       print("Step ",len(timearray)-1," has been completed. Time consumed: ",round(timearray[-1]-timearray[-2],3),"s")
    return timearray

def cal(func,data):
    with Pool() as pool:
        return pool.map(func,data,100)

def readatomnumber(bondfilename):
    with open(bondfilename,buffering=4096) as file: #read Atom number
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
        if i<b:
            bo=(i,b,l)
        else:
            bo=(b,i,l)
        if not bo in bondlist:
            bondlist.append(bo)
        if not done[b]:
            molecule,done,bondlist=mo(b,bond,level,molecule,done,bondlist)
    return molecule,done,bondlist

def openbondfile(bondfilename,moleculetempfilename,N):
    step=0
    realstep=0
    iscontent=False
    isfirststep=True
    timestep=[]
    atomtype=[0 for x in range(N+1)]
    d={}
    with open(bondfilename,buffering=4096) as file:
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
    with open(moleculetempfilename,'wt',buffering=4096) as f:
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

def gethmm(ori,model):
    o = np.array([ori]).T
    logprob, h = model.decode(o, algorithm="viterbi")
    hmmlist=list(map(lambda x: states[x], h))
    return hmmlist

def produce(semaphore, from_file):
    with open(from_file) as reader:
        for line in reader:
            # Reduce Semaphore by 1 or wait if 0
            semaphore.acquire()
            # Now deliver an item to the caller (pool)
            yield line

def getoriginandhmm(item):
    list=item.split()
    value=[int(x) for x in list[-1].split(",")]
    origin= [1 if i in value else 0 for i in range(1, step+1)]
    global model
    hmm=gethmm(origin,model)
    return origin,hmm,item

def calhmm(originfilename,hmmfilename,moleculetempfilename,moleculetemp2filename):
    with open(originfilename, 'wt',buffering=4096) as fo,open(hmmfilename, 'wt',buffering=4096) as fh,open(moleculetemp2filename,'wt',buffering=4096) as ft,Pool(maxtasksperchild=100) as pool:
        semaphore = Semaphore(512)
        results=pool.imap_unordered(getoriginandhmm,produce(semaphore,moleculetempfilename),10)
        for originsignal,hmmsignal,mlist in results:
            print("".join([str(i) for i in originsignal]), file=fo)
            print("".join([str(i) for i in hmmsignal]), file=fh)
            print(mlist,end='',file=ft)
            semaphore.release()

def getatomroute(i):
    global step,atomeach,atomtype,atomname
    route=[]
    routestrarr=[]
    molecule=0
    for j in range(0,step):
        if atomeach[i][j]>0 and atomeach[i][j]!=molecule:
            routestrarr.append(mname[atomeach[i][j]-1] + " ("+ str(atomeach[i][j])+" step "+str(timestep[j])+")")
            molecule=atomeach[i][j]
            route.append(molecule)
    routestr="Atom "+str(i)+" "+atomname[atomtype[i]-1]+": "+" -> ".join(routestrarr)
    return route,routestr

def getallmoleculeroute(i):
    global atomroute
    moleculeroute=[]
    for j in range(len(atomroute[i])-1):
        left=[atomroute[i][j]]
        right=[atomroute[i][j+1]]                                                     
        if not [set(left),set(right)] in moleculeroute:
            moleculeroute.append([set(left),set(right)])
    return moleculeroute

def printatomroute(atomroutefilename,N):
    with open(atomroutefilename, 'wt',buffering=4096) as f,Pool() as pool:
        atomroute=[]
        for route in pool.imap(getatomroute,range(1,N+1)):
            routearray,routestr=route
            print(routestr, file=f)
            atomroute.append(routearray)
    return atomroute

def printmoleculename(moleculefilename,moleculetempfilename,atomname,atomtype):
    mname=[]
    d={}
    em = iso.numerical_edge_match(['atom','level'], [-1,1])
    with open(moleculefilename, 'wt',buffering=4096) as fm,open(moleculetempfilename,buffering=4096) as ft:
        for line in ft:
            list=line.split()
            key1=[int(x) for x in list[0].split(",")]
            if len(list)==3:
                key2=[tuple(int(y) for y in x.split(",")) for x in list[1].split(";")]
            else:
                key2=[]
            typenumber=[0 for i in range(len(atomname))]
            G=nx.Graph()
            for line in key2:
                G.add_edge(line[0],line[1],level=line[2])
            for atomnumber in key1:
                G.add_node(atomnumber, atom=atomtype[atomnumber])
                typenumber[atomtype[atomnumber]-1]+=1
            name="".join([atomname[i]+(str(typenumber[i] if typenumber[i]>1 else "")) if typenumber[i]>0 else "" for i in range(0,len(atomname))])                
            if name in d:
                for j in range(len(d[name])):
                    if nx.is_isomorphic(G,d[name][j],em):
                        if(j>0):
                            name+="_"+str(j+1)
                        break
                else:
                    d[name].append(G)
                    name+="_"+str(len(d[name]))
            else:
                d[name]=[G]
            mname.append(name)
            print(name,key1,key2, file=fm)
    return mname

def getatomeach(hmmfilename,moleculetemp2filename,atomfilename):
    atomeach=[[0 for j in range(0,step)] for i in range(0,N+1)]
    with open(hmmfilename,buffering=4096) as fh,open(moleculetemp2filename,buffering=4096) as ft:
        i=0
        for lineh,linet in zip(fh,ft):
            list=linet.split()
            key1=[int(x) for x in list[0].split(",")]
            for j in range(0,step):#j is step
                if lineh[j]=="1":
                    for a in key1:
                        atomeach[a][j]=i+1
            i+=1
    with open(atomfilename, 'wt',buffering=4096) as f:
        for atom in atomeach[1:]:
            print(atom, file=f)
    return atomeach

def getallroute(reactionfilename,allmoleculeroute,mname):
    allroute={}
    for moleculeroute in allmoleculeroute:
        left=right=""
        leftlist=[mname[m-1] for m in moleculeroute[0]]
        rightlist=[mname[m-1] for m in moleculeroute[1]]
        leftlist.sort()
        rightlist.sort()
        for l in leftlist[:]:
            if l in rightlist:
                rightlist.remove(l)
                leftlist.remove(l)
        if not (leftlist and rightlist):
            continue
        equation="+".join(leftlist)+"->"+"+".join(rightlist)
        if equation in allroute:
            allroute[equation]+=1
        else:
            allroute[equation]=1            
    return allroute

def printtable(tablefilename,allroute):
    species=[]
    table=[[0 for x in range(100)] for x in range(100)]
    with open(reactionfilename,'wt',buffering=4096) as f:
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
    with open(tablefilename,'wt',buffering=4096) as f:
        print("\t"+"\t".join(species),file=f)
        for i in range(len(species)):
            print(species[i],end='\t',file=f)
            for j in range(len(species)):
                print(table[i][j],end='\t',file=f)
            print(file=f)
  

######## run ########
if __name__ == '__main__':
    timearray=printtime([])
    ######## step 1 ##### 
    N=readatomnumber(bondfilename)
    atomtype,step,timestep=openbondfile(bondfilename,moleculetempfilename,N)
    gc.collect()
    timearray=printtime(timearray)
    ######## step 2 ##### 
    model=initHMM(states,observations,p,a,b)
    gc.collect()
    calhmm(originfilename,hmmfilename,moleculetempfilename,moleculetemp2filename)
    del model
    gc.collect()
    timearray=printtime(timearray)
    ######## step 3 ##### 
    mname=printmoleculename(moleculefilename,moleculetemp2filename,atomname,atomtype)
    gc.collect()
    atomeach=getatomeach(hmmfilename,moleculetemp2filename,atomfilename)
    gc.collect()
    atomroute=printatomroute(atomroutefilename,N)
    del atomeach,atomtype,timestep
    gc.collect()
    allmoleculeroute=reduce(lambda x,y:x if y in x else x + [y], [[], ] +sum(cal(getallmoleculeroute,range(len(atomroute))),[]))
    del atomroute
    gc.collect()
    allroute=getallroute(reactionfilename,allmoleculeroute,mname)
    del allmoleculeroute
    gc.collect()
    printtable(tablefilename,allroute)
    del allroute
    gc.collect()
    timearray=printtime(timearray)
    ####### end #######
    print()
    print("Time consumed:")
    for i in range(1,len(timearray)):
        print("Step ",i," consumed: ",round(timearray[i]-timearray[i-1],3),"s")
    print("Total time:",round(timearray[-1]-timearray[0],3),"s")
