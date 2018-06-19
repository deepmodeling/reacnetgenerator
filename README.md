# Reaction Network Generator (ReacNetGenerator)
[![Anaconda-Server Badge](https://anaconda.org/njzjz/reacnetgenerator/badges/installer/conda.svg)](https://conda.anaconda.org/njzjz)

An automatic generator of reaction network for reactive molecular dynamics simulation.
## Features
- Processing of MD trajectory containing atomic coordinates or bond orders
- Hidden Markov Model (HMM) based noise filtering
- Isomers identifying accoarding to SMILES
- Generation of reaction network for visualization using force-directed algorithm
- Parallel computing
## Installation
### Linux 64 Under anaconda python (fastest install)
Before installation you need to [get conda](https://conda.io/docs/user-guide/install/index.html) first.
```sh
$ conda install -c njzjz -c openbabel -c rdkit -c omnia reacnetgenerator
```
### With pip
If you install ReacNetGenerator with pip, you must install [RDKit](https://github.com/rdkit/rdkit) and [OpenBabel](https://github.com/openbabel/openbabel) on your own.
```sh
$ pip install reacnetgenerator
```
## Simple example
Process a LAMMPS bond file named bonds.reaxc.
```python
>>> import ReacNetGenerator
>>> ReacNetGenerator.runanddraw(inputfiletype="lammpsbondfile",inputfilename="bonds.reaxc",atomname=["C","H","O"])
```
## Reference
### ReacNetGenerator.runanddraw
```python
runanddraw(inputfiletype='lammpsbondfile',inputfilename='bonds.reaxc',atomname=['C','H','O'],originfilename='originsignal.txt',hmmfilename='hmmsignal.txt',atomfilename='atom.txt',moleculefilename='moleculename.txt',atomroutefilename='atomroute.txt',reactionfilename='reaction.txt',tablefilename='table.txt',moleculetempfilename='moleculetemp.txt',moleculetemp2filename='moleculetemp2.txt',moleculestructurefilename='moleculestructure.txt',imagefilename='image.svg',stepinterval=1,p=[0.5,0.5],a=[[0.999,0.001],[0.001,0.999]],b=[[0.6, 0.4],[0.4, 0.6]],runHMM=True,SMILES=True,getoriginfile=False,printfiltersignal=False,species={},node_size=200,font_size=6,widthcoefficient=1,show=False,maxspecies=20,n_color=256,nolabel=False,filter=[],node_color=[135/256,206/256,250/256],pos={},showid=True,k=None)
```
**Parameters:**
- **inputfiletype**(_string (optional, default='lammpsbondfile')_)- Type of the input file. 'lammpsbondfile' or 'lammpsdumpfile' is acceptable.
- **inputfilename**(_string (optional, default='bonds.reaxc')_)- Path of the input file.
- **atomname**(_list (optional, default=['C','H','O'])_)- Names of atom types, corresponding atom type 1, atom type 2 and so on.
- **originfilename**(_string (optional, default='originsignal.txt')_)- Path of the original signal file.
- **hmmfilename**(_string (optional, default='hmmsignal.txt')_)- Path of the HMM signal file.
- **atomfilename**(_string (optional, default='atom.txt')_)- Path of the atom file, containing molecule IDs of atoms.
- **moleculefilename**(_string (optional, default='moleculename.txt')_)- Path of the molecule name and structure file.
- **atomroutefilename**(_string (optional, default='atomroute.txt')_)- Path of the atom route file.
- **reactionfilename**(_string (optional, default='reaction.txt')_)- Path of the reaction file, containing reactions and numbers of them.
- **tablefilename**(_string (optional, default='table.txt')_)- Path of the reaction matrix file.
- **moleculetempfilename**(_string (optional, default='moleculetemp.txt')_)- Path of the temp file.
- **moleculetemp2filename**(_string (optional, default='moleculetemp2.txt')_)- Path of another temp file.
- **moleculestructurefilename**(_string (optional, default='moleculestructure.txt')_)- Path of the species structure file if `SMILES=False`.
- **imagefilename**(_string (optional, default='image.svg')_)- Path of the reaction network image file.
- **stepinterval**(_int (optional, default=1)_)- Processing a timestep every N timesteps.
- **p**(_list (optional, default=[0.5,0.5])_)- The initial state vector of HMM.
- **a**(_list (optional, default=[[0.999,0.001],[0.001,0.999]])_)- The transition matrix of HMM.
- **b**(_list (optional, default=[[0.6, 0.4],[0.4, 0.6]])_)- The emission matrix of HMM.
- **runHMM**(_boolen (optional, default=True)_)- If True use HMM for noise filtering.
- **SMILES**(_boolen (optional, default=True)_)- If True use SMILES to identify isomers. Otherwise use VF2 algorithm.
- **getoriginfile**(_boolen (optional, default=False)_)- If True print original signals.
- **printfiltersignal**(_boolen (optional, default=False)_)- If True print filtered signals.
- **species**(_dictionary (optional, default={})_)- Species shown in the reaction network. If a dictionary is empty, several of the most species will be taken.
- **node_size**(_int (optional, default=200)_)- Size of nodes in the reaction network.
- **font_size**(_int (optional, default=6)_)- Size of fonts in the reaction network.
- **widthcoefficient**(_int (optional, default=1)_)- Coefficient of widths in the reaction network.
- **show**(_boolen (optional, default=False)_)- If True show the reaction network while running the script.
- **maxspecies**(_int (optional, default=20)_)- Max number of species in the reaction network if `species={}`.
- **n_color**(_int (optional, default=256)_)- Number of arrow colors in the reaction network.
- **nolabel**(_boolen (optional, default=False)_)- If True generate a reaction network without node labels.
- **filter**(_list (optional, default=[])_)- Filtered species in the reaction network.
- **node_color**(_list (optional, default=[135/256,206/256,250/256])_)- The RGB color in the reaction network.
- **pos**(_dictionary (optional, default={})_)- Positions of nodes in the reaction network.
- **showid**(_string (_boolen, default=True)_)- If True label species ID in nodes. Otherwise label names of them.
- **k**(_int or None (optional, defaultNone)_)- k for position of nodes.
