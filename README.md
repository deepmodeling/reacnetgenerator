# Reaction Network Generator (ReacNetGenerator)
[![python3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://badge.fury.io/py/ReacNetGenerator)[![pypi](https://badge.fury.io/py/ReacNetGenerator.svg)](https://badge.fury.io/py/ReacNetGenerator)[![Anaconda-Server Badge](https://anaconda.org/njzjz/reacnetgenerator/badges/version.svg)](https://conda.anaconda.org/njzjz/reacnetgenerator)

An automatic generator of reaction network for reactive molecular dynamics simulation.

**Author**: Jinzhe Zeng

**Email**: njzjz@qq.com 10154601140@stu.ecnu.edu.cn

[Research Group](http://computchem.cn/)
## Features
- Processing of MD trajectory containing atomic coordinates or bond orders
- Hidden Markov Model (HMM) based noise filtering
- Isomers identifying accoarding to SMILES
- Generation of reaction network for visualization using force-directed algorithm
- Parallel computing

## Requirements
* Python 3 (Note: Python 2 is not supported!)
* [numpy](https://github.com/numpy/numpy)
* [scipy](https://github.com/scipy/scipy)
* [networkx](https://github.com/networkx/networkx)
* [scikit-learn](https://github.com/scikit-learn/scikit-learn)
* [matplotlib](https://github.com/matplotlib/matplotlib)
* [hmmlearn](https://github.com/hmmlearn/hmmlearn)
* [OpenBabel](https://github.com/openbabel/openbabel)
* [RDKit](https://github.com/rdkit/rdkit)

## Installation
### Linux 64 under anaconda python (fastest install)
Before installation you need to [get conda](https://conda.io/docs/user-guide/install/index.html) first.
```sh
$ conda install -c njzjz -c bioconda reacnetgenerator
```
### With pip
If you install ReacNetGenerator with pip, you must install [RDKit](https://github.com/rdkit/rdkit) and [OpenBabel](https://github.com/openbabel/openbabel) on your own.
```sh
$ pip install reacnetgenerator
```
### Build from source
If you build ReacNetGenerator from source, you must install [RDKit](https://github.com/rdkit/rdkit) and [OpenBabel](https://github.com/openbabel/openbabel) on your own.
```sh
$ git clone https://github.com/njzjz/ReacNetGenerator.git
$ cd ReacNetGenerator/
$ python3 setup.py install
```
## Simple example
Process a [LAMMPS bond file](http://lammps.sandia.gov/doc/fix_reax_bonds.html) named bonds.reaxc.
```python
>>> from ReacNetGenerator import ReacNetGenerator
>>> ReacNetGenerator(inputfiletype="lammpsbondfile",inputfilename="bonds.reaxc",atomname=["C","H","O"]).runanddraw()
```
## Reference
### ReacNetGenerator.ReacNetGenerator
```python
class ReacNetGenerator(inputfiletype='lammpsbondfile',inputfilename='bonds.reaxc',atomname=['C','H','O'],originfilename=None,hmmfilename=None,atomfilename=None,moleculefilename=None,atomroutefilename=None,reactionfilename=None,tablefilename=None,moleculetempfilename=None,moleculetemp2filename=None,moleculestructurefilename=None,imagefilename=None,stepinterval=1,p=[0.5,0.5],a=[[0.999,0.001],[0.001,0.999]],b=[[0.6, 0.4],[0.4, 0.6]],runHMM=True,SMILES=True,getoriginfile=False,printfiltersignal=False,species={},node_size=200,font_size=6,widthcoefficient=1,show=False,maxspecies=20,n_color=256,nolabel=False,filter=[],node_color=[78/256,196/256,238/256],pos={},showid=True,k=None)
```
**Parameters:**
- **inputfiletype**(_string (optional, default='lammpsbondfile')_)- Type of the input file. 'lammpsbondfile' or 'lammpsdumpfile' is acceptable.
- **inputfilename**(_string (optional, default='bonds.reaxc')_)- Path of the input file.
- **atomname**(_list (optional, default=['C','H','O'])_)- Names of atom types, corresponding atom type 1, atom type 2 and so on.
- **originfilename**(_string (optional, default=None)_)- Path of the original signal file. If None it will be set to inputfilename+'.origin'.
- **hmmfilename**(_string (optional, default=None)_)- Path of the HMM signal file. If None it will be set to inputfilename+'.hmm'.
- **atomfilename**(_string (optional, default=None)_)- Path of the atom file, containing molecule IDs of atoms. If None it will be set to inputfilename+'.atom'.
- **moleculefilename**(_string (optional, default=None)_)- Path of the molecule name and structure file. If None it will be set to inputfilename+'.moname'.
- **atomroutefilename**(_string (optional, default=None)_)- Path of the atom route file. If None it will be set to inputfilename+'.route'.
- **reactionfilename**(_string (optional, default=None)_)- Path of the reaction file, containing reactions and numbers of them. If None it will be set to inputfilename+'.reaction'.
- **tablefilename**(_string (optional, default=None)_)- Path of the reaction matrix file. If None it will be set to inputfilename+'.table'.
- **moleculetempfilename**(_string (optional, default=None)_)- Path of the temp file. If None it will be set to inputfilename+'.temp'.
- **moleculetemp2filename**(_string (optional, default=None)_)- Path of another temp file. If None it will be set to inputfilename+'.temp2'.
- **moleculestructurefilename**(_string (optional, default=None)_)- Path of the species structure file if `SMILES=False`. If None it will be set to inputfilename+'.structure'.
- **imagefilename**(_string (optional, default=None)_)- Path of the reaction network image file. If None it will be set to inputfilename+'.svg'.
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
- **node_color**(_list (optional, default=[78/256,196/256,238/256])_)- The RGB color in the reaction network.
- **pos**(_dictionary (optional, default={})_)- Positions of nodes in the reaction network.
- **showid**(_string (_boolen, default=True)_)- If True label species ID in nodes. Otherwise label names of them.
- **k**(_int or None (optional, defaultNone)_)- k for position of nodes.

**Methods:**
- **ReacNetGenerator.runanddraw**() - Process the input file and draw the network.
- **ReacNetGenerator.run**() - Only process the input file.
- **ReacNetGenerator.draw**() - Only draw the network.
