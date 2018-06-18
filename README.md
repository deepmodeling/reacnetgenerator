# Reaction Network Generator (ReacNetGenerator)
[![Anaconda-Server Badge](https://anaconda.org/njzjz/reacnetgenerator/badges/installer/conda.svg)](https://conda.anaconda.org/njzjz)

An automatic generation of reaction network for reactive molecular dynamics simulation.
## Features
- Processing of MD trajectory containing atomic coordinates or bond orders
- Hidden Markov Model (HMM) based noise filtering
- Isomers identifying accoarding to SMILES
- Generation of reaction network for visualization using force-directed algorithm
## Installation
### Linux 64 Under anaconda python (fastest install)
Before installation you need to [get conda](https://conda.io/docs/user-guide/install/index.html) first.
```
$ conda install -c njzjz -c openbabel -c rdkit -c omnia reacnetgenerator
```
### With pip
If you install ReacNetGenerator with pip, you must install [RDKit](https://github.com/rdkit/rdkit) and [OpenBabel](https://github.com/openbabel/openbabel) on your own.
```
$ pip install reacnetgenerator
```
## Simple example
Process a LAMMPS bond file named bonds.reaxc.
```
>>> import ReacNetGenerator
>>> ReacNetGenerator.runanddraw(inputfiletype="lammpsbondfile",inputfilename="bonds.reaxc",atomname=["C","H","O"])
```
