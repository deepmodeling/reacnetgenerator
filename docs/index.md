# Reaction Network Generator (ReacNetGenerator)
An automatic generation of reaction network for reactive molecular dynamics simulation.
## Features
- Processing of MD trajectory containing atomic coordinates or bond orders
- Hidden Markov Model (HMM) based noise filtering
- Isomers identifying accoarding to SMILES
- Generation of reaction network for visualization using force-directed algorithm
## Environment
Please install Anaconda3 first.
```
$ conda install openbabel -c openbabel
$ conda install rdkit -c rdkit
$ pip install numpy scipy networkx sklearn matplotlib hmmlearn
```
## Simple example
Process a LAMMPS bond file named bonds.reaxc.
```
>>> import ReacNetGenerator
>>> ReacNetGenerator.runanddraw(inputfiletype="lammpsbondfile",inputfilename="bonds.reaxc",atomname=["C","H","O"])
```
