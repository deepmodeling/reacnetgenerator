# ReacNetGenerator
![python3.7](https://img.shields.io/badge/python-3.7-blue.svg)

An automatic generator of reaction network for reactive molecular dynamics simulation.

**Please cite:** J. Zeng, L. Cao, J.Z.H. Zhang, C.H. Chin, T. Zhu: ReacNetGen: an Automatic Reaction Network Generator for Reactive Molecular Dynamic Simulations, 2018, doi: [10.26434/chemrxiv.7421534](https://dx.doi.org/10.26434/chemrxiv.7421534)

**Author**: [Jinzhe Zeng](https://cv.njzjz.win), Liqun Cao, [John ZH Zhang](https://research.shanghai.nyu.edu/centers-and-institutes/chemistry/people/john-zenghui-zhang), Chih-Hao Chin, [Tong Zhu](http://computchem.cn/people/)

**Email**: tzhu@lps.ecnu.edu.cn, jzzeng@stu.ecnu.edu.cn

[Research Group](http://computchem.cn/)

## Features
- Processing of MD trajectory containing atomic coordinates or bond orders
- Hidden Markov Model (HMM) based noise filtering
- Isomers identifying accoarding to SMILES
- Generation of reaction network for visualization using force-directed algorithm
- Parallel computing

## Requirements
* Python >= 3.6 (**Note:** Python 2 is not supported!)
* Python packages: [numpy](https://github.com/numpy/numpy), [scipy](https://github.com/scipy/scipy), [networkx](https://github.com/networkx/networkx), [scikit-learn](https://github.com/scikit-learn/scikit-learn), [matplotlib](https://github.com/matplotlib/matplotlib), [hmmlearn](https://github.com/hmmlearn/hmmlearn), [htmlmin](https://github.com/mankyd/htmlmin/), [ASE](https://gitlab.com/ase/ase), [scour](https://github.com/scour-project/scour)
* Extra packages: [OpenBabel](https://github.com/openbabel/openbabel), [RDKit](https://github.com/rdkit/rdkit)

## Installation
1. [Get conda](https://conda.io/docs/user-guide/install/index.html) to install Python 3.

2. Use conda to install extra packages:
```sh
conda install -c openbabel openbabel 
conda install -c rdkit rdkit
```

3. Use pip to install required packages: 
```sh
$ pip install numpy scipy networkx scikit-learn matplotlib hmmlearn htmlmin ase scour
```

4. Download ReacNetGenerator and build it from source:
```sh
$ cd ReacNetGenerator/
$ python3 setup.py install
```

## Simple example
Prepare a [LAMMPS bond file](http://lammps.sandia.gov/doc/fix_reax_bonds.html) named bonds.reaxc, then run the script:

```sh
$ reacnetgenerator -i bonds.reaxc -a C H O
```

where C, H, and O are atomic names in the input file. [Analysis report](https://njzjz.github.io/reacnetgenerator/report.html) will be generated automatically.  

A [LAMMPS dump file](https://lammps.sandia.gov/doc/dump.html) is also supported. You can prepare it by running "dump 1 all custom 100 dump.reaxc id type x y z" in LAMMPS.

```sh
$ reacnetgenerator --dump -i dump.reaxc -a C H O
```

You can running the following script for help:

```sh
$ reacnetgenerator -h
```

## GUI version
You can open a GUI version for ReacNetGenerator by typing:

```sh
$ reacnetgeneratorgui
```
