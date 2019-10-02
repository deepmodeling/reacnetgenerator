# ReacNetGenerator

[![DOI:10.26434/chemrxiv.7421534](https://zenodo.org/badge/DOI/10.26434/chemrxiv.7421534.svg)](https://doi.org/10.26434/chemrxiv.7421534)
[![Research Group](https://img.shields.io/website-up-down-green-red/http/computchem.cn.svg?label=Research%20Group)](http://computchem.cn) [![Greenkeeper badge](https://badges.greenkeeper.io/njzjz/reacnetgenerator.svg?token=f40c024fe456fbe93b531cfb0fc30315211fd10d04eadb05f9f969da0a893a03&ts=1569994170255)](https://greenkeeper.io/)

An automatic generator of reaction network for reactive molecular dynamics simulation.

ReacNetGenerator: an Automatic Reaction Network Generator for Reactive Molecular Dynamic Simulations, 2019, doi: [10.26434/chemrxiv.7421534](https://dx.doi.org/10.26434/chemrxiv.7421534)

jinzhe.zeng@rutgers.edu (Jinzhe Zeng), tzhu@lps.ecnu.edu.cn (Tong Zhu)

## Features

-   Processing of MD trajectory containing atomic coordinates or bond orders
-   Hidden Markov Model (HMM) based noise filtering
-   Isomers identifying accoarding to SMILES
-   Generation of reaction network for visualization using force-directed algorithm
-   Parallel computing

## Installation

First, you need to download the source code on [our group website](http://computchem.cn/reacnetgenerator/) or email us to get the newest one. Then install ReacNetGenerator with one of the following guides:

### Building a conda package
1. [Install Anaconda or Miniconda](https://conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html) to obtain conda.
2. Decompress reacnetgenerator.zip and build in the main directory of ReacNetGenerator:

```bash
conda config --add channels conda-forge
conda build conda/recipe
conda install reacnetgenerator --use-local
reacnetgenerator -h
```

### Building a Docker Image
1. [Install Docker](https://docs.docker.com/install/).
2. Decompress reacnetgenerator.zip and build in the main directory of ReacNetGenerator:

```bash
docker build . -t njzjz/reacnetgenerator
docker run njzjz/reacnetgenerator reacnetgenerator -h
```

## Usage

### Command line

ReacNetGenerator can process any kind of trajectory files containing atomic coordinates, e.g. a LAMMPS dump file prepared by running “dump 1 all custom 100 dump.reaxc id type x y z” in LAMMPS:

```bash
reacnetgenerator --dump -i dump.reaxc -a C H O
```
where C, H, and O are atomic names in the input file. [Analysis report](https://njzjz.github.io/reacnetgenerator/report.html?jdata=https%3A%2F%2Fgist.githubusercontent.com%2Fnjzjz%2Fe9a4b42ceb7d2c3c7ada189f38708bf3%2Fraw%2F83d01b9ab1780b0ad2d1e7f934e61fa113cb0f9f%2Fmethane.json) will be generated automatically.

Also, ReacNetGenerator can process files containing bond information, e.g. LAMMPS bond file:

```bash
reacnetgenerator -i bonds.reaxc -a C H O
```

You can running the following script for help:

```bash
reacnetgenerator -h
```

### GUI version

You can open a GUI version for ReacNetGenerator by typing:

```bash
reacnetgeneratorgui
```

## Awards
* The First Prize in 2019 (the 11th Session) Shanghai Computer Application Competition for College Students
* The First Prize in 2019 (the 12th Session) Chinese Computer Design Competition for College Students
