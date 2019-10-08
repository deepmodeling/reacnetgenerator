---
home: true
heroImage: /reacnetgen.svg
heroText: ReacNetGenerator
tagline: An automatic generator of reaction network for reactive molecular dynamics simulation
actionText: Download
actionLink: https://computchem.cn/reacnetgenerator/
features:
- title: Trajectory
  details: Processing of MD trajectory containing atomic coordinates or bond orders
- title: Filtering
  details: Hidden Markov Model (HMM) based noise filtering
- title: Isomers
  details: Isomers identifying accoarding to SMILES
- title: Network
  details: Generation of reaction network for visualization using force-directed algorithm
- title: HTML5
  details: Showing an interactive web page
- title: Fast
  details: Parallel computing and performance optimization
  footer: Copyright © 2018-2019 East China Normal University
---

# Citation and contact

ReacNetGenerator: an Automatic Reaction Network Generator for Reactive Molecular Dynamic Simulations, 2019, doi: [10.26434/chemrxiv.7421534](https://dx.doi.org/10.26434/chemrxiv.7421534)

jinzhe.zeng@rutgers.edu (Jinzhe Zeng), tzhu@lps.ecnu.edu.cn (Tong Zhu)

# Installation

First, you need to download the source code on [our group website](https://computchem.cn/reacnetgenerator/) or email us to get the newest one. Then install ReacNetGenerator with one of the following guides:

## Building a conda package
1. [Install Anaconda or Miniconda](https://conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html) to obtain conda.
2. Decompress reacnetgenerator.zip and build in the main directory of ReacNetGenerator:

```bash
conda config --add channels conda-forge
conda build conda/recipe
conda install reacnetgenerator --use-local
reacnetgenerator -h
```

## Building a Docker Image
1. [Install Docker](https://docs.docker.com/install/).
2. Decompress reacnetgenerator.zip and build in the main directory of ReacNetGenerator:

```bash
docker build . -t njzjz/reacnetgenerator
docker run njzjz/reacnetgenerator reacnetgenerator -h
```

# Usage

## Command line

ReacNetGenerator can process any kind of trajectory files containing atomic coordinates, e.g. a LAMMPS dump file prepared by running “dump 1 all custom 100 dump.reaxc id type x y z” in LAMMPS:

```bash
reacnetgenerator --dump -i dump.reaxc -a C H O
```
where C, H, and O are atomic names in the input file. <a href="/report.html?jdata=https%3A%2F%2Fgist.githubusercontent.com%2Fnjzjz%2Fe9a4b42ceb7d2c3c7ada189f38708bf3%2Fraw%2F83d01b9ab1780b0ad2d1e7f934e61fa113cb0f9f%2Fmethane.json" target="_blank">Analysis report</a> will be generated automatically.

Also, ReacNetGenerator can process files containing bond information, e.g. LAMMPS bond file:

```bash
reacnetgenerator -i bonds.reaxc -a C H O
```

You can running the following script for help:

```bash
reacnetgenerator -h
```

## GUI version

You can open a GUI version for ReacNetGenerator by typing:

```bash
reacnetgeneratorgui
```

# Awards
* The First Prize in 2019 (the 11th Session) Shanghai Computer Application Competition for College Students
* The First Prize in 2019 (the 12th Session) Chinese Computer Design Competition for College Students
