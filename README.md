# <img src=docs/.vuepress/public/reacnetgen.svg height=40/>  ReacNetGenerator

[![DOI:10.1039/C9CP05091D](https://zenodo.org/badge/DOI/10.1039/C9CP05091D.svg)](https://doi.org/10.1039/C9CP05091D)
[![Research Group](https://img.shields.io/website-up-down-green-red/https/computchem.cn.svg?label=Research%20Group)](https://computchem.cn)

An automatic reaction network generator for reactive molecular dynamics simulation.

ReacNetGenerator: an automatic reaction network generator for reactive molecular dynamic simulations, Phys. Chem. Chem. Phys., 2020, 22 (2): 683–691, doi: [10.1039/C9CP05091D](https://dx.doi.org/10.1039/C9CP05091D)

jinzhe.zeng@rutgers.edu (Jinzhe Zeng), tzhu@lps.ecnu.edu.cn (Tong Zhu)

## Features

-   Processing of MD trajectory containing atomic coordinates or bond orders
-   Hidden Markov Model (HMM) based noise filtering
-   Isomers identifying accoarding to SMILES
-   Generation of reaction network for visualization using force-directed algorithm
-   Parallel computing

## Installation

You can [install Anaconda or Miniconda](https://conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html) to obtain conda, and install ReacNetGenerator easily with conda:

```bash
conda install reacnetgenerator -c conda-forge
```

See [the build guide](guide/build.md) if you want to build ReacNetGenerator by yourself. 

## Usage

### Command line

ReacNetGenerator can process any kind of trajectory files containing atomic coordinates, e.g. a LAMMPS dump file prepared by running “dump 1 all custom 100 dump.reaxc id type x y z” in LAMMPS:

```bash
reacnetgenerator --dump -i dump.reaxc -a C H O
```
where C, H, and O are atomic names in the input file. [Analysis report](https://reacnetgenerator.njzjz.win/report.html?jdata=https%3A%2F%2Fgist.githubusercontent.com%2Fnjzjz%2Fe9a4b42ceb7d2c3c7ada189f38708bf3%2Fraw%2F83d01b9ab1780b0ad2d1e7f934e61fa113cb0f9f%2Fmethane.json) will be generated automatically.

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

## Acknowledge
* National Natural Science Foundation of China (Grants No. 91641116)
* National Innovation and Entrepreneurship Training Program for Undergraduate (201910269080)
* ECNU Multifunctional Platform for Innovation (No. 001)
