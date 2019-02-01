# ReacNetGenerator

[![python version](https://img.shields.io/pypi/pyversions/reacnetgenerator.svg?logo=python&logoColor=white)](https://pypi.org/project/reacnetgenerator)
[![PyPI](https://img.shields.io/pypi/v/reacnetgenerator.svg)](https://pypi.org/project/reacnetgenerator)
[![Build Status](https://travis-ci.com/njzjz/reacnetgenerator.svg?branch=master)](https://travis-ci.com/njzjz/reacnetgenerator)
[![Build status](https://ci.appveyor.com/api/projects/status/t92hud34lgel1eu6?svg=true)](https://ci.appveyor.com/project/jzzeng/reacnetgenerator)
[![Coverage Status](https://coveralls.io/repos/github/njzjz/reacnetgenerator/badge.svg?branch=master)](https://coveralls.io/github/njzjz/reacnetgenerator?branch=master)
[![codecov](https://codecov.io/gh/njzjz/reacnetgenerator/branch/master/graph/badge.svg)](https://codecov.io/gh/njzjz/reacnetgenerator)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/b2336e2a2ff04aceab42604792c1c3e1)](https://www.codacy.com/app/jzzeng/reacnetgenerator?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=njzjz/reacnetgenerator&amp;utm_campaign=Badge_Grade)
[![DOI:10.26434/chemrxiv.7421534](https://zenodo.org/badge/DOI/10.26434/chemrxiv.7421534.svg)](https://doi.org/10.26434/chemrxiv.7421534)

An automatic generator of reaction network for reactive molecular dynamics simulation.

**Please cite:** J. Zeng, L. Cao, J.Z.H. Zhang, C.H. Chin, T. Zhu: ReacNetGen: an Automatic Reaction Network Generator for Reactive Molecular Dynamic Simulations, 2018, doi: [10.26434/chemrxiv.7421534](https://dx.doi.org/10.26434/chemrxiv.7421534)

**Author**: [Jinzhe Zeng](https://cv.njzjz.win), Liqun Cao, [John ZH Zhang](https://research.shanghai.nyu.edu/centers-and-institutes/chemistry/people/john-zenghui-zhang), Chih-Hao Chin, [Tong Zhu](http://computchem.cn/people/)

**Email**: tzhu@lps.ecnu.edu.cn, jzzeng@stu.ecnu.edu.cn

[![Research Group](https://img.shields.io/website-up-down-green-red/http/computchem.cn.svg?label=Research%20Group)](http://computechem.cn)

## Features

-   Processing of MD trajectory containing atomic coordinates or bond orders
-   Hidden Markov Model (HMM) based noise filtering
-   Isomers identifying accoarding to SMILES
-   Generation of reaction network for visualization using force-directed algorithm
-   Parallel computing

## Requirements

-   Python 3.6 - 3.7 (**Note:** Python &lt;= 3.5 is not supported!)
-   Python packages: [numpy](https://github.com/numpy/numpy), [scipy](https://github.com/scipy/scipy), [pandas](https://github.com/pandas-dev/pandas), [networkx](https://github.com/networkx/networkx), [scikit-learn](https://github.com/scikit-learn/scikit-learn), [matplotlib](https://github.com/matplotlib/matplotlib), [hmmlearn](https://github.com/hmmlearn/hmmlearn), [ASE](https://gitlab.com/ase/ase), [scour](https://github.com/scour-project/scour), [tqdm](https://github.com/tqdm/tqdm), [jinja2](https://github.com/pallets/jinja), [coloredlogs](https://github.com/xolox/python-coloredlogs), [htmlmin](https://github.com/mankyd/htmlmin/), [jsmin](https://github.com/tikitu/jsmin/), [cssmin](https://github.com/zacharyvoase/cssmin)
-   Extra libraries: [OpenBabel](https://github.com/openbabel/openbabel), [RDKit](https://github.com/rdkit/rdkit)
-   Javascript and CSS frameworks (already packaged): [jQuery](https://github.com/jquery/jquery), [jQuery Easing Plugin](https://github.com/gdsmith/jquery.easing), [Magnific Popup](https://github.com/dimsemenov/Magnific-Popup), [ScrollReveal](https://github.com/scrollreveal/scrollreveal), [Bootstrap](https://github.com/twbs/bootstrap), [Start Bootstrap - Creative](https://github.com/BlackrockDigital/startbootstrap-creative), [D3](https://github.com/d3/d3), [JSNetworkX](https://github.com/fkling/JSNetworkX)

## Installation

1.  [Get conda](https://conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html) and install Anaconda or Miniconda.
2.  Use conda to create environment:

```sh
conda create -q -n reacnetgenerator python=3.7 openbabel rdkit hmmlearn -c openbabel -c conda-forge
source activate reacnetgenerator # for Windows, just use: activate reacnetgenerator
```

3.  Download ReacNetGenerator and install it from source:

```sh
git clone https://github.com/njzjz/reacnetgenerator
cd reacnetgenerator
pip install .
```

You can test whether ReacNetGenerator is running normally:

```sh
python setup.py pytest
```

## Simple example

Prepare a [LAMMPS bond file](http://lammps.sandia.gov/doc/fix_reax_bonds.html) named bonds.reaxc, then run the script:

```sh
reacnetgenerator -i bonds.reaxc -a C H O
```

where C, H, and O are atomic names in the input file. [Analysis report](report.html) will be generated automatically.  

A [LAMMPS dump file](https://lammps.sandia.gov/doc/dump.html) is also supported. You can prepare it by running "dump 1 all custom 100 dump.reaxc id type x y z" in LAMMPS.

```sh
reacnetgenerator --dump -i dump.reaxc -a C H O
```

You can running the following script for help:

```sh
reacnetgenerator -h
```

## GUI version

You can open a GUI version for ReacNetGenerator by typing:

```sh
reacnetgeneratorgui
```
