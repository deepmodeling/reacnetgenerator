---
home: true
heroImage: /reacnetgen.svg
heroText: ReacNetGenerator
tagline: An automatic generator of reaction network for reactive molecular dynamics simulation
actionText: Download
actionLink: http://computchem.cn/reacnetgenerator/
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

**Please cite:** J. Zeng, L. Cao, J.Z.H. Zhang, C.H. Chin, T. Zhu: ReacNetGen: an Automatic Reaction Network Generator for Reactive Molecular Dynamic Simulations, 2018, doi: [10.26434/chemrxiv.7421534](https://dx.doi.org/10.26434/chemrxiv.7421534)

**Author**:
[Jinzhe Zeng](https://cv.njzjz.win),
[Liqun Cao](http://computchem.cn/people/),
[John ZH Zhang](https://research.shanghai.nyu.edu/centers-and-institutes/chemistry/people/john-zenghui-zhang),
[Chih-Hao Chin](http://computchem.cn/people/),
[Tong Zhu](http://computchem.cn/people/)

**Email**: tzhu@lps.ecnu.edu.cn, jzzeng@stu.ecnu.edu.cn

# Installation

1. Download the source code on [our group website](http://computchem.cn/reacnetgenerator/).
2. [Install Anaconda or Miniconda](https://conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html) to obtain conda.
3. Decompress reacnetgenerator.zip and build in the main directory:

```bash
conda config --add channels conda-forge
conda build conda/recipe
conda install reacnetgenerator --use-local
```

# Usage

## Command line

ReacNetGenerator can process any kind of trajectory files containing atomic coordinates, e.g. a LAMMPS dump file prepared by running “dump 1 all custom 100 dump.reaxc id type x y z” in LAMMPS:

```bash
reacnetgenerator --dump -i dump.reaxc -a C H O
```
where C, H, and O are atomic names in the input file. [Analysis report](https://njzjz.github.io/reacnetgenerator/report.html) will be generated automatically.

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

# Dependencies

-   Python >= 3.6
-   Python packages:
    [numpy](https://github.com/numpy/numpy),
    [scipy](https://github.com/scipy/scipy),
    [pandas](https://github.com/pandas-dev/pandas),
    [networkx](https://github.com/networkx/networkx),
    [scikit-learn](https://github.com/scikit-learn/scikit-learn),
    [matplotlib](https://github.com/matplotlib/matplotlib),
    [hmmlearn](https://github.com/hmmlearn/hmmlearn),
    [ASE](https://gitlab.com/ase/ase),
    [scour](https://github.com/scour-project/scour),
    [tqdm](https://github.com/tqdm/tqdm),
    [jinja2](https://github.com/pallets/jinja),
    [coloredlogs](https://github.com/xolox/python-coloredlogs),
    [htmlmin](https://github.com/mankyd/htmlmin/),
    [lz4](https://github.com/python-lz4/python-lz4),
    [pybase64](https://github.com/mayeut/pybase64)
-   Extra libraries:
    [Yarn](https://github.com/yarnpkg/yarn),
    [OpenBabel](https://github.com/openbabel/openbabel),
    [RDKit](https://github.com/rdkit/rdkit)
-   npm packages:
    [jQuery](https://github.com/jquery/jquery),
    [jQuery Easing Plugin](https://github.com/gdsmith/jquery.easing),
    [Magnific Popup](https://github.com/dimsemenov/Magnific-Popup),
    [Bootstrap](https://github.com/twbs/bootstrap),
    [Start Bootstrap - Creative](https://github.com/BlackrockDigital/startbootstrap-creative),
    [D3](https://github.com/d3/d3),
    [JSNetworkX](https://github.com/fkling/JSNetworkX)
