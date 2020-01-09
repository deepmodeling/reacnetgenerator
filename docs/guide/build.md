# Build ReacNetGenerator

To build ReacNetGenerator by yourself, first, you need to download the source code from [the Releases page](https://github.com/tongzhugroup/reacnetgenerator/releases). Then install ReacNetGenerator with one of the following guides:

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

## Installing via pip
1. Install [OpenBabel](https://github.com/openbabel), [RDKit](https://github.com/rdkit/rdkit), and [Yarn](https://github.com/yarnpkg/yarn).
2. Decompress reacnetgenerator.zip and use `pip` to install in the main directory of ReacNetGenerator. Note that a C/C++ compiler must be installed.
```bash
pip install .
```

