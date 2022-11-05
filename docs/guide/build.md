# Build ReacNetGenerator

To build ReacNetGenerator by yourself, first, you need to download the source code from [the Releases page](https://github.com/tongzhugroup/reacnetgenerator/releases) or use `git` to clone the latest code:

```bash
git clone https://github.com/tongzhugroup/reacnetgenerator
cd reacnetgenerator
```

Then build ReacNetGenerator with one of the following guides:

## Installing via pip

Use `pip` to install in the main directory of ReacNetGenerator. Note that a C/C++ compiler must be installed.
```bash
# upgrade pip as old pip may not be supported
pip install -U pip
pip install .
```

Test installation by
```bash
reacnetgenerator -h
```

## Building a conda package

[Install Anaconda or Miniconda](https://conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html) to obtain conda and build in the main directory of ReacNetGenerator:

```bash
conda config --add channels conda-forge
conda build conda/recipe
conda install reacnetgenerator --use-local
reacnetgenerator -h
```

## Building a Docker Image

[Install Docker](https://docs.docker.com/install/) and build in the main directory of ReacNetGenerator:

```bash
docker build . -t njzjz/reacnetgenerator
docker run njzjz/reacnetgenerator reacnetgenerator -h
```
