# Installation

## Install via conda

You can [install Anaconda or Miniconda](https://conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html) to obtain conda, and install ReacNetGenerator easily with conda:

```bash
conda install reacnetgenerator -c conda-forge
```

## Install via pip

[![Nbviewer](https://img.shields.io/badge/render-nbviewer-orange)](https://nbviewer.jupyter.org/github/tongzhugroup/reacnetgenerator/blob/master/tutorial/install.ipynb?flush_cache=false)
[![Colab](https://images.weserv.nl/?url=colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/tongzhugroup/reacnetgenerator/blob/master/tutorial/install.ipynb)

```bash
pip install reacnetgenerator
```

## Docker images

If you have [installed Docker](https://docs.docker.com/install/), an official Docker image is provided: 

```bash
docker run ghcr.io/tongzhugroup/reacnetgenerator reacnetgenerator -h
```

When analyzing trajectories, you need to [monut](https://docs.docker.com/storage/bind-mounts/) the local directory into the container.

If your HPC node has installed [Singularity](https://sylabs.io/docs/), you can also used it:

```bash
singularity run docker://ghcr.io/tongzhugroup/reacnetgenerator reacnetgenerator -h
```

:::{note}
See [the build guide](build.md) if you want to build ReacNetGenerator from source. 
:::
