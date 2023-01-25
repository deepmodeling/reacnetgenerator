# Installation

:::{note}
The latest version requires Python 3.7 or later.
:::

## Install via conda

You can [install Anaconda or Miniconda](https://conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html) to obtain conda, and install ReacNetGenerator easily with conda:

```bash
conda install reacnetgenerator -c conda-forge
```

## Install via pip

```bash
# upgrade pip as old pip may not be supported
pip install -U pip
pip install reacnetgenerator
```

See also [Tutorial: Installation](../tutorial/install.ipynb).

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
