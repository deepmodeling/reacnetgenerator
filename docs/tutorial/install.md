# Install ReacNetGenerator

[![Nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/tongzhugroup/reacnetgenerator/blob/master/tutorial/install.ipynb?flush_cache=false)
[![Colab](https://img.njzjz.win/?url=colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/tongzhugroup/reacnetgenerator/blob/master/tutorial/install.ipynb)

In this tutorial, we will use conda to install ReacNetGenerator. This is the eastiest way without compiling the package. You can click the above badget to repeat it.

## Install Miniconda

Firstly, we download the latest version of Miniconda 3 and then install it in `/usr/local`.

```bash
MINICONDA_INSTALLER_SCRIPT=Miniconda3-latest-Linux-x86_64.sh
MINICONDA_PREFIX=/usr/local
wget https://repo.continuum.io/miniconda/$MINICONDA_INSTALLER_SCRIPT --no-verbose
chmod +x $MINICONDA_INSTALLER_SCRIPT
./$MINICONDA_INSTALLER_SCRIPT -b -f -p $MINICONDA_PREFIX
```

## Install ReacNetGenerator

Then, we use `conda` to install the ReacNetGenerator from the [conda-forge](https://github.com/conda-forge) channel.

```bash
conda install reacnetgenerator -c conda-forge -y
```

Now, ReacNetGenerator has been already installed in the environment. Try it now!

```bash
reacnetgenerator -h
```
