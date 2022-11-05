# <img src=docs/_static/reacnetgen.svg height=40/>  ReacNetGenerator

[![DOI:10.1039/C9CP05091D](https://img.shields.io/badge/DOI-10.1039%2FC9CP05091D-blue)](https://doi.org/10.1039/C9CP05091D)
[![Citations](https://citations.njzjz.win/10.1039/C9CP05091D)](https://doi.org/10.1039/C9CP05091D)
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

## Guide and Tutorial

The latest version requires Python 3.7 or later.
You can install ReacNetGenerator with `conda`:

```sh
conda install reacnetgenerator -c conda-forge
reacnetgenerator -h
```

See [the guide](https://reacnetgenerator.njzjz.win/guide/) to learn how to install and use ReacNetGenerattor. We also provide [a series of tutorials](https://reacnetgenerator.njzjz.win/tutorial/) to help you learn ReacNetGenerator.

## Awards
* The First Prize in 2019 (the 11th Session) Shanghai Computer Application Competition for College Students
* The First Prize in 2019 (the 12th Session) Chinese Computer Design Competition for College Students

## Acknowledge
* National Natural Science Foundation of China (Grants No. 91641116)
* National Innovation and Entrepreneurship Training Program for Undergraduate (201910269080)
* ECNU Multifunctional Platform for Innovation (No. 001)
