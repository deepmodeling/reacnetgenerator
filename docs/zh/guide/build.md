# 构建ReacNetGenerator

首先在[Releases](https://github.com/tongzhugroup/reacnetgenerator/releases)下载reacnetgenerator.zip。接着选择一项安装ReacNetGenerator：

## 构建conda包
1. [从清华大学开源软件镜像站下载并安装 Anaconda 或 Miniconda](https://mirror.tuna.tsinghua.edu.cn/help/anaconda/) ；
2. 解压reacnetgenerator.zip，并在ReacNetGenerator主目录编译：

```bash
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda build conda/recipe
conda install reacnetgenerator --use-local
```

## 构建Docker镜像
1. [安装Docker](https://mirror.tuna.tsinghua.edu.cn/help/docker-ce/)；
2. 解压reacnetgenerator.zip，并在ReacNetGenerator主目录编译：

```bash
docker build . -t njzjz/reacnetgenerator
docker run njzjz/reacnetgenerator reacnetgenerator -h
```

## 用pip安装
1. 安装[OpenBabel](https://github.com/openbabel)、[RDKit](https://github.com/rdkit/rdkit)和[Yarn](https://github.com/yarnpkg/yarn)。
2. 解压reacnetgenerator.zip，并在ReacNetGenerator主目录用`pip`安装。注意系统中必须有C/C++编译器。
```bash
pip install .
```

