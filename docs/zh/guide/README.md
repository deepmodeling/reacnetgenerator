# 指南

ReacNetGenerator是基于反应分子动力学模拟轨迹的反应网络自动构建和可视化软件。本页面提供了快速上手的指南。

## 安装
[从清华大学开源软件镜像站下载并安装 Anaconda 或 Miniconda](https://mirror.tuna.tsinghua.edu.cn/help/anaconda/)，然后用conda安装ReacNetGenerator:

```bash
conda install reacnetgenerator -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
```

如果你想自己构建ReacNetGenerator，参见[构建教程](guide/build.md)。

## 使用

### 命令行

ReacNetGenerator可以处理任意类型含有原子坐标的轨迹文件，例如LAMMPS dump文件，可以通过在 LAMMPS 中执行dump 1 all custom 100 dump.reaxc id type x y z来获得：

```bash
reacnetgenerator --dump -i dump.reaxc -a C H O
```

其中，C、H、O 是轨迹中的原子种类。软件将自动生成<a href="/report.html?jdata=https%3A%2F%2Fgist.githubusercontent.com%2Fnjzjz%2Fe9a4b42ceb7d2c3c7ada189f38708bf3%2Fraw%2F83d01b9ab1780b0ad2d1e7f934e61fa113cb0f9f%2Fmethane.json" target="_blank">分析结果</a>。

软件也可以处理含有键级信息的文件，例如LAMMPS键级文件：

```bash
reacnetgenerator -i bonds.reaxc -a C H O
```

运行以下命令查看帮助：

```bash
reacnetgenerator -h
```

### 图形界面

你可以打开 ReacNetGenerator 的图形界面：

```bash
reacnetgeneratorgui
```