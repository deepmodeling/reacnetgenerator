---
home: true
heroImage: /reacnetgen.svg
heroText: ReacNetGenerator
tagline: 反应分子动力学模拟的反应网络自动生成器
actionText: 分析结果
actionLink: report.html
features:
- title: 轨迹
  details: 处理包含原子坐标或键级的轨迹
- title: 过滤
  details: 基于隐马尔科夫模型的噪声过滤
- title: 同分异构体
  details: 根据 SMILES 识别同分异构体
- title: 网络
  details: 使用力引导算法生成反应网络
- title: HTML5
  details: 展示可交互的网页结果
- title: 快速
  details: 包含并行计算和性能优化
  footer: Copyright © 2018-2019 East China Normal University
---

# 引用和联系方式

**请引用：** J. Zeng, L. Cao, J.Z.H. Zhang, C.H. Chin, T. Zhu: ReacNetGen: an Automatic Reaction Network Generator for Reactive Molecular Dynamic Simulations, 2018, doi: [10.26434/chemrxiv.7421534](https://dx.doi.org/10.26434/chemrxiv.7421534)

**作者**：
[Jinzhe Zeng](https://cv.njzjz.win),
[Liqun Cao](http://computchem.cn/people/),
[John ZH Zhang](https://research.shanghai.nyu.edu/centers-and-institutes/chemistry/people/john-zenghui-zhang),
[Chih-Hao Chin](http://computchem.cn/people/),
[Tong Zhu](http://computchem.cn/people/)

**Email**： tzhu@lps.ecnu.edu.cn, jzzeng@stu.ecnu.edu.cn

## 安装

[从清华大学开源镜像站下载 Anaconda 或 Miniconda](https://mirror.tuna.tsinghua.edu.cn/help/anaconda/) 并安装，然后：

```bash
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda build conda/recipe
conda install reacnetgenerator --use-local
```

## 使用

### 命令行

准备名为 bonds.reaxc 的 [LAMMPS 键级文件](http://lammps.sandia.gov/doc/fix_reax_bonds.html)，输入：

```bash
reacnetgenerator -i bonds.reaxc -a C H O
```

这里，C、H、O 是轨迹中的原子种类。软件将自动生成[分析结果](report.html)。 

软件也支持 [LAMMPS dump 文件](https://lammps.sandia.gov/doc/dump.html)。你可以在 LAMMPS 中执行`dump 1 all custom 100 dump.reaxc id type x y z`来获得。

```bash
reacnetgenerator --dump -i dump.reaxc -a C H O
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

## 依赖

-   Python 3.6 - 3.7
-   Python 包：
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
-   额外的软件依赖：
    [Yarn](https://github.com/yarnpkg/yarn),
    [OpenBabel](https://github.com/openbabel/openbabel),
    [RDKit](https://github.com/rdkit/rdkit)
-   npm 包：
    [jQuery](https://github.com/jquery/jquery),
    [jQuery Easing Plugin](https://github.com/gdsmith/jquery.easing),
    [Magnific Popup](https://github.com/dimsemenov/Magnific-Popup),
    [ScrollReveal](https://github.com/scrollreveal/scrollreveal),
    [Bootstrap](https://github.com/twbs/bootstrap),
    [Start Bootstrap - Creative](https://github.com/BlackrockDigital/startbootstrap-creative),
    [D3](https://github.com/d3/d3),
    [JSNetworkX](https://github.com/fkling/JSNetworkX)