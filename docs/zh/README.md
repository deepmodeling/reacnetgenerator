---
home: true
heroImage: /reacnetgen.svg
heroText: ReacNetGenerator
tagline: 基于反应分子动力学模拟轨迹的反应网络自动构建和可视化软件
actionText: 下载
actionLink: http://computchem.cn/reacnetgenerator
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

1. 在[我们的课题组网站](http://computchem.cn/reacnetgenerator)下载reacnetgenerator.zip；
2. [从清华大学开源软件镜像站下载并安装 Anaconda 或 Miniconda](https://mirror.tuna.tsinghua.edu.cn/help/anaconda/) ；
3. 解压reacnetgenerator.zip，并在软件主目录编译：

```bash
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda build conda/recipe
conda install reacnetgenerator --use-local
```

## 使用

### 命令行

ReacNetGenerator可以处理任意类型含有原子坐标的轨迹文件，例如LAMMPS dump文件，可以通过在 LAMMPS 中执行dump 1 all custom 100 dump.reaxc id type x y z来获得：

```bash
reacnetgenerator --dump -i dump.reaxc -a C H O
```

其中，C、H、O 是轨迹中的原子种类。软件将自动生成[分析结果](https://njzjz.github.io/reacnetgenerator/report.html?jdata=https%3A%2F%2Fgist.githubusercontent.com%2Fnjzjz%2Fe9a4b42ceb7d2c3c7ada189f38708bf3%2Fraw%2F83d01b9ab1780b0ad2d1e7f934e61fa113cb0f9f%2Fmethane.json)。

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

## 依赖

-   Python >= 3.6
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
    [coloredlogs](https://github.com/xolox/python-coloredlogs),
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
    [Bootstrap](https://github.com/twbs/bootstrap),
    [Start Bootstrap - Creative](https://github.com/BlackrockDigital/startbootstrap-creative),
    [D3](https://github.com/d3/d3),
    [JSNetworkX](https://github.com/fkling/JSNetworkX)

# 荣誉
* 2019年（第十一届）上海市大学生计算机应用能力大赛一等奖
* 2019年（第12届）中国大学生计算机设计大赛一等奖