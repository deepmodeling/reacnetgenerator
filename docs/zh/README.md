---
home: true
heroImage: /reacnetgen.svg
heroText: ReacNetGenerator
tagline: 基于反应分子动力学模拟轨迹的反应网络自动构建和可视化软件
actionText: 下载
actionLink: https://github.com/tongzhugroup/reacnetgenerator/releases
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

<video width="512" height="397.6" controls>
  <source src="http://www.rsc.org/suppdata/c9/cp/c9cp05091d/c9cp05091d2.mp4" type="video/mp4">
</video> 

# 引用和联系方式

ReacNetGenerator: an automatic reaction network generator for reactive molecular dynamic simulations, Phys. Chem. Chem. Phys., 2020, 22 (2): 683–691, doi: [10.1039/C9CP05091D](https://dx.doi.org/10.1039/C9CP05091D)

jinzhe.zeng@rutgers.edu (Jinzhe Zeng), tzhu@lps.ecnu.edu.cn (Tong Zhu)

# 安装
[从清华大学开源软件镜像站下载并安装 Anaconda 或 Miniconda](https://mirror.tuna.tsinghua.edu.cn/help/anaconda/)，然后用conda安装ReacNetGenerator:

```bash
conda install reacnetgenerator -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
```

如果你想自己构建ReacNetGenerator，参见[构建教程](guide/build.md)。

# 使用

## 命令行

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

## 图形界面

你可以打开 ReacNetGenerator 的图形界面：

```bash
reacnetgeneratorgui
```

# 荣誉
* 2019年（第十一届）上海市大学生计算机应用能力大赛一等奖
* 2019年（第12届）中国大学生计算机设计大赛一等奖

# Acknowledge
* 国家自然科学基金（编号91641116）
* 国家大学生创新创业训练计划项目（201910269080）
* 华东师范大学仪器共享平台（001号）
