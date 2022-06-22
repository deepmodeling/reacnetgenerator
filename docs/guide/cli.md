# Command line interface

ReacNetGenerator can process any kind of trajectory files containing atomic coordinates, e.g. a LAMMPS dump file prepared by running “dump 1 all custom 100 dump.reaxc id type x y z” in LAMMPS:

```bash
reacnetgenerator --type lammpsdumpfile -i dump.reaxc -a C H O --nohmm
```
where C, H, and O are atomic names in the input file. <a href="/report.html?jdata=https%3A%2F%2Fgist.githubusercontent.com%2Fnjzjz%2Fe9a4b42ceb7d2c3c7ada189f38708bf3%2Fraw%2F83d01b9ab1780b0ad2d1e7f934e61fa113cb0f9f%2Fmethane.json" target="_blank">Analysis report</a> will be generated automatically.

Also, ReacNetGenerator can process files containing bond information, e.g. LAMMPS bond file:

```bash
reacnetgenerator --type lammpsbondfile -i bonds.reaxc -a C H O --nohmm
```

You can running the following script for help:

```bash
reacnetgenerator -h
```

```{argparse}
---
module: reacnetgenerator.commandline
func: main_parser
prog: reacnetgenerator
---
```

