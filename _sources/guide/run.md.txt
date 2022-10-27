# Analysis

After a [trajectory file](format.md) is prepared, you can process the file using the command line:

```bash
reacnetgenerator --type dump -i dump.reaxc -a C H O --nohmm
```
where C, H, and O are atomic names in the input file.
`--type` decides the [format](format.md) of the trajectory file `dump.reaxc`.
`--nohmm` controls whether [HMM filter](hmm.md) is enabled.

For example, if you want to process a LAMMPS bond file instead:

```bash
reacnetgenerator --type bond -i bonds.reaxc -a C H O --nohmm
```

[A serial of files](report.md) will be generated. You can start with <a href="/report.html?jdata=https%3A%2F%2Fgist.githubusercontent.com%2Fnjzjz%2Fe9a4b42ceb7d2c3c7ada189f38708bf3%2Fraw%2F83d01b9ab1780b0ad2d1e7f934e61fa113cb0f9f%2Fmethane.json" target="_blank">a web page for Analysis report</a>, which ends with the `.html` suffix.

You can run the following script for help:

```bash
reacnetgenerator -h
```

See [here](cli.md) for the usage of the command line.
