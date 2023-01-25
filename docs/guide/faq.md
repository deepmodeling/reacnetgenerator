# FAQ

## Can ReacNetGenerator process my trajectory?

Generally ReacNetGenerator has no limitation on any specific elements.

When perceiving bond orders from atomic coordinates, ReacNetGenerator will call [Open Babel][openbabel] which reads [element parameters](https://github.com/openbabel/openbabel/blob/2f34bda337d7ddefa8f2bebfc23931a63e45241f/src/elementtable.h).{cite:p}`O'Boyle_JCheminform_2011_v3_p33`
The parameters may not fit the system.
If you have new ideas about parameters, you can report to [Open Babel][openbabel] and recompile Open Babel from the new source code.

When processing a ReaxFF bond file, bond orders are directly provided by ReaxFF with decimal rounding.{cite:p}`Aktulga_ParallelComputing_2012_v38_p245`
The accuracy of the bond orders depends on the accuracy of the force field.

You can refer the list of [publications driven by ReacNetGenerator](https://njzjz.win/reacnetgenerator/) to see other researchers' applications.

## Out of memory (OOM)

When processing a large trajectory, you may get different OOM errors such as `Memory Error`, or subprocesses are directly killed by the system with `broken pipe` thrown.
If you have a machine that has more memory, just use it.
Otherwise, try to reduce the size of the trajectory or split the trajectory into multiple files or increase the number given by `--stepinterval`.
It is also helpful to recuding the number of processes.

If you are using a Windows OS, it's known that the program may consume large memory through multiprocessing.
In this situation, it's suggested to use [Windows Subsystem Linux (WSL)](https://docs.microsoft.com/windows/wsl).

[openbabel]: https://github.com/openbabel/openbabel
