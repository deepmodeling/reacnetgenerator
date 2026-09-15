# Command line interface

```{argparse}
---
module: reacnetgenerator.commandline
func: main_parser
prog: reacnetgenerator
---
```

## Explicit output directory

Pass --output-dir to keep all default-generated artifacts for one run in an
independent directory. The input basename is not used to distinguish outputs.

```
reacnetgenerator --type dump -i trajectory.lammpstrj -a C H O \
    --nohmm --output-dir artifacts
```

The option changes only default output paths; explicit filename arguments remain
authoritative for backward compatibility.
