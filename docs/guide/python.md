# Python interface

## Running the ReacNetGenerator

You can use the Python interface:

```py
from reacnetgenerator import ReacNetGenerator
ReacNetGenerator(
  inputfiletype="dump",
  inputfilename="dump.ch4",
  atomname=['C', 'H', 'O'],
  ).runanddraw()
```

See {class}`ReacNetGenerator <reacnetgenerator.ReacNetGenerator>` class for detailed parameters.

## Calculate rate constants

An effiective tool is provided in {meth}`reacnetgenerator.tools.calculate_rate <reacnetgenerator.tools.calculate_rate>` to calculate rate constants.

## Explicit output directory

The convenience API returns semantic artifact paths and writes its default
outputs below the requested directory.

    from reacnetgenerator import run

    artifacts = run(
        input_path="trajectory.lammpstrj",
        output_dir="artifacts",
        input_type="dump",
        atomname=["C", "H", "O"],
        items=("species", "reactions", "network", "report"),
        runHMM=False,
    )
    print(artifacts["species"])

The existing ReacNetGenerator class accepts the same output_dir keyword and
exposes the mapping as generator.artifacts.
