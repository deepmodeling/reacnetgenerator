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

```
from reacnetgenerator import run

result = run(
    input_path="trajectory.lammpstrj",
    output_dir="artifacts",
    input_type="dump",
    atomname=["C", "H", "O"],
    items=("species", "reactions", "network", "report"),
    runHMM=False,
)
print(result["artifacts"]["species"])
print(result["provenance"]["parameters"])
```

The existing ReacNetGenerator class accepts the same output_dir keyword and
exposes the mapping as generator.artifacts. The convenience result also records
the normalized/defaulted parameters, explicitly supplied parameter names, and
requested items as JSON-serializable provenance.

## Time-resolved output

Setting `printmoleculetime=True` or `printreactionevent=True` writes the existing
CSV output and a schema 1.1 HDF5 timeline. The default HDF5 path is
`<input>.timeline.h5`, or `timeline.h5` below `output_dir`. Pass
`timed_output="custom.h5"` to override it or to request a timeline without the
CSV switches. See [HDF5 timeline](timed-output-schema.md) for the public readers
and validator.
