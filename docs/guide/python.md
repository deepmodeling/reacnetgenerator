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
