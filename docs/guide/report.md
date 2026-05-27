# Generated files

The following files are reported:

## Web page

suffix: `.html`

The web page including analysis.
<a href="../report.html?jdata=https%3A%2F%2Fgist.githubusercontent.com%2Fnjzjz%2Fe9a4b42ceb7d2c3c7ada189f38708bf3%2Fraw%2F83d01b9ab1780b0ad2d1e7f934e61fa113cb0f9f%2Fmethane.json" target="_blank">Here</a> is an example.
You can open it using a modern browser.
Note that $\ce{A + B -> C + D}$ information may be not accurate when [HMM](hmm.md) is enabled.

## Data file

suffix: `.json`

The JSON file storing necessary data for the web page.
You can load it through the <a href="../report.html" target="_blank">web page loader</a>.

## Species file

suffix: `.species`

This text file stores the number of each species in each time step.
One can read this file through the Python method {meth}`reacnetgenerator.tools.read_species <reacnetgenerator.tools.read_species>`.
Note that when [HMM](hmm.md) is enabled, information in this file may not be accurate.

## Molecule file

suffix: `.moname`

This file contains information of each molecule.
In each line, the first column is its SMILES.
The second column is the atomic index (starts from 0) of atoms.
The last column shows all the bonds in the molecule.
This file always keeps the historical three-column format.

## Molecule timeline file

suffix: `.molecules.csv`

This optional CSV file is written when `--show-molecule-time`,
`--molecule-frame`, or `--molecule-timestep` is enabled. It contains one row per
molecule occurrence in an analyzed frame, with the following columns:

```text
Timestep,Species,AtomIDs,BondIDs
```

Rows are sorted by `Timestep`. `Timestep` is the original timestep value from
the input trajectory. `AtomIDs` and `BondIDs` are semicolon-separated lists.
Bond IDs use `atom1-atom2-order` labels. The
`--molecule-frame` and `--molecule-timestep` filters limit the CSV rows to
selected frames or timesteps.

## Route file

suffix: `.route`

This file contains the route of each atom.
It shows which species an atom is inside thorugh the whole trajectory in the following format:

```
Atom {idx}: {time} {SMILES} -> {time} {SMILES} -> ...
```

## Reaction files

suffix: `.reaction`, `.reactionabcd`

`.reaction` shows the frequency of the reaction $\ce{A -> B}$ while `.reactionabcd` shows the frequency of the reaction $\ce{A + B -> C + D}$.
One can read these files through the Python method {meth}`reacnetgenerator.tools.read_reactions <reacnetgenerator.tools.read_reactions>`.
Note that $\ce{A + B -> C + D}$ information may be not accurate when [HMM](hmm.md) is enabled.

## Reaction event file

suffix: `.reactionevent.csv`

This optional CSV file is written when `--reaction-event` is enabled. It contains
one row per reaction event detected between two adjacent analyzed frames:

```text
Timestep_Index,Reactant,Product
```

`Timestep_Index` is the analyzed transition index. A value of `0` means the
reaction was detected between analyzed frames 0 and 1. `Reactant` and `Product`
use the same species names as the reaction summary files.

## Reaction matrix file

suffix: `.table`

This file shows the reaction matrix.
To save space, only first top 100 species are printed.
You can process [reaction files](#reaction-files) manually for other information.

## Reaction network

suffix: `.svg`

A SVG file to show reaction networks.
