---
name: reacnetgenerator
description: Analyze reactive molecular dynamics trajectories and generate preliminary reaction networks with ReacNetGenerator.
compatibility: Requires uv and access to the internet.
metadata:
  author: your-name
  version: '1.0'
---

# reacnetgenerator

This skill uses ReacNetGenerator to analyze reactive molecular dynamics trajectories, identify species and reaction events, and generate preliminary reaction networks.

## Agent responsibilities

1. Determine the trajectory file and infer the correct input type (`dump`, `lammpsbondfile`, or `xyz`).
2. Obtain the ordered element list for `-a/--atomname` from a trusted source such as a LAMMPS data file or mapping file.
3. Choose simple defaults unless the user explicitly requests otherwise:
   - keep HMM enabled by default
   - `--stepinterval 1`
   - `--split 1`
   - `--maxspecies 20`
4. Build the ReacNetGenerator command yourself (do **not** ask the user to write the command).
5. Run the command with `uvx --from reacnetgenerator reacnetgenerator ...`.
6. Report what was analyzed, which parameters were used, and where outputs were written.

## Ask the user in plain language

Only ask if key information is missing. Keep questions simple, for example:

- Which trajectory file should be analyzed?
- Is this a LAMMPS dump file, a bond file, or an XYZ trajectory?
- Do you already have a LAMMPS data file or mapping file so I can determine the atom order?
- Do you want the default HMM-filtered analysis, or do you want raw unfiltered events?
- Do you want the whole network, or only a smaller network with fewer species?
- For XYZ trajectories, should periodic boundary conditions be used, and if so, what is the cell?

## Generate the command from user input

Translate user input into:

- input trajectory file
- `--type`
- ordered `-a/--atomname`
- optional flags such as `--nohmm`, `--split`, `--maxspecies`, or `-c`

### Preferred defaults

When the user asks for a standard trajectory analysis, prefer these defaults:

- keep HMM enabled
- `--stepinterval 1`
- `--split 1`
- `--maxspecies 20`

### Important rule for `-a/--atomname`

The ordered element list must come from a trusted upstream source and must preserve atom-type order.

Do **not** sort elements alphabetically.  
Do **not** guess from a documentation example like `C H O Cl`.

## Required commands

```bash
# 1) Check CLI help
uvx --from reacnetgenerator reacnetgenerator --help

# 2) Run a standard analysis
uvx --from reacnetgenerator reacnetgenerator --type dump -i dump.lammpstrj -a C H O Cl

# 3) Run without HMM
uvx --from reacnetgenerator reacnetgenerator --type dump -i dump.lammpstrj -a C H O Cl --nohmm