"""Useful methods to futhur process ReacNetGenerator results."""
from collections import defaultdict
from typing import Tuple, Dict, List

import numpy as np


def read_species(specfile: str) -> Tuple[List[int], Dict[str, np.ndarray]]:
    """Read species from the species file (ends with .species).

    Parameters
    ----------
    specfile : str
        The species file.

    Returns
    -------
    step_idx : np.ndarray
        The index of the step.
    n_species : Dict[str, np.ndarray]
        The number of species in each step. The dict key is the species SMILES.

    Examples
    --------
    Plot the number of methane in each step.

    >>> from reacnetgenerator.tools import read_species
    >>> import matplotlib.pyplot as plt
    >>> step_idx, n_species = read_species('methane.species')
    >>> plt.plot(step_idx, n_species['[H]C([H])([H])[H]'])
    >>> plt.savefig("methane.svg")
    """
    step_idx = []
    n_species = defaultdict(lambda: defaultdict(int))
    with open(specfile) as f:
        for ii, line in enumerate(f):
            s = line.split()
            step_idx.append(int(s[1].strip(':')))
            for ss, nn in zip(s[2::2], [int(x) for x in s[3::2]]):
                n_species[ss][ii] = nn
        else:
            nsteps = ii + 1
    n_species2 = {}
    for ss in n_species:
        n_species2[ss] = np.array([n_species[ss][ii] for ii in range(nsteps)], dtype=int)
    return np.array(step_idx, dtype=int), n_species2
