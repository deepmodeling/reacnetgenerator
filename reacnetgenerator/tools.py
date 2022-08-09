"""Useful methods to futhur process ReacNetGenerator results."""
from collections import defaultdict
from typing import Tuple, Dict, List
from collections import Counter

import numpy as np
import ase


def read_species(specfile: str) -> Tuple[List[int], Dict[str, np.ndarray]]:
    """Read species from the species file (ends with .species).

    For accuracy, HMM filter should be disabled.

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


def read_reactions(reacfile) -> List[Tuple[int, Counter, str]]:
    """Read reactions from the reactions file (ends with .reactionsabcd).

    For accuracy, HMM filter should be disabled.
    
    Parameters
    ----------
    reacfile : str
        The reactions file.

    Returns
    -------
    occs : List[Tuple[int, Counter, str]]
        The number of occurences of each reaction. The tuple is (occurence, counter_reactants, reaction).
    """
    occs = []
    with open(reacfile) as f:
        for line in f:
            s = line.split()
            occs.append((int(s[0]), Counter(s[1].split('->')[0].split('+')), s[1]))
    return occs


def calculate_rate(specfile: str, reacfile: str, cell: np.ndarray, timestep: float) -> Dict[str, float]:
    """Calculate the rate constant of each reaction.

    The rate constants are calculated by the method developed in [1].

    Parameters
    ----------
    specfile : str
        The species file.
    reacfile : str
        The reactions file.
    cell : np.ndarray
        The cell with the shape (3, 3). Unit: Angstrom.
    timestep : float
        The time step. Unit: femtosecond.

    Returns
    -------
    rates : Dict[str, float]
        The rate of each reaction. The dict key is the reaction SMILES.
        The value is in unit of [(cm^3/mol)s^(-1)].

    References
    ----------
    .. [1] J Comput Chem 40, 16, 1586-1592.
    """
    cell = ase.geometry.Cell(cell)
    
    timestep *= 10**-15 # fs to s
    #N, step_tot =read_species(specfile)
    step_idx, n_species = read_species(specfile)
    occs = read_reactions(reacfile)

    # total time during simulation
    time_tot = (step_idx[-1] - step_idx[0]) * timestep
    # volume of the cell
    volume = cell.volume
    volume *= 10**-24 # Ang^3 to cm^3
    volume_times_na = volume * ase.units.mol # V * NA
    
    rates = {}
    for occ, reacts, reactions in occs:
        # k = occ_tot / ( V * time_tot * c_tot )
        # c_tot = N_tot / (V * NA)
        n_react = np.array([n_species[kk] for kk in reacts.keys()])
        nu = np.array(list(reacts.values()))
        c_po = np.power(n_react / volume_times_na, np.repeat(nu, n_react.shape[1]).reshape(n_react.shape))
        c_tot = np.sum(np.prod(c_po, axis=0))
        k = occ / (volume_times_na * time_tot * c_tot)
        rates[reactions] = k
    return rates
