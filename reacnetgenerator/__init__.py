# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright 2018-2022, East China Normal University
"""ReacNetGenerator is an automatic reaction network generator for
reactive molecular dynamics simulation.[1]_.

References
----------
.. [1] Jinzhe Zeng, Liqun Cao, Chih-Hao Chin, Haisheng Ren, John Z. H.
   Zhang, Tong Zhu, ReacNetGenerator: an automatic reaction network
   generator for reactive molecular dynamic simulations, Phys. Chem.
   Chem. Phys., 2020, 22 (2): 683-691, doi: 10.1039/C9CP05091D.
"""

__date__ = "2018-03-11"
__author__ = "Jinzhe Zeng"
__email__ = "jinzhe.zeng@rutgers.edu"
__credits__ = ["Jinzhe Zeng", "Tong Zhu", "Liqun Cao", "Chih-Hao Chin", "John ZH Zhang"]
__copyright__ = "Copyright 2018-2022, East China Normal University"

import matplotlib as mpl

from ._version import __version__

mpl.use("svg")

from .reacnetgen import ReacNetGenerator

__all__ = ["ReacNetGenerator", "__version__"]
