# Copyright 2018-2019, East China Normal University
"""ReacNetGen."""

__date__ = '2018-03-11'
__update__ = '2019-02-10'
__author__ = 'Jinzhe Zeng'
__email__ = 'jzzeng@stu.ecnu.edu.cn'
__credits__ = ['Jinzhe Zeng', 'Tong Zhu',
               'Liqun Cao', 'Chih-Hao Chin', 'John ZH Zhang']
__copyright__ = 'Copyright 2018-2019, East China Normal University'


import logging

import coloredlogs
import matplotlib as mpl
mpl.use("svg")  # noqa
from pkg_resources import DistributionNotFound, get_distribution

from .reacnetgen import ReacNetGenerator


__all__ = ['ReacNetGenerator']

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    __version__ = ''

coloredlogs.install(
    fmt=f'%(asctime)s - ReacNetGen {__version__} - %(levelname)s: %(message)s',
    level=logging.INFO, milliseconds=True)
