# Copyright 2018-2019, East China Normal University
"""ReacNetGen."""

__date__ = '2018-03-11'
__update__ = '2019-07-29'
__author__ = 'Jinzhe Zeng'
__email__ = 'jzzeng@stu.ecnu.edu.cn'
__credits__ = ['Jinzhe Zeng', 'Tong Zhu',
               'Liqun Cao', 'Chih-Hao Chin', 'John ZH Zhang']
__copyright__ = 'Copyright 2018-2019, East China Normal University'

import matplotlib as mpl
mpl.use("svg")  # noqa
import networkx  # avoid qhull library error

from . import _logging
from ._version import __version__
from .reacnetgen import ReacNetGenerator

__all__ = ['ReacNetGenerator']

# fix coverage for multiprocessing
try:
    from pytest_cov.embed import cleanup_on_sigterm
except ImportError:
    pass
else:
    cleanup_on_sigterm()
