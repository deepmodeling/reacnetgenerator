# Copyright 2018-2019, East China Normal University
"""ReacNetGen."""


import logging
import coloredlogs
import matplotlib as mpl
mpl.use("svg")  # noqa

from .reacnetgen import ReacNetGenerator, __version__

__all__ = ['ReacNetGenerator']

coloredlogs.install(
    fmt=f'%(asctime)s - ReacNetGen {__version__} - %(levelname)s: %(message)s',
    level=logging.INFO, milliseconds=True)
