# Copyright 2018, East China Normal University
"""ReacNetGen"""


import logging
import coloredlogs
import matplotlib as mpl
mpl.use("svg")

from .reacnetgen import ReacNetGenerator, __version__

coloredlogs.install(
    fmt=f'%(asctime)s - ReacNetGen {__version__} - %(levelname)s: %(message)s', level=logging.INFO, milliseconds=True)
