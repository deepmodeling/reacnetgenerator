# cython: language_level=3
# cython: linetrace=True
"""Init logging."""
import logging
import coloredlogs
from ._version import __version__
coloredlogs.install(
    fmt=f'%(asctime)s - ReacNetGen {__version__} - %(levelname)s: %(message)s',
    level=logging.INFO, milliseconds=True)
