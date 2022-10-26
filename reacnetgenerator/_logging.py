# cython: language_level=3
# cython: linetrace=True
"""This module provides a logger for logging."""
import logging
import coloredlogs
from ._version import __version__

logger = logging.getLogger(__name__)
coloredlogs.install(
    fmt=f'%(asctime)s - ReacNetGenerator {__version__} - %(levelname)s: %(message)s',
    level=logging.INFO, milliseconds=True, logger=logger)

__all__ = ['logger']
