# cython: language_level=3
# cython: linetrace=True
"""Obtain the version."""
from pkg_resources import DistributionNotFound, get_distribution
try:
    from ._version2 import version as __version__
except ImportError:
    try:
        __version__ = get_distribution('reacnetgenerator').version
    except DistributionNotFound:
        __version__ = 'Unknown'
