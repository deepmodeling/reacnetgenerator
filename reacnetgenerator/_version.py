# SPDX-License-Identifier: LGPL-3.0-or-later
# cython: language_level=3
# cython: linetrace=True
"""Obtain the version."""
try:
    from ._version2 import version as __version__
except ImportError:
    __version__ = "Unknown"
