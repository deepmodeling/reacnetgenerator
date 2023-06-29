# SPDX-License-Identifier: LGPL-3.0-or-later
# cython: language_level=3
# cython: linetrace=True
"""Download trajectory before running."""

from typing import List

from .utils import SharedRNGData, download_multifiles


class DownloadData(SharedRNGData):
    urls: List[dict]

    def __init__(self, rng):
        SharedRNGData.__init__(self, rng, ["urls"], [])

    def download_files(self):
        download_multifiles(self.urls)
