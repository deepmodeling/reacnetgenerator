# cython: language_level=3
# cython: linetrace=True
"""Download trajectory before running."""

from .utils import SharedRNGData, download_multifiles

class DownloadData(SharedRNGData):
    def __init__(self, rng):
        SharedRNGData.__init__(self, rng, ["urls"], [])

    def download_files(self):
        download_multifiles(self.urls)
