"""Download trajectory before running."""

from .utils import SharedRNGData, download_file

class DownloadData(SharedRNGData):
    def __init__(self, rng):
        SharedRNGData.__init__(self, rng, ["urls"], [])
    
    def download_files(self):
        for jdata in self.urls:
            download_file(jdata["url"], jdata["fn"], jdata.get("sha256", None))