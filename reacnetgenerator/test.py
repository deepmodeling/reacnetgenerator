"""Test ReacNetGen."""


import hashlib
import json
import logging
import math
import os
import unittest

import pkg_resources
import reacnetgenerator
import reacnetgenerator.gui
import requests
from tqdm import tqdm


class TestReacNetGen(unittest.TestCase):
    """Test ReacNetGenerator."""

    def test_reacnetgen(self):
        """Test main process of ReacNetGen."""
        testparms = json.load(
            pkg_resources.resource_stream(__name__, 'test.json'))

        for testparm in testparms:
            pathfilename = os.path.join(
                testparm['folder'], testparm['filename'])

            self._download_file(
                testparm['url'], pathfilename, testparm['sha256'])

            r = reacnetgenerator.ReacNetGenerator(
                inputfilename=pathfilename, atomname=testparm['atomname'],
                SMILES=testparm['smiles'],
                inputfiletype=testparm['inputfiletype'],
                runHMM=testparm['hmm'],
                speciescenter=testparm['speciescenter']
                if 'speciescenter' in testparm else None)
            r.runanddraw()

            logging.info("Here are reactions:")
            with open(r.reactionfilename) as f:
                for line in f:
                    print(line.strip())
            self.assertTrue(os.path.exists(r.resultfilename))

    def _download_file(self, url, pathfilename, sha256):
        # download if not exists
        while not os.path.isfile(pathfilename) or self._checksha256(
                pathfilename) != sha256:
            try:
                os.makedirs(os.path.split(pathfilename)[0])
            except OSError:
                pass

            # from https://stackoverflow.com/questions/16694907/how-to-download-large-file-in-python-with-requests-py
            r = requests.get(url, stream=True)
            total_size = int(r.headers.get('content-length', 0))
            block_size = 1024
            with open(pathfilename, 'wb') as f:
                for chunk in tqdm(
                        r.iter_content(chunk_size=1024),
                        total=math.ceil(total_size // block_size),
                        unit='KB', unit_scale=True,
                        desc=f"Downloading {pathfilename}..."):
                    if chunk:
                        f.write(chunk)
        return pathfilename

    @staticmethod
    def _checksha256(filename):
        if not os.path.isfile(filename):
            return
        h = hashlib.sha256()
        b = bytearray(128*1024)
        mv = memoryview(b)
        with open(filename, 'rb', buffering=0) as f:
            for n in iter(lambda: f.readinto(mv), 0):
                h.update(mv[:n])
        sha256 = h.hexdigest()
        logging.info(f"SHA256 of {filename}: {sha256}")
        return sha256


if __name__ == '__main__':
    unittest.main()
