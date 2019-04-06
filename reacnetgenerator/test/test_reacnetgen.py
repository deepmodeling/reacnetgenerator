"""Test ReacNetGen."""


import hashlib
import json
import logging
import math
import os
import tempfile
import pytest
from tkinter import TclError

import pkg_resources
import reacnetgenerator
import reacnetgenerator.gui
import requests
from tqdm import tqdm


class TestReacNetGen:
    """Test ReacNetGenerator."""

    @pytest.fixture(params=json.load(
        pkg_resources.resource_stream(
            __name__, 'test.json')))
    def reacnetgen(self, request):
        folder = tempfile.mkdtemp(prefix='testfiles-', dir='.')
        logging.info(f'Folder: {folder}:')
        os.chdir(folder)

        testparm = request.param
        self._download_file(
            testparm['url'], testparm['filename'], testparm['sha256'])

        return reacnetgenerator.ReacNetGenerator(
            inputfilename=testparm['filename'], atomname=testparm['atomname'],
            SMILES=testparm['smiles'],
            inputfiletype=testparm['inputfiletype'],
            runHMM=testparm['hmm'],
            speciescenter=testparm['speciescenter']
            if 'speciescenter' in testparm else None,
            split=testparm['split'] if 'split' in testparm else 1,
            )

    def test_reacnetgen(self, reacnetgen):
        """Test main process of ReacNetGen."""
        reacnetgen.runanddraw()

        logging.info("Here are reactions:")
        with open(reacnetgen.reactionfilename) as f:
            print(f.read())
        assert os.path.exists(reacnetgen.resultfilename)

    def test_gui(self):
        """Test GUI of ReacNetGen."""
        try:
            gui = reacnetgenerator.gui.GUI()
            gui.root.after(1000, gui.root.destroy)
            gui.gui()
        except TclError:
            logging.warning("No display for GUI.")

    def _download_file(self, urls, pathfilename, sha256):
        times = 0
        # download if not exists
        while times < 3:
            if os.path.isfile(pathfilename) and self._checksha256(
                    pathfilename, sha256):
                break

            # from https://stackoverflow.com/questions/16694907
            if not isinstance(urls, list):
                urls = [urls]
            for url in urls:
                try:
                    logging.info(f"Try to download {pathfilename} from {url}")
                    r = requests.get(url, stream=True)
                    break
                except requests.exceptions.RequestException as e:
                    logging.warning(e)
                    logging.warning("Request Error.")
            else:
                logging.error(f"Cannot download {pathfilename}.")
                raise IOError(f"Cannot download {pathfilename}.")

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
        else:
            logging.error(f"Retry too much times.")
            raise IOError(f"Retry too much times.")
        return pathfilename

    @staticmethod
    def _checksha256(filename, sha256_check):
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
        if sha256 == sha256_check:
            return True
        logging.warning("SHA256 is not correct.")
        return False
