# cython: language_level=3
"""Test ReacNetGen."""


import hashlib
import json
import logging
import math
import os
import tempfile
import pytest
import shutil
from tkinter import TclError

import pkg_resources
import reacnetgenerator
import reacnetgenerator.gui
import requests

class TestReacNetGen:
    """Test ReacNetGenerator."""

    @pytest.fixture(params=json.load(
        pkg_resources.resource_stream(
            __name__, 'test.json')))
    def reacnetgen(self, request, tmp_path):
        os.chdir(tmp_path)

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
            ), testparm

    def test_reacnetgen(self, reacnetgen):
        """Test main process of ReacNetGen."""
        rngclass, testparm = reacnetgen
        rngclass.runanddraw()

        assert self._checksha256(reacnetgen.reactionfilename, testparm['reaction_sha256'])
        assert os.path.exists(reacnetgen.reactionfilename)
        assert os.path.exists(reacnetgen.resultfilename)
        

    def test_gui(self):
        """Test GUI of ReacNetGen."""
        try:
            gui = reacnetgenerator.gui.GUI()
            gui.root.after(1000, gui.root.destroy)
            gui.gui()
        except TclError:
            pytest.skip("No display for GUI.")

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
                    with requests.get(url, stream=True) as r, open(pathfilename, 'wb') as f:
                        shutil.copyfileobj(r.raw, f)
                    break
                except requests.exceptions.RequestException as e:
                    logging.warning("Request Error.", exc_info=e)
            else:
                raise IOError(f"Cannot download {pathfilename}.")                
        else:
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
