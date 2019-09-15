# cython: language_level=3
"""Test ReacNetGen."""


import json
import os
import pytest
from tkinter import TclError

import pkg_resources
from reacnetgenerator import ReacNetGenerator
from reacnetgenerator.gui import GUI
from reacnetgenerator.utils import checksha256


class TestReacNetGen:
    """Test ReacNetGenerator."""

    @pytest.fixture(autouse=True)
    def chdir(self, tmp_path):
        start_direcroty = os.getcwd()
        os.chdir(tmp_path)
        yield
        os.chdir(start_direcroty)

    @pytest.fixture(params=json.load(pkg_resources.resource_stream(__name__, 'test.json')))
    def reacnetgen(self, request, tmp_path):
        testparm = request.param
        rngclass = ReacNetGenerator(**testparm['rngparams'])
        yield rngclass
        assert checksha256(rngclass.reactionfilename,
                           testparm['reaction_sha256'])
        assert os.path.exists(rngclass.reactionfilename)
        assert os.path.exists(rngclass.resultfilename)

    def test_reacnetgen(self, reacnetgen):
        """Test main process of ReacNetGen."""
        reacnetgen.runanddraw()

    def test_reacnetgen_step(self, reacnetgen):
        reacnetgen.run()
        reacnetgen.draw()
        reacnetgen.report()

    def test_gui(self):
        """Test GUI of ReacNetGen."""
        try:
            gui = GUI()
            gui.root.after(1000, gui.root.destroy)
            gui.gui()
        except TclError:
            pytest.skip("No display for GUI.")

    def test_commandline_help(self, script_runner):
        ret = script_runner.run('reacnetgenerator', '-h')
        assert ret.success

    def test_commandline_run(self, script_runner):
        ret = script_runner.run('reacnetgenerator', '-i', 'dump.reaxc', '-a', 'C', 'H', 'O', '--dump', '-s', 'C', '--nohmm',
                                '--urls', 'dump.reaxc', 'https://drive.google.com/uc?authuser=0&id=1-MZZEpTj71JJn4JfKPh5yb_lD2V7NS-Y&export=download')
        assert ret.success

    @pytest.mark.xfail
    def test_unsupported_filetype(self):
        ReacNetGenerator(inputfilename="xx", inputfiletype="abc", atomname=[])
