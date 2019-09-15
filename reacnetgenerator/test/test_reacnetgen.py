# cython: language_level=3
"""Test ReacNetGen."""


import json
import os
import pytest
from tkinter import END, TclError

import pkg_resources
from reacnetgenerator import ReacNetGenerator
from reacnetgenerator.gui import GUI
from reacnetgenerator.utils import checksha256, download_file


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

    def reacnetgengui(self):
        try:
            gui = GUI()
            yield gui
            gui.quit()
        except TclError:
            pytest.skip("No display for GUI.")

    def test_gui(self, reacnetgengui):
        """Test GUI of ReacNetGen."""
        reacnetgengui.root.after(100, gui.root.destroy)
        reacnetgengui.gui()
    
    def test_gui_openandrun(self, reacnetgengui, mocker):
        mocker.patch("tkinter.filedialog.askopenfilename", return_value="dump.reaxc")
        download_file('https://drive.google.com/uc?authuser=0&id=1-MZZEpTj71JJn4JfKPh5yb_lD2V7NS-Y&export=download', 'dump.reaxc', None)
        reacnetgengui._atomnameet.delete(0, END)
        reacnetgengui._atomnameet.insert(0, "H O")
        reacnetgengui._openbtn.invoke()
        reacnetgengui._runbtn.invoke()

    def test_commandline_help(self, script_runner):
        ret = script_runner.run('reacnetgenerator', '-h')
        assert ret.success

    def test_commandline_run(self, script_runner):
        ret = script_runner.run('reacnetgenerator', '-i', 'dump.reaxc', '-a', 'H', 'O', '--dump', '-s', 'H', '--nohmm',
                                '--urls', 'dump.reaxc', 'https://drive.google.com/uc?authuser=0&id=1-MZZEpTj71JJn4JfKPh5yb_lD2V7NS-Y&export=download')
        assert ret.success

    @pytest.mark.xfail
    def test_unsupported_filetype(self):
        ReacNetGenerator(inputfilename="xx", inputfiletype="abc", atomname=[])
