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
    def reacnetgen_param(self, request):
        return request.param
        
    @pytest.fixture()
    def reacnetgen(self, reacnetgen_param):
        rngclass = ReacNetGenerator(**reacnetgen_param['rngparams'])
        yield rngclass
        assert checksha256(rngclass.reactionfilename,
                           reacnetgen_param['reaction_sha256'])
        assert os.path.exists(rngclass.reactionfilename)
        assert os.path.exists(rngclass.resultfilename)

    def test_reacnetgen(self, reacnetgen):
        """Test main process of ReacNetGen."""
        reacnetgen.runanddraw()

    def test_reacnetgen_step(self, reacnetgen):
        reacnetgen.run()
        reacnetgen.draw()
        reacnetgen.report()

    @pytest.fixture()
    def reacnetgengui(self):
        try:
            gui = GUI()
            yield gui
            gui.quit()
        except TclError:
            pytest.skip("No display for GUI.")

    def test_gui(self, reacnetgengui):
        """Test GUI of ReacNetGen."""
        reacnetgengui.root.after(100, reacnetgengui.root.destroy)
        reacnetgengui.gui()
    
    def test_gui_openandrun(self, reacnetgengui, mocker, reacnetgen_param):
        mocker.patch("tkinter.filedialog.askopenfilename", return_value="dump.reaxc")
        pp = reacnetgen_param['rngparams']
        download_file(pp['urls'][0]['url'][0], pp['urls'][0]['fn'], None)
        reacnetgengui._atomnameet.delete(0, END)
        reacnetgengui._atomnameet.insert(0, " ".join(pp['atomname']))
        reacnetgengui._filetype.set(pp['inputfiletype'])
        reacnetgengui._openbtn.invoke()
        reacnetgengui._runbtn.invoke()

    def test_commandline_help(self, script_runner):
        ret = script_runner.run('reacnetgenerator', '-h')
        assert ret.success

    def test_commandline_run(self, script_runner, reacnetgen_param):
        pp = reacnetgen_param['rngparams']
        cc_hmm = '' if pp['runHMM'] else '--nohmm'
        cc_dump = '--dump' if pp['inputfiletype'] == 'lammpsdumpfile' else ''
        cc_atomname = ' '.join(pp['atomname'])
        ret = script_runner.run('reacnetgenerator', '-i', pp['inputfilename'], '-a', cc_atomname, cc_dump,
                                '-s', pp['atomname'][0], cc_hmm, '--urls', pp['urls'][0]['fn'], pp['urls'][0]['url'][0])
        assert ret.success

    @pytest.mark.xfail
    def test_unsupported_filetype(self):
        ReacNetGenerator(inputfilename="xx", inputfiletype="abc", atomname=[])
