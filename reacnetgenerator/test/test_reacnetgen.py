# cython: language_level=3
"""Test ReacNetGen."""


import fileinput
import itertools
import json
import os
from tkinter import END, TclError

import numpy as np
import pkg_resources
import pytest

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._detect import _Detect
from reacnetgenerator._hmmfilter import _HMMFilter
from reacnetgenerator.gui import GUI
from reacnetgenerator.utils import checksha256, download_multifiles, listtobytes
from reacnetgenerator.commandline import parm2cmd


class TestReacNetGen:
    """Test ReacNetGenerator."""

    @pytest.fixture(autouse=True)
    def chdir(self, tmp_path):
        start_direcroty = os.getcwd()
        os.chdir(tmp_path)
        yield
        os.chdir(start_direcroty)

    @pytest.fixture(params=[pytest.param(param, marks=(pytest.mark.xfail if param.get("xfail", False) else ()))
                            for param in json.load(pkg_resources.resource_stream(__name__, 'test.json'))])
    def reacnetgen_param(self, request):
        return request.param

    @pytest.fixture()
    def reacnetgen(self, reacnetgen_param):
        rngclass = ReacNetGenerator(**reacnetgen_param['rngparams'])
        yield rngclass
        assert checksha256(rngclass.reactionfilename,
                           reacnetgen_param.get('reaction_sha256', []))
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
            try:
                gui.root.destroy()
            except TclError:
                pass
        except TclError:
            pytest.skip("No display for GUI.")

    def test_gui(self, reacnetgengui):
        """Test GUI of ReacNetGen."""
        reacnetgengui.root.after(100, reacnetgengui.root.destroy)
        reacnetgengui.gui()

    def test_gui_openandrun(self, reacnetgengui, mocker, reacnetgen_param):
        pp = reacnetgen_param['rngparams']
        mocker.patch("tkinter.filedialog.askopenfilename",
                     return_value=pp['inputfilename'])
        mocker.patch("tkinter.messagebox.showerror")
        download_multifiles(pp.get('urls', []))
        reacnetgengui._atomnameet.delete(0, END)
        reacnetgengui._atomnameet.insert(0, " ".join(pp['atomname']))
        if pp['inputfiletype'] in ['lammpsbondfile', 'lammpsdumpfile']:
            reacnetgengui._filetype.set(pp['inputfiletype'])
        reacnetgengui._openbtn.invoke()
        reacnetgengui._runbtn.invoke()

    def test_commandline_help(self, script_runner):
        ret = script_runner.run('reacnetgenerator', '-h')
        assert ret.success

    def test_commandline_run(self, script_runner, reacnetgen_param):
        ret = script_runner.run(*parm2cmd(reacnetgen_param['rngparams']))
        assert ret.success

    def test_benchmark_detect(self, benchmark, reacnetgen_param):
        reacnetgen = ReacNetGenerator(**reacnetgen_param['rngparams'])
        reacnetgen._process((reacnetgen.Status.DOWNLOAD,))
        detectclass = _Detect.gettype(reacnetgen)
        with fileinput.input(files=detectclass.inputfilename) as f:
            nlines = detectclass._readNfunc(f)
        with fileinput.input(files=detectclass.inputfilename) as f:
            lines = next(itertools.zip_longest(*[f] * nlines))
        benchmark(detectclass._readstepfunc, (0, lines))

    def test_benchmark_hmm(self, benchmark, reacnetgen_param):
        reacnetgen = ReacNetGenerator(**reacnetgen_param['rngparams'])
        reacnetgen.step = 250000
        hmmclass = _HMMFilter(reacnetgen)
        if hmmclass.runHMM:
            hmmclass._initHMM()
        index = np.sort(np.random.choice(
            hmmclass.step, hmmclass.step//2, replace=False))
        compressed_bytes = [listtobytes((5, 6)), listtobytes(
            ((5, 6, 1),)), listtobytes(index)]
        benchmark(hmmclass._getoriginandhmm, compressed_bytes)
