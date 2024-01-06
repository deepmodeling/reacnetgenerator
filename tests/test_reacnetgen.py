# SPDX-License-Identifier: LGPL-3.0-or-later
# cython: language_level=3
"""Test ReacNetGen."""


import fileinput
import itertools
import json
import os
from tkinter import END, TclError

import numpy as np
import pytest

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._detect import _Detect
from reacnetgenerator._hmmfilter import _HMMFilter
from reacnetgenerator._path import _CollectSMILESPaths
from reacnetgenerator.commandline import parm2cmd
from reacnetgenerator.gui import GUI
from reacnetgenerator.utils import checksha256, download_multifiles, listtobytes

with open(os.path.join(os.path.dirname(__file__), "test.json")) as f:
    test_data = json.load(f)


class TestReacNetGen:
    """Test ReacNetGenerator."""

    @pytest.fixture(autouse=True)
    def chdir(self, tmp_path):
        """Change directory to tmp_path."""
        start_direcroty = os.getcwd()
        os.chdir(tmp_path)
        yield
        os.chdir(start_direcroty)

    @pytest.fixture(
        params=[
            pytest.param(
                param, marks=(pytest.mark.xfail if param.get("xfail", False) else ())
            )
            for param in test_data
        ]
    )
    def reacnetgen_param(self, request):
        """Fixture for ReacNetGenerator parameters."""
        return request.param

    @pytest.fixture()
    def reacnetgen(self, reacnetgen_param):
        """Fixture for ReacNetGenerator."""
        rngclass = ReacNetGenerator(**reacnetgen_param["rngparams"])
        yield rngclass
        assert checksha256(
            rngclass.reactionfilename, reacnetgen_param.get("reaction_sha256", [])
        )
        assert os.path.exists(rngclass.reactionfilename)
        assert os.path.exists(rngclass.resultfilename)

    def test_reacnetgen(self, reacnetgen):
        """Test main process of ReacNetGen."""
        reacnetgen.runanddraw()

    def test_reacnetgen_step(self, reacnetgen):
        """Test main process of ReacNetGen."""
        reacnetgen.run()
        reacnetgen.draw()
        reacnetgen.report()

    @pytest.fixture()
    def reacnetgengui(self):
        """Fixture for GUI test."""
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
        """Test GUI of ReacNetGenerator."""
        pp = reacnetgen_param["rngparams"]
        mocker.patch(
            "tkinter.filedialog.askopenfilename", return_value=pp["inputfilename"]
        )
        mocker.patch("tkinter.messagebox.showerror")
        download_multifiles(pp.get("urls", []))
        reacnetgengui._atomnameet.delete(0, END)
        reacnetgengui._atomnameet.insert(0, " ".join(pp["atomname"]))
        if pp["inputfiletype"] in ["lammpsbondfile", "lammpsdumpfile"]:
            reacnetgengui._filetype.set(pp["inputfiletype"])
        reacnetgengui._openbtn.invoke()
        reacnetgengui._runbtn.invoke()

    def test_commandline_help(self, script_runner):
        """Test commandline of ReacNetGenerator."""
        ret = script_runner.run(["reacnetgenerator", "-h"])
        assert ret.success

    def test_commandline_version(self, script_runner):
        """Test commandline of ReacNetGenerator."""
        ret = script_runner.run(["reacnetgenerator", "--version"])
        assert ret.success

    def test_commandline_run(self, script_runner, reacnetgen_param):
        """Test commandline of ReacNetGenerator."""
        ret = script_runner.run(parm2cmd(reacnetgen_param["rngparams"]))
        assert ret.success

    def test_benchmark_detect(self, benchmark, reacnetgen_param):
        """Benchmark _readstepfunc."""
        reacnetgen = ReacNetGenerator(**reacnetgen_param["rngparams"])
        reacnetgen._process((reacnetgen.Status.DOWNLOAD,))
        detectclass = _Detect.gettype(reacnetgen)
        with fileinput.input(files=detectclass.inputfilename) as f:
            nlines = detectclass._readNfunc(f)
        with fileinput.input(files=detectclass.inputfilename) as f:
            lines = next(itertools.zip_longest(*[f] * nlines))
        benchmark(detectclass._readstepfunc, (0, lines))

    def test_benchmark_hmm(self, benchmark, reacnetgen_param):
        """Benchmark _getoriginandhmm."""
        reacnetgen = ReacNetGenerator(**reacnetgen_param["rngparams"])
        reacnetgen.step = 250000
        hmmclass = _HMMFilter(reacnetgen)
        if hmmclass.runHMM:
            hmmclass._initHMM()
        rng = np.random.default_rng()
        index = np.sort(rng.choice(hmmclass.step, hmmclass.step // 2, replace=False))
        compressed_bytes = [
            listtobytes((5, 6)),
            listtobytes(((5, 6, 1),)),
            listtobytes(index),
        ]
        benchmark(hmmclass._getoriginandhmm, compressed_bytes)

    def test_re(self, reacnetgen_param):
        """Test regular expression of _HTMLResult."""
        reacnetgen = ReacNetGenerator(**reacnetgen_param["rngparams"])
        r = _CollectSMILESPaths(reacnetgen)
        r.atomname = ["C", "H", "O", "Na", "Cl"]
        assert r._re("C"), "[C]"
        assert r._re("[C]"), "[C]"
        assert r._re("[CH]"), "[CH]"
        assert r._re("Na"), "[Na]"
        assert r._re("[H]c(Cl)C([H])Cl"), "[H][c]([Cl])[C]([H])[Cl]"
