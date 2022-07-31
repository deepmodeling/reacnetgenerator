import os
from pathlib import Path

import numpy as np
import pytest
from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._detect import _Detect

p_inputs = Path(__file__).parent / 'inputs'


class TestDetect:
    """Test different detect format.
    
    All systems contain a single water molecule: H, H, O.
    """

    @pytest.fixture(autouse=True)
    def chdir(self, tmp_path):
        start_direcroty = os.getcwd()
        os.chdir(tmp_path)
        yield
        os.chdir(start_direcroty)

    @pytest.fixture(params=[
        # inputfiletype, inputfilename
        ("lammpsdumpfile", p_inputs / "water.dump"),
        ("lammpsbondfile", p_inputs / "water.bond"),
        ("xyz", p_inputs / "water.xyz"),
    ])
    def reacnetgen_param(self, request):
        return request.param

    @pytest.fixture()
    def reacnetgen(self, reacnetgen_param):
        rngclass = ReacNetGenerator(
            inputfiletype=reacnetgen_param[0],
            inputfilename=reacnetgen_param[1],
            atomname=["H", "O"],
            pbc=False,
        )
        yield rngclass

    def test_reacnetgen(self, reacnetgen):
        """Test main process of ReacNetGen."""
        _Detect.gettype(reacnetgen).detect()
        assert reacnetgen.N == 3
        np.testing.assert_array_equal(reacnetgen.atomtype, np.array([0, 0, 1], dtype=int))
