import subprocess
import sys

import reacnetgenerator

def test_module_script():
    expected_version = reacnetgenerator.__version__
    output = subprocess.check_output([sys.executable, "-m", "reacnetgenerator", "--version"]).decode("ascii")
    assert output.splitlines()[0] == "ReacNetGenerator v%s" % expected_version