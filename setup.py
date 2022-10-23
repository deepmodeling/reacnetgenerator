"""Welcome to install ReacNetGenerator.

Just use `pip install .` to install.
"""
import os
from pathlib import Path

import yaml
from distutils import log
from distutils.file_util import copy_file
from setuptools import setup, Extension
import setuptools.command.build_ext
try:
    from nodejs.node import call as node_call
except ModuleNotFoundError:
    log.info("nodejs-bin is not installed, try system node")
    import subprocess
    def node_call(args, **kwargs):
        return subprocess.call(['node', *args], **kwargs)


class BuildExtCommand(setuptools.command.build_ext.build_ext):

    def run(self):
        assert __doc__ is not None
        log.info(__doc__)
        log.info('Prepare JavaScript files with webpack...')
        this_directory = Path(__file__).parent
        webpack_dir = this_directory / "reacnetgenerator" / "static" / "webpack"
        with open(webpack_dir / ".yarnrc.yml") as f:
            yarn_path = str(Path(yaml.load(f, Loader=yaml.Loader)["yarnPath"]))
        node_call([yarn_path], cwd=webpack_dir)
        node_call([yarn_path, "start"], cwd=webpack_dir)
        try:
            assert (webpack_dir / "bundle.html").exists()
        except AssertionError:
            raise RuntimeError(
                "Fail to build bundle.html with Yarn, please retry.")
        # copy files
        try:
            os.makedirs(os.path.join(self.build_lib, 'reacnetgenerator', 'static', 'webpack'))
        except OSError:
            pass
        copy_file(
            os.path.join(this_directory, 'reacnetgenerator', 'static', 'webpack', 'bundle.html'),
            os.path.join(self.build_lib, 'reacnetgenerator', 'static', 'webpack', 'bundle.html'),
            verbose=self.verbose,
            dry_run=self.dry_run
        )
        # Add numpy headers to include_dirs
        import numpy as np
        self.include_dirs.append(np.get_include())
        setuptools.command.build_ext.build_ext.run(self)

if __name__ == '__main__':
    define_macros = []
    if os.environ.get("DEBUG", 0):
         define_macros.extend((('CYTHON_TRACE', '1'), ('CYTHON_TRACE_NOGIL', '1')))

    ext_modules = [
        Extension("reacnetgenerator.dps", sources=[
            "reacnetgenerator/dps.pyx", "reacnetgenerator/c_stack.cpp"],
            language="c++", define_macros=define_macros,
        ),
        Extension("reacnetgenerator.utils_np", sources=[
            "reacnetgenerator/utils_np.pyx"],
            language="c", define_macros=define_macros,
        ),
    ]

    setup(
          ext_modules=ext_modules,
          cmdclass={
              "build_ext": BuildExtCommand,
          },
          )
