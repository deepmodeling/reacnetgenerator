# SPDX-License-Identifier: LGPL-3.0-or-later
"""Welcome to install ReacNetGenerator.

Just use `pip install .` to install.
"""
import logging
import os
import subprocess
from pathlib import Path

import setuptools.command.build_ext
import yaml
from setuptools import Extension, setup

log = logging.getLogger(__name__)


def run_node_command(args, cwd):
    """Call node with subprocess."""
    return subprocess.call(["node", *args], cwd=cwd)


class BuildExtCommand(setuptools.command.build_ext.build_ext):
    """A custom command to build the extension."""

    def run(self):
        """Run the command."""
        assert __doc__ is not None
        log.info(__doc__)
        log.info("Preparing JavaScript files with webpack...")

        this_directory = Path(__file__).parent
        webpack_dir = this_directory / "reacnetgenerator" / "static" / "webpack"

        node_call = run_node_command

        with open(webpack_dir / ".yarnrc.yml") as f:
            yarn_path = str(Path(yaml.load(f, Loader=yaml.Loader)["yarnPath"]))

        node_call([yarn_path], cwd=webpack_dir)
        node_call([yarn_path, "start"], cwd=webpack_dir)

        bundle_html_path = webpack_dir / "bundle.html"

        if not bundle_html_path.exists():
            raise RuntimeError("Failed to build bundle.html with Yarn, please retry.")

        # Copy files to build_lib
        build_lib_dir = os.path.join(
            self.build_lib, "reacnetgenerator", "static", "webpack"
        )
        os.makedirs(build_lib_dir, exist_ok=True)

        self.copy_file(
            str(bundle_html_path),
            os.path.join(build_lib_dir, "bundle.html"),
        )

        # Add numpy headers to include_dirs
        import numpy as np

        self.include_dirs.append(np.get_include())
        super().run()


if __name__ == "__main__":
    define_macros = []
    if os.environ.get("DEBUG", 0):
        define_macros.extend((("CYTHON_TRACE", "1"), ("CYTHON_TRACE_NOGIL", "1")))

    ext_modules = [
        Extension(
            "reacnetgenerator.dps",
            sources=["reacnetgenerator/dps.pyx", "reacnetgenerator/c_stack.cpp"],
            language="c++",
            define_macros=define_macros,
        ),
        Extension(
            "reacnetgenerator.utils_np",
            sources=["reacnetgenerator/utils_np.pyx"],
            language="c",
            define_macros=define_macros,
        ),
    ]

    setup(
        ext_modules=ext_modules,
        cmdclass={
            "build_ext": BuildExtCommand,
        },
    )
