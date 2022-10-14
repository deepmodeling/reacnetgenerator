"""Welcome to install ReacNetGenerator.

Just use `pip install .` to install.
"""
import os
from pathlib import Path
import fnmatch

import yaml
from nodejs import node
from distutils import log
from distutils.file_util import copy_file
from setuptools import setup, find_packages, Extension
import setuptools.command.build_py
import setuptools.command.build_ext


class BuildExtCommand(setuptools.command.build_ext.build_ext):

    def run(self):
        log.info(__doc__)
        log.info('Prepare JavaScript files with webpack...')
        this_directory = Path(__file__).parent
        webpack_dir = this_directory / "reacnetgenerator" / "static" / "webpack"
        with open(webpack_dir / ".yarnrc.yml") as f:
            yarn_path = Path(yaml.load(f, Loader=yaml.Loader)["yarnPath"])
        node.call([yarn_path], cwd=webpack_dir)
        node.call([yarn_path, "start"], cwd=webpack_dir)
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

class BuildPyCommand(setuptools.command.build_py.build_py):

    def find_package_modules(self, package, package_dir):
        modules = super().find_package_modules(package, package_dir)
        return [(pkg, mod, file, ) for (pkg, mod, file, ) in modules
                if not any(fnmatch.fnmatchcase(pkg + '.' + mod, pat=pattern)
                           for pattern in encrypted_python_files)]


def readme():
    this_directory = Path(__file__).parent
    with open(os.path.join(this_directory, 'README.md'), encoding="utf8") as f:
        return f.read()


if __name__ == '__main__':
    encrypted_python_files = [
    ]

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
    # encrypt python files
    ext_modules.extend([Extension(encrypted_python_file, sources=[
            os.path.join(*encrypted_python_file.split('.')) + os.path.extsep + "py"],
            language="c", define_macros=define_macros,
        ) for encrypted_python_file in encrypted_python_files])

    tests_require = ['pytest-sugar', 'pytest-cov<5', 'cython',
                     'pytest-xvfb', "codecov>=1.4.0", "pytest-console-scripts",
                     "pytest-mock", "pytest-benchmark",
                    ],
    no_conda_require = []
    if not os.environ.get("CONDA_BUILD", 0):
        no_conda_require.extend(["openbabel-wheel", "rdkit"])
    setup(name='reacnetgenerator',
          description='Reaction Network Generator',
          keywords="reaction network",
          url='https://njzjz.github.io/reacnetgenerator/',
          author='Jinzhe Zeng',
          author_email='jzzeng@stu.ecnu.edu.cn',
          packages=find_packages(),
          python_requires='~=3.6',
          install_requires=[
              'numpy', 'scipy>=0.20.1', 'networkx',
              'matplotlib', 'hmmlearn>=0.2.1',
              'ase', 'scour', 'tqdm',
              'coloredlogs',
              'pandas', 'pybase64', 'lz4',
              'requests',
          ] + no_conda_require,
          entry_points={'console_scripts': [
              'reacnetgenerator=reacnetgenerator.commandline:_commandline',
              'reacnetgeneratorgui=reacnetgenerator.gui:gui'
          ]
          },
          tests_require=tests_require,
          extras_require={
              "test": tests_require,
              "docs": [
                'sphinx',
                'sphinx-press-theme',
                'numpydoc',
                'sphinx-argparse',
                'myst_parser',
                'sphinx-favicon',
              ],
          },
          use_scm_version=True,
          setup_requires=[
              'setuptools>=18.0',
              'setuptools_scm',
              'pytest-runner',
              'cython>=0.16',
              'numpy',
          ],
          package_data={
              'reacnetgenerator': ['static/webpack/bundle.html',
                                   'static/img-title.png',
                                   ],
          },
          long_description=readme(),
          long_description_content_type='text/markdown',
          classifiers=[
              "Natural Language :: English",
              "Operating System :: POSIX :: Linux",
              "Operating System :: Microsoft :: Windows",
              "Programming Language :: Python :: 3.6",
              "Programming Language :: Python :: 3.7",
              "Programming Language :: Python :: 3.8",
              "Programming Language :: JavaScript",
              "Topic :: Scientific/Engineering :: Chemistry",
              "Topic :: Scientific/Engineering :: Visualization",
              "Topic :: Software Development :: Libraries :: Python Modules",
              "Topic :: Software Development :: Version Control :: Git",
          ],
          zip_safe=True,
          ext_modules=ext_modules,
          cmdclass={
              "build_py": BuildPyCommand,
              "build_ext": BuildExtCommand,
          },
          )
