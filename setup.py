"""Welcome to install ReacNetGenerator.

Just use `pip install .` to install.
Note you should install Yarn, OpenBabel, and RDkit first:
conda install python=3 yarn openbabel rdkit compilers -c conda-forge
"""
import subprocess as sp
import os
import shutil
import fnmatch

from distutils import log
from distutils.file_util import copy_file
from setuptools import setup, find_packages, Extension
import setuptools.command.build_py
import setuptools.command.build_ext


class BuildExtCommand(setuptools.command.build_ext.build_ext):

    def run(self):
        log.info(__doc__)
        log.info('Prepare JavaScript files with webpack...')
        yarn = shutil.which('yarn')
        if yarn is None:
            raise RuntimeError(
                "Yarn is not installed. Plase install it by `conda install yarn`.")
        try:
            sp.run(yarn, check=True, cwd=os.path.join(
                this_directory, 'reacnetgenerator', 'static', 'webpack'))
            sp.run([yarn, 'start'], check=True, cwd=os.path.join(
                this_directory, 'reacnetgenerator', 'static', 'webpack'))
            assert os.path.exists(os.path.join(
                this_directory, 'reacnetgenerator', 'static', 'webpack', 'bundle.html'))
        except (sp.CalledProcessError, AssertionError):
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
    with open(os.path.join(this_directory, 'README.md'), encoding="utf8") as f:
        return f.read()


if __name__ == '__main__':
    this_directory = os.path.abspath(os.path.dirname(__file__))
    encrypted_python_files = [
        "reacnetgenerator._detect",
        "reacnetgenerator._draw",
        "reacnetgenerator._hmmfilter",
        "reacnetgenerator._logging",
        "reacnetgenerator._matrix",
        "reacnetgenerator._path",
        "reacnetgenerator._reachtml",
        "reacnetgenerator._version",
        "reacnetgenerator.gui",
        "reacnetgenerator.reacnetgen",
        "reacnetgenerator._reaction",
        "reacnetgenerator.commandline",
        "reacnetgenerator.utils",
        "reacnetgenerator._download",
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

    tests_require = ['pytest-sugar', 'pytest-cov', 'cython',
                     'pytest-xvfb', "codecov>=1.4.0", "pytest-console-scripts",
                     "pytest-mock", "pytest-benchmark",
                    ],
    setup(name='reacnetgenerator',
          description='Reaction Network Generator',
          keywords="reaction network",
          url='https://njzjz.github.io/reacnetgenerator/',
          author='Jinzhe Zeng',
          author_email='jzzeng@stu.ecnu.edu.cn',
          packages=find_packages(),
          python_requires='~=3.6',
          install_requires=[
              'numpy>=1.15', 'scipy>=0.20.1', 'networkx',
              'matplotlib', 'hmmlearn>=0.2.1',
              'ase', 'scour', 'tqdm',
              'coloredlogs',
              'pandas', 'pybase64', 'lz4',
              'requests',
          ],
          entry_points={'console_scripts': [
              'reacnetgenerator=reacnetgenerator.commandline:_commandline',
              'reacnetgeneratorgui=reacnetgenerator.gui:gui'
          ]
          },
          test_suite='reacnetgenerator.test',
          tests_require=tests_require,
          extras_require={
              "test": tests_require,
              "docs": ['sphinx-markdown-builder'],
          },
          use_scm_version=os.path.exists(os.path.join(this_directory, ".git")),
          setup_requires=[
              'setuptools>=18.0',
              'setuptools_scm',
              'pytest-runner',
              'cython>=0.16',
              'numpy>=1.15',
          ],
          package_data={
              'reacnetgenerator': ['static/webpack/bundle.html',
                                   'static/img-title.png',
                                   'test/test.json',
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
          version="1.0.0",
          )
