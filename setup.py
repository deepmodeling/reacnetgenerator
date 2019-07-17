"""Welcome to install ReacNetGenerator.

Just use `pip install .` to install.
Note you should install Yarn, OpenBabel, and RDkit first:
conda install python=3 yarn openbabel rdkit compilers -c conda-forge
"""
import subprocess as sp
from os import path
import shutil
import fnmatch

from setuptools import setup, find_packages, Extension
import setuptools.command.build_py


class BuildCommand(setuptools.command.build_py.build_py):

    def run(self):
        try:
            print(__doc__)
            print('Prepare JavaScript files with webpack...')
            yarn = shutil.which('yarn')
            sp.run(yarn, check=True, cwd=path.join(
                this_directory, 'reacnetgenerator', 'static', 'webpack'))
            sp.run([yarn, 'start'], check=True, cwd=path.join(
                this_directory, 'reacnetgenerator', 'static', 'webpack'))
            assert path.exists(path.join(
                this_directory, 'reacnetgenerator', 'static', 'webpack', 'bundle.js'))
        except sp.CalledProcessError:
            raise ImportError(
                "Maybe you didn't install yarn? Plase install it by `conda install yarn`.")
        except AssertionError:
            raise OSError("No bundle.js found, please retry.")
        setuptools.command.build_py.build_py.run(self)

    def find_package_modules(self, package, package_dir):
        modules = super().find_package_modules(package, package_dir)
        return [(pkg, mod, file, ) for (pkg, mod, file, ) in modules
                if not any(fnmatch.fnmatchcase(pkg + '.' + mod, pat=pattern)
                           for pattern in encrypted_python_files)]


def readme():
    with open(path.join(this_directory, 'README.md'), encoding="utf8") as f:
        return f.read()


if __name__ == '__main__':
    this_directory = path.abspath(path.dirname(__file__))
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
    ]

    ext_modules = [
        Extension("reacnetgenerator.dps", sources=[
            "reacnetgenerator/dps.pyx", "reacnetgenerator/c_stack.cpp"], language="c++", define_macros=[('CYTHON_TRACE', '1')]),
    ]
    # encrypt python files
    ext_modules.extend([Extension(encrypted_python_file, sources=[
                       f"{path.join(*encrypted_python_file.split('.'))}{path.extsep}py"],
                       language="c", define_macros=[('CYTHON_TRACE', '1')]) for encrypted_python_file in encrypted_python_files])

    tests_require = ['requests', 'pytest-sugar', 'pytest-cov', 'cython'],
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
              'scikit-learn', 'matplotlib', 'hmmlearn>=0.2.1',
              'htmlmin', 'ase', 'scour', 'tqdm',
              'jinja2', 'coloredlogs',
              'pandas', 'pybase64', 'lz4'
          ],
          entry_points={'console_scripts': [
              'reacnetgenerator=reacnetgenerator.reacnetgen:_commandline',
              'reacnetgeneratorgui=reacnetgenerator.gui:gui'
          ]
          },
          test_suite='reacnetgenerator.test',
          tests_require=tests_require,
          extras_require={
              "test": tests_require,
          },
          use_scm_version=path.exists(path.join(this_directory, ".git")),
          setup_requires=[
              'setuptools>=18.0',
              'setuptools_scm',
              'pytest-runner',
              'cython',
          ],
          package_data={
              'reacnetgenerator': ['static/webpack/bundle.html',
                                   'static/img-title.png',
                                   'test/test.json',
                                   'tox.ini'
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
              "Programming Language :: JavaScript",
              "Topic :: Scientific/Engineering :: Chemistry",
              "Topic :: Scientific/Engineering :: Visualization",
              "Topic :: Software Development :: Libraries :: Python Modules",
              "Topic :: Software Development :: Version Control :: Git",
          ],
          zip_safe=True,
          ext_modules=ext_modules,
          cmdclass={"build_py": BuildCommand},
          version="1.0.0",
          )
