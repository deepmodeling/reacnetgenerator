"""Welcome to install ReacNetGenerator.

Just use `pip install .` to install.
Note you should install Yarn, OpenBabel and RDkit first:
conda install python=3 yarn openbabel rdkit -c openbabel -c conda-forge
"""
import subprocess as sp
from os import path
import shutil

from setuptools import setup, find_packages, Extension
import setuptools.command.build_py


class BuildCommand(setuptools.command.build_py.build_py):

    def run(self):
        try:
            print('Prepare JavaScript files with webpack...')
            yarn = shutil.which('yarn')
            sp.run(yarn, check=True, cwd=path.join(
                this_directory, 'reacnetgenerator', 'static', 'webpack'))
            sp.run([yarn, 'start'], check=True, cwd=path.join(
                this_directory, 'reacnetgenerator', 'static', 'webpack'))
            assert path.exist(path.join(this_directory, 'reacnetgenerator', 'static', 'webpack', 'bundle.js'))
        except sp.CalledProcessError:
            raise ImportError(
                "Maybe you didn't install yarn? Plase install it by `conda install yarn`.")
        except AssertionError:
            raise OSError("No bundle.js found, please retry.")
        setuptools.command.build_py.build_py.run(self)


if __name__ == '__main__':
    print(__doc__)
    this_directory = path.abspath(path.dirname(__file__))
    with open(path.join(this_directory, 'docs', 'README.md')) as f:
        long_description = f.read()

    tests_require = ['requests', 'pytest-sugar', 'pytest-cov'],
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
          use_scm_version=True,
          setup_requires=[
              'setuptools>=18.0',
              'setuptools_scm',
              'pytest-runner',
              'cython',
          ],
          package_data={
              'reacnetgenerator': ['static/template.html',
                                   'static/webpack/bundle.js',
                                   'static/img-title.png',
                                   'test/test.json',
                                   ],
          },
          long_description=long_description,
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
          ext_modules=[
              Extension("reacnetgenerator.dps", sources=[
                        "reacnetgenerator/dps.pyx", "reacnetgenerator/c_stack.cpp"], language="c++"),
          ],
          cmdclass={"build_py": BuildCommand},
          )
