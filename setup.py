"""Welcome to install ReacNetGenerator.

Just use the following command to install:
$ pip install .
Note you should install OpenBabel and RDkit first:
$ conda create -q -n reacnetgenerator python=3.7 openbabel rdkit -c openbabel -c conda-forge
$ source activate reacnetgenerator
"""
from os import path

from setuptools import setup, find_packages, Extension


if __name__ == '__main__':
    print(__doc__)
    this_directory = path.abspath(path.dirname(__file__))
    with open(path.join(this_directory, 'docs', 'README.md')) as f:
        long_description = f.read()

    extensions = [
        Extension("reacnetgenerator.dps", ["reacnetgenerator/dps.pyx"]),
    ]

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
              'jinja2', 'coloredlogs', 'jsmin', 'cssmin',
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
              'reacnetgenerator': ['static/html/*.html', 'static/js/*.js',
                                   'static/css/*.css', 'static/img/*.png',
                                   'static/js/vendor/*.js',
                                   'static/css/vendor/*.css',
                                   'test/test.json',
                                   '*.pyx', '*.c', '*.cpp',
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
            Extension("reacnetgenerator.dps", sources=["reacnetgenerator/dps.pyx"]),
            ],
          )
