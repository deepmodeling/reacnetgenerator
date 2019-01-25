from setuptools import setup
from os import path

if __name__=='__main__':
    this_directory = path.abspath(path.dirname(__file__))
    with open(path.join(this_directory, 'docs', 'README.md'), encoding='utf-8') as f:
        long_description = f.read()

    setup(name='reacnetgenerator',
          description='Reaction Network Generator',
          keywords="reaction network",
          url='https://njzjz.github.io/reacnetgenerator/',
          author='Jinzhe Zeng',
          author_email='njzjz@qq.com',
          packages=['reacnetgenerator'],
          python_requires='~=3.6',
          install_requires=['numpy', 'scipy>=0.20.1', 'networkx',
                            'scikit-learn', 'matplotlib', 'hmmlearn>=0.2.1',
                            'htmlmin', 'ase', 'scour', 'tqdm',
                            'jinja2',
                            ],
          entry_points={
              'console_scripts': ['reacnetgenerator=reacnetgenerator.reacnetgen:_commandline',
                                  'reacnetgeneratorgui=reacnetgenerator.gui:gui']
          },
          test_suite='reacnetgenerator.test',
          tests_require=['requests'],
          use_scm_version=True,
          setup_requires=['setuptools_scm'],
          package_data={
              'reacnetgenerator': ['static/html/*.html', 'static/js/*.js',
                                   'static/css/*.css', 'static/img/*.png', 'test.json'],
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
          )
