from setuptools import setup
setup(name='reacnetgenerator',
      description='Reaction Network Generator',
      keywords="reaction network",
      url='https://njzjz.github.io/reacnetgenerator/',
      author='Jinzhe Zeng',
      author_email='njzjz@qq.com',
      packages=['reacnetgenerator'],
      python_requires='~=3.6.0',
      install_requires=['numpy', 'scipy>=0.20.1', 'networkx',
                        'scikit-learn', 'matplotlib', 'hmmlearn>=0.2.1', 'htmlmin', 'ase', 'scour'],
      entry_points={
          'console_scripts': ['reacnetgenerator=reacnetgenerator.reacnetgen:_commandline', 'reacnetgeneratorgui=reacnetgenerator.gui:gui']
      },
      test_suite='reacnetgenerator.test',
      tests_require=['requests'],
      use_scm_version=True,
      setup_requires=['setuptools_scm'],
      )
