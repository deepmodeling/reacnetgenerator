from setuptools import setup
setup(name='reacnetgenerator',
      version='1.2.20',
      description='Reaction Network Generator',
      keywords="reaction network",
      url='https://njzjz.github.io/reacnetgenerator/',
      author='Jinzhe Zeng',
      author_email='njzjz@qq.com',
      packages=['reacnetgenerator'],
      install_requires=['numpy', 'scipy>=0.20.1', 'networkx',
                        'scikit-learn', 'matplotlib', 'hmmlearn>=0.2.1', 'htmlmin', 'ase', 'scour'],
      entry_points={
          'console_scripts': ['reacnetgenerator=reacnetgenerator.reacnetgen:_commandline', 'reacnetgeneratorgui=reacnetgenerator.gui:gui']
      },
      test_suite='reacnetgenerator.test',
      tests_require=['requests'],
      )
