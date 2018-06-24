from setuptools import setup
setup(name='ReacNetGenerator',
      version='1.1.8',
      description='Reaction Network Generator',
      keywords="reaction network",
      url='https://github.com/njzjz/ReacNetGenerator',
      author='Jin-Zhe Zeng',
      author_email='njzjz@qq.com',
      packages=['ReacNetGenerator'],
      install_requires=['numpy','scipy','networkx','sklearn','matplotlib','hmmlearn'])
