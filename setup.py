from setuptools import setup
setup(name='ReacNetGenerator',
      version='1.2.10',
      description='Reaction Network Generator',
      keywords="reaction network",
      url='https://github.com/njzjz/ReacNetGenerator',
      author='Jinzhe Zeng',
      author_email='njzjz@qq.com',
      packages=['ReacNetGenerator'],
      install_requires=['numpy','scipy','networkx','scikit-learn','matplotlib','hmmlearn','htmlmin'])
