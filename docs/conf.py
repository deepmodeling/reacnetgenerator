# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
#import os
#import sys
#sys.path.insert(0, os.path.abspath('..'))
from datetime import datetime


# -- Project information -----------------------------------------------------

project = 'ReacNetGenerator'
copyright = '2019-%d, East China Normal University' % datetime.now().year
author = 'Jinzhe Zeng'

# The full version, including alpha/beta/rc tags
#release = 'reacnetgenerator'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinxarg.ext',
    'myst_parser',
    'numpydoc',
]

# Add any paths that contain templates here, relative to this directory.
#templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'press'
#html_logo = '_static/reacnetgen.svg'
html_static_path = ['_static']
html_js_files = ['https://unpkg.com/bilitube@0/dist/bilitube.min.js']
html_css_files = ['css/custom.css']
html_extra_path = ['report.html', 'fire.png', 'bundle.js', 'bundle.css']

html_theme_options = {
  "external_links": [
      ("Paper", "https://doi.org/10.1039/C9CP05091D"),
      ("Group", "https://computchem.cn"),
      ("GitHub", "https://github.com/tongzhugroup/reacnetgenerator"),
  ]
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']

intersphinx_mapping = {
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "python": ("https://docs.python.org/", None),
    "ase": ("https://wiki.fysik.dtu.dk/ase/", None),
}

def run_apidoc(_):
    from sphinx.ext.apidoc import main
    import os
    import sys
    sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    module = os.path.join(cur_dir, "..", "reacnetgenerator")
    main(['-M', '--tocfile', 'api', '-H', 'Python API', '-o', os.path.join(cur_dir, "api"), module, '--force'])

def copy_report(_):
    import subprocess as sp
    import os
    import shutil

    cur_dir = os.path.abspath(os.path.dirname(__file__))
    webpack = os.path.join(cur_dir, "..", "reacnetgenerator", "static", "webpack")
    sp.check_output(["yarn"], cwd=webpack)
    sp.check_output(["yarn", "start"], cwd=webpack, env={**os.environ, "REACNETGENERATOR_BUILDWEB": "1"})
    shutil.copy(os.path.join(webpack, "bundle.html"), os.path.join(cur_dir, "report.html"))
    shutil.copy(os.path.join(webpack, "bundle.js"), os.path.join(cur_dir, "bundle.js"))
    shutil.copy(os.path.join(webpack, "bundle.css"), os.path.join(cur_dir, "bundle.css"))
    shutil.copy(os.path.join(webpack, "fire.png"), os.path.join(cur_dir, "fire.png"))


def setup(app):
    app.connect('builder-inited', run_apidoc)
    app.connect('builder-inited', copy_report)