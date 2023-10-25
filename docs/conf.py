# SPDX-License-Identifier: LGPL-3.0-or-later
"""Configuration file for the Sphinx documentation builder."""
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('..'))
from datetime import datetime

# -- Project information -----------------------------------------------------

project = "ReacNetGenerator"
copyright = "2019-%d, East China Normal University" % datetime.now().year
author = "Jinzhe Zeng"

# The full version, including alpha/beta/rc tags
# release = 'reacnetgenerator'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinxarg.ext",
    "myst_nb",
    "numpydoc",
    "sphinx_favicon",
    "deepmodeling_sphinx",
    "sphinxcontrib.bibtex",
]

# Add any paths that contain templates here, relative to this directory.
# templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "pydata_sphinx_theme"
html_logo = "_static/reacnetgen.svg"
html_static_path = ["_static"]
html_js_files = ["https://unpkg.com/bilitube@0/dist/bilitube.min.js"]
html_css_files = ["css/custom.css"]
html_extra_path = ["report.html", "fire.png", "bundle.js", "bundle.css"]

html_theme_options = {
    "use_edit_page_button": True,
    "icon_links": [
        {
            "name": "Paper",
            "url": "https://doi.org/10.1039/C9CP05091D",
            "icon": "fa-solid fa-book-open",
            "type": "fontawesome",
        },
    ],
    "external_links": [{"name": "Report Loader", "url": "/report.html"}],
    "github_url": "https://github.com/tongzhugroup/reacnetgenerator",
    "logo": {
        "text": "ReacNetGenerator",
        "alt_text": "ReacNetGenerator",
    },
}

html_context = {
    "github_user": "tongzhugroup",
    "github_repo": "reacnetgenerator",
    "github_version": "master",
    "doc_path": "docs",
}

myst_heading_anchors = 3

favicons = [
    {
        "rel": "icon",
        "static-file": "reacnetgen.svg",
        "type": "image/svg+xml",
    },
]

enable_deepmodeling = False

myst_enable_extensions = [
    "dollarmath",
    "colon_fence",
]
mathjax_path = (
    "https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.min.js"
)
mathjax_options = {
    "integrity": "sha512-6FaAxxHuKuzaGHWnV00ftWqP3luSBRSopnNAA2RvQH1fOfnF/A1wOfiUWF7cLIOFcfb1dEhXwo5VG3DAisocRw==",
    "crossorigin": "anonymous",
}
mathjax3_config = {
    "loader": {"load": ["[tex]/mhchem"]},
    "tex": {"packages": {"[+]": ["mhchem"]}},
}

execution_mode = "off"
numpydoc_show_class_members = False

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

intersphinx_mapping = {
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "python": ("https://docs.python.org/", None),
    "ase": ("https://wiki.fysik.dtu.dk/ase/", None),
}

bibtex_bibfiles = ["refs.bib"]


def run_apidoc(_):
    """Run sphinx-apidoc on every build."""
    import os
    import sys

    from sphinx.ext.apidoc import main

    sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    module = os.path.join(cur_dir, "..", "reacnetgenerator")
    main(
        [
            "-M",
            "--tocfile",
            "api",
            "-H",
            "Python API",
            "-o",
            os.path.join(cur_dir, "api"),
            module,
            "--force",
        ]
    )


def run_node_command(args, **kwargs):
    """Call node with subprocess."""
    import subprocess

    return subprocess.call(["node", *args], **kwargs)


def copy_report(app):
    """Copy report.html to the build directory."""
    import os
    from pathlib import Path

    import yaml
    from sphinx.util.fileutil import copy_asset_file

    cur_dir = Path(__file__).parent.absolute()
    webpack = Path(cur_dir) / ".." / "reacnetgenerator" / "static" / "webpack"
    with open(webpack / ".yarnrc.yml") as f:
        yarn_path = Path(yaml.load(f, Loader=yaml.Loader)["yarnPath"])
    run_node_command([yarn_path], cwd=webpack)
    run_node_command(
        [yarn_path, "start"],
        cwd=webpack,
        env={**os.environ, "REACNETGENERATOR_BUILDWEB": "1"},
    )
    outdir = app.builder.outdir
    # first create the directory..
    os.makedirs(outdir, exist_ok=True)
    copy_asset_file(
        os.path.join(webpack, "bundle.html"), os.path.join(outdir, "report.html")
    )
    copy_asset_file(
        os.path.join(webpack, "bundle.js"), os.path.join(outdir, "bundle.js")
    )
    copy_asset_file(
        os.path.join(webpack, "bundle.css"), os.path.join(outdir, "bundle.css")
    )
    copy_asset_file(os.path.join(webpack, "fire.png"), os.path.join(outdir, "fire.png"))


def setup(app):
    """Set up sphinx.

    Parameters
    ----------
    app : sphinx.application.Sphinx
        The sphinx application.
    """
    app.connect("builder-inited", run_apidoc)
    app.connect("builder-inited", copy_report)
