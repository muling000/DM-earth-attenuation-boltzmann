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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import sphinx_rtd_theme

# -- Project information -----------------------------------------------------

# other project information is handled by cmake
copyright = '2021, Chen Xia'
author = 'Chen Xia'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.viewcode",
    "sphinx.ext.graphviz",
    "breathe",
    "myst_parser",
    "sphinx_rtd_theme"
]
autosummary_generate = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

math_number_all = True
numfig = True

source_suffix = {'.rst': 'restructuredtext', '.md': 'markdown'}

math_eqref_format = "Eq. ({number})"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


def setup(app):
    app.add_css_file('custom.css')


html_last_updated_fmt = "%Y-%m-%d"

# -- Options for LaTeX output -------------------------------------------------

latex_engine = 'pdflatex'
latex_documents = [('index', 'darkprop-manual.tex', 'DarkProp Manual',
                    author, 'manual', True)]
latex_elements = {
    'papersize': 'a4paper',
    'pointsize': '10pt',
    'extraclassoptions': 'openany,oneside',
    'preamble': r'\setcounter{tocdepth}{1}',
    'releasename': 'Version'
}
