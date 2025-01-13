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
import os
import sys

sys.path.insert(0, os.path.abspath("../"))

# -- Project information -----------------------------------------------------

project = "x3d2"
copyright = "2025, Xcompact3d Team"
author = "Xcompact3d Team"
master_doc = "index"

# -- Get version information and date from Git ----------------------------

try:
    from subprocess import check_output

    release = check_output(["git", "describe", "--tags", "--always"])
    release = release.decode().strip()
    today = check_output(["git", "show", "-s", "--format=%ad", "--date=short"])
    today = today.decode().strip()
except Exception:
    release = "<unknown>"
    today = "<unknown date>"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx_rtd_theme',
    'sphinx.ext.graphviz'
    #,
#    'sphinx_click.ext'
    #,
#    'recommonmark',
#    'nbsphinx',
#    'sphinx_rtd_theme'
]

source_suffix = [".rst", ".md"]

#intersphinx_mapping = {
#    "toolbox": ("https://xcompact3d.readthedocs.io/en/stable/", None,),
#}

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]