"""Sphinx configuration for stangene documentation."""

import os
import sys

# Add source to path for autodoc
sys.path.insert(0, os.path.abspath("../src"))

project = "stangene"
copyright = "2026, Sijie Chen"
author = "Sijie Chen"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
    "myst_parser",
]

# MyST (Markdown) settings
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "tasklist",
]
myst_heading_anchors = 3

# Autodoc settings
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
autodoc_member_order = "bysource"
autosummary_generate = True

# Napoleon settings (Google/NumPy docstrings)
napoleon_google_docstring = True
napoleon_numpy_docstring = True

# Intersphinx
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "anndata": ("https://anndata.readthedocs.io/en/stable/", None),
}

# Theme
html_theme = "furo"
html_title = "stangene"
html_logo = "_static/logo.svg"
html_favicon = "_static/logo.svg"
html_theme_options = {
    "source_repository": "https://github.com/chansigit/stangene",
    "source_branch": "master",
    "source_directory": "docs/",
}

# Source file suffixes
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Suppress warnings for missing type stubs
suppress_warnings = ["autodoc.import_error"]

# Templates and static
templates_path = ["_templates"]
exclude_patterns = ["_build", "superpowers"]
html_static_path = ["_static"]
