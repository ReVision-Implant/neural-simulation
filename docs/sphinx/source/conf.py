import os
import sys
sys.path.insert(0, os.path.join(os.path.abspath(__file__),r'../../../../v1_Anke/components/'))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'neural-simulation'
copyright = '2023, Nils Van Rompaey'
author = 'Nils Van Rompaey'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    # 'autoapi.extension',
]

autoapi_type = 'python'
autoapi_dirs = ['../../../v1_Anke/components']
templates_path = ['_templates']
exclude_patterns = []
autosummary_generate = True


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

import sphinx_rtd_theme
html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_static_path = ['_static']
html_css_files = ['my_theme.css']
autoclass_content = 'both'
autoapi_root = 'code'

# Automatic generation of autosummary
from sphinx.ext.autosummary import Autosummary
from sphinx.ext.autosummary import get_documenter
from docutils.parsers.rst import directives
from sphinx.util.inspect import safe_getattr
import re