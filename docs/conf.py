"""Sphinx configuration for the Phyling documentation."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

project = "Phyling"
author = "Jason Stajich, Cheng-Hung Tsai"
copyright = "2026, Jason Stajich and Cheng-Hung Tsai"

try:
    release = version("phyling")
except PackageNotFoundError:
    release = "unknown"
version = ".".join(release.split(".")[:2])

extensions = [
    "myst_parser",
    "sphinx_copybutton",
    "sphinx.ext.graphviz",
]

source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
root_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

myst_enable_extensions = ["colon_fence", "deflist", "tasklist", "attrs_inline"]
myst_heading_anchors = 3

graphviz_output_format = "svg"

html_theme = "furo"
html_title = f"Phyling {release}"
html_static_path = ["_static"]
html_theme_options = {
    "source_repository": "https://github.com/stajichlab/Phyling/",
    "source_branch": "main",
    "source_directory": "docs/",
}
