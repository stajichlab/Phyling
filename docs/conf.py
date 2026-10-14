"""Sphinx configuration for the Phyling documentation."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from string import Template

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
exclude_patterns = ["_build", "_diagrams", "_generated", "Thumbs.db", ".DS_Store"]

myst_enable_extensions = ["colon_fence", "deflist", "tasklist", "attrs_inline"]
myst_heading_anchors = 3

graphviz_output_format = "svg"

# Render each diagram template once per Furo color scheme. The pages include both and Furo's only-light / only-dark
# classes show the one matching the active theme. The SVG background is transparent so it blends into either theme.
DIAGRAM_PALETTES = {
    "light": {
        "text": "#1f2933",
        "line": "#5b6b7f",
        "cluster": "#8a99ad",
        "node_fill": "#f5f7fa",
        "decision_fill": "#e8eef7",
        "final_fill": "#dff0e3",
        "final_line": "#4f8a5b",
    },
    "dark": {
        "text": "#e3e8ef",
        "line": "#9aa8b8",
        "cluster": "#6c7a8c",
        "node_fill": "#262c35",
        "decision_fill": "#25344a",
        "final_fill": "#1f3a29",
        "final_line": "#6fbf80",
    },
}

_docs_dir = Path(__file__).parent
(_docs_dir / "_generated").mkdir(exist_ok=True)
for _template in (_docs_dir / "_diagrams").glob("*.dot.in"):
    _name = _template.name.removesuffix(".dot.in")
    for _scheme, _palette in DIAGRAM_PALETTES.items():
        _content = Template(_template.read_text()).substitute(_palette)
        _out = _docs_dir / "_generated" / f"{_name}-{_scheme}.dot"
        # Only rewrite on change so Sphinx does not rebuild the diagrams every time
        if not _out.exists() or _out.read_text() != _content:
            _out.write_text(_content)

html_theme = "furo"
html_title = f"Phyling {release}"
html_static_path = ["_static"]
html_theme_options = {
    "source_repository": "https://github.com/stajichlab/Phyling/",
    "source_branch": "main",
    "source_directory": "docs/",
}
