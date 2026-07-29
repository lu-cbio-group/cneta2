"""Sphinx configuration for the CNETA documentation site.

Build locally with:

    conda activate cneta-docs
    make -C docs live      # live-reload preview
    make -C docs strict    # what CI runs

See docs/developer-guide/documentation.md for the full authoring guide.
"""

import os
import subprocess
from pathlib import Path

# ---------------------------------------------------------------------------
# Project information
# ---------------------------------------------------------------------------

project = "CNETA"
author = "CNETA team"
copyright = "2026, CNETA contributors"

# Bump on release; also drives the version shown in the sidebar.
# Kept in sync with the `project(cneta VERSION ...)` call in code/CMakeLists.txt.
release = "2.0"
version = release

SOURCE_DIR = Path(__file__).parent
DOCS_DIR = SOURCE_DIR.parent
DOXYGEN_XML = DOCS_DIR / "build" / "doxygen" / "xml"

# ---------------------------------------------------------------------------
# Extensions
# ---------------------------------------------------------------------------

extensions = [
    "myst_parser",              # Markdown authoring
    "sphinxcontrib.mermaid",    # workflow diagrams
    "breathe",                  # C++ API via Doxygen XML
    "sphinx_copybutton",        # copy button on code blocks
    "sphinx_design",            # tab-sets used in developer-guide/documentation.md
    "sphinx-tabs",              # tabbed content (e.g. macOS / Linux / HPC install)
    "sphinx-togglebutton",      # collapsible admonitions
    "sphinx-notfound-page",     # correct 404 page for a project served under /REPO/
    "sphinx-sitemap",           # sitemap.xml; requires html_baseurl to be set
]

exclude_patterns = ["build", "Thumbs.db", ".DS_Store"]

# `make strict` runs with -n (nitpicky), which turns every unresolved
# cross-reference into a warning-as-error. Breathe generates a cpp:identifier
# cross-reference for every type in every function signature it renders, and
# with EXTRACT_ALL=YES and only two headers exposed via {doxygenfile} so far
# (see api/cpp.md), most of those types are never independently documented —
# STL containers, types from headers not yet exposed, forward-declared
# classes. Failing the build on those would block every doc build until the
# whole `code/` surface is documented, which defeats the point of shipping
# "structural" API coverage early. Narrowed to cpp:identifier only, so a
# genuine mistake in a MyST cross-reference (e.g. a broken relative link)
# still fails strict as intended.
nitpick_ignore_regex = [
    ("cpp:identifier", r".*"),
    # Breathe also emits internal :ref: links to the compound pages of
    # types/files it mentions but that aren't exposed via their own
    # {doxygenfile}/{doxygenstruct} directive yet (e.g. a struct defined in
    # a header not shown in api/cpp.md). Doxygen's mangled anchor names are
    # distinctive enough to target narrowly: `..._8hpp`/`..._8cpp` for files,
    # `struct...`/`class...` for types (CASE_SENSE_NAMES=YES in Doxyfile
    # means no underscore directly after the struct/class keyword, e.g.
    # "structCN__CHANGE" rather than "struct_c_n___c_h_a_n_g_e").
    ("std:ref", r".*_8[a-z]+$"),
    ("std:ref", r"^(struct|class).*"),
]

# ---------------------------------------------------------------------------
# Markdown (MyST)
# ---------------------------------------------------------------------------

myst_enable_extensions = [
    "colon_fence",   # ::: fences, friendlier than ``` for directives
    "deflist",       # definition lists, handy for option/parameter tables
    "linkify",       # bare URLs become links
    "substitution",  # reusable values via myst_substitutions
    "dollarmath",    # $...$ inline maths
]

# Auto-generate anchors for ##, ###, so cross-page links to headings work.
myst_heading_anchors = 3

# The important one: a plain ```mermaid fenced block is treated as the
# {mermaid} directive. This means diagrams render both on GitHub (in the
# raw Markdown) and here, with no rewriting.
myst_fence_as_directive = ["mermaid"]

# Values usable in pages as {{ repo_url }} etc.
myst_substitutions = {
    "repo_url": "https://github.com/ORG/REPO",  # TODO: set to the real repo
}

# ---------------------------------------------------------------------------
# HTML output
# ---------------------------------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "navigation_depth": 4,
    "collapse_navigation": False,
    "titles_only": False,
    "style_external_links": True,
}
html_static_path = ["_static"]
html_title = f"{project} {release}"

# TODO: set once GitHub Pages is live, e.g.
#   https://ORG.github.io/REPO/
# Used for canonical links and sitemap generation.
html_baseurl = os.environ.get("DOCS_BASEURL", "")

# ---------------------------------------------------------------------------
# C++ API (Doxygen + Breathe)
# ---------------------------------------------------------------------------

breathe_projects = {"cneta": str(DOXYGEN_XML)}
breathe_default_project = "cneta"
breathe_default_members = ("members",)


def _maybe_run_doxygen() -> bool:
    """Generate Doxygen XML if it is missing, or if FORCE_DOXYGEN is set.

    Kept conditional so that ``make live`` does not re-run Doxygen over the
    whole of ``code/`` on every keystroke. Use ``make api`` to refresh it
    deliberately; CI always runs it because ``build/`` starts empty.

    Returns whether the XML is present and usable afterwards.
    """
    if DOXYGEN_XML.exists() and not os.environ.get("FORCE_DOXYGEN"):
        return True
    try:
        # Doxygen creates at most one missing path segment for
        # OUTPUT_DIRECTORY, so on a fresh checkout (no build/ yet) it fails
        # to create the nested "build/doxygen" and silently produces no XML.
        DOXYGEN_XML.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["doxygen", "Doxyfile"], cwd=DOCS_DIR, check=True)
        return True
    except FileNotFoundError:
        print(
            "[conf.py] Doxygen skipped: the 'doxygen' binary is not on PATH. "
            "API pages will be empty. `pip install -r requirements.txt` does "
            "NOT install it (it isn't a Python package) — either "
            "`conda activate cneta-docs` (installs it via environment.yml), "
            "or install it yourself (`brew install doxygen` / "
            "`apt install doxygen`) if you're on the venv/uv workflow."
        )
        return False
    except subprocess.CalledProcessError as exc:
        print(f"[conf.py] Doxygen failed ({exc}); API pages will be empty.")
        return False


# If Doxygen didn't run (not installed — e.g. the wrong conda env is active
# — or it failed), there is no XML for Breathe to read. Breathe does not
# degrade gracefully in that case: it raises breathe.file_state_cache.
# MTimeError deep inside a parallel Sphinx worker, which surfaces as an
# opaque SphinxParallelError and crashes the whole build. Dropping "breathe"
# from the extension list instead means {doxygenfile} in api/cpp.md is just
# an unrecognized directive — a normal (non-fatal outside `make strict`)
# warning, and the rest of the site still builds.
if not _maybe_run_doxygen():
    extensions.remove("breathe")
