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
from urllib.parse import urlparse

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

# Site language: "en" (default) or "zh_CN". Each language is a separate build
# into build/html/<language>/ (see `make html-all`). Read from DOCS_LANG
# rather than `-D language=...` because html_baseurl below depends on it, and
# conf.py cannot see -D overrides. Translations live in locale/<language>/
# as .po files; untranslated strings fall back to the English source.
language = os.environ.get("DOCS_LANG", "en")
locale_dirs = ["locale/"]
gettext_compact = False     # one .po file per page, not one per directory

SOURCE_DIR = Path(__file__).parent
DOCS_DIR = SOURCE_DIR.parent
DOXYGEN_XML = DOCS_DIR / "build" / "doxygen" / "xml"

# ---------------------------------------------------------------------------
# Extensions
# ---------------------------------------------------------------------------

extensions = [
    "myst_parser",              # Markdown authoring
    "sphinxcontrib.mermaid",    # workflow diagrams
    "sphinxcontrib.images",     # click-to-enlarge lightbox for figures
    "breathe",                  # C++ API via Doxygen XML
    "sphinx_copybutton",        # copy button on code blocks
    "sphinx_design",            # tab-sets used in developer-guide/documentation.md
    "sphinx_togglebutton",      # collapsible admonitions
    "notfound.extension",     # correct 404 page for a project served under /REPO/
    "sphinx_sitemap",           # sitemap.xml; requires html_baseurl to be set
]

exclude_patterns = ["build", "locale", "Thumbs.db", ".DS_Store"]

# _templates/layout.html overrides the theme layout to add a "View on
# GitHub" link pinned to the bottom of the sidebar, for visitors who land
# on the docs site directly rather than coming from the repo.
templates_path = ["_templates"]

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
    "repo_url": "https://github.com/lu-cbio-group/cneta2",
}


# ---------------------------------------------------------------------------
# HTML output
# ---------------------------------------------------------------------------

html_theme = "furo"

# Logo and favicon live in assets/logo/ at the repo root, shared with
# README.md; the exports there are committed, so no TeX is needed to build the
# docs. Furo has a light and a dark mode, so it takes both logo variants by
# file name from the static path; see assets/logo/README.md.
LOGO_DIR = DOCS_DIR.parent / "assets" / "logo"
html_static_path = ["_static", str(LOGO_DIR / "svg")]
html_favicon = str(LOGO_DIR / "favicon.ico")

html_theme_options = {
    "light_logo": "cneta-logo.svg",
    "dark_logo": "cneta-logo-dark.svg",
    # The logo already contains the wordmark, so don't repeat the project name.
    "sidebar_hide_name": True,
    # No "view source" / "edit this page" icons; the top-right row holds the
    # language switch, GitHub link and light/dark toggle (_templates/page.html).
    "top_of_page_buttons": [],
    "footer_icons": [
        {
            "name": "GitHub",
            "url": myst_substitutions["repo_url"],
            "html": (
                '<svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">'
                '<path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.01 8.01 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>'
                "</svg>"
            ),
            "class": "",
        },
    ],
}

# Exposed to templates (_templates/header-links.html).
html_context = {
    "repo_url": myst_substitutions["repo_url"],
}
html_css_files = ["custom.css"]
html_title = f"{project} {release}"

# Used for canonical links and sitemap generation. Each language build is
# served from its own subdirectory, so the base URL includes it.
SITE_URL = os.environ.get("DOCS_BASEURL", "https://lu-cbio-group.github.io/cneta2/")
html_baseurl = f"{SITE_URL.rstrip('/')}/{language}/"
# sphinx-sitemap's default scheme prepends {lang}{version}, which would
# duplicate the language already in html_baseurl.
sitemap_url_scheme = "{link}"

# sphinx-not-found-page settings
# for more, see https://sphinx-notfound-page.readthedocs.io/en/latest/configuration.html
# The 404 page is served for any missing URL, so its links must be absolute
# from the site root: e.g. /cneta2/en/_static/... on GitHub Pages.
notfound_urls_prefix = urlparse(html_baseurl).path


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
