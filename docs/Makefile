# Makefile for the CNETA documentation site
# -----------------------------------------
# Common targets:
#   make html        → build the English site into build/html/en
#   make html-zh     → build the Chinese site into build/html/zh_CN
#   make html-all    → build every language
#   make live        → live-reload server (English) at http://localhost:8080
#   make live-zh     → live-reload server (Chinese) at http://localhost:8080
#   make serve       → serve build/html without rebuilding (language switch works here)
#   make api         → regenerate the C++ API XML via Doxygen
#   make linkcheck   → check external links
#   make strict      → build every language with -n -W (nitpicky, warnings-as-errors)
#   make all         → clean + api + html-all
#   make clean       → remove build/
#
# Translation (see developer-guide/documentation.md):
#   make intl-update → extract English strings and update source/locale/*/LC_MESSAGES/*.po
#
# Local setup:
#   conda env create -f environment.yml
#   conda activate cneta-docs
#   make live

SHELL := /bin/bash

# You can set these from the command line or env
SPHINXBUILD ?= sphinx-build
SPHINXOPTS  ?= -j auto
SOURCEDIR   ?= source
BUILDDIR    ?= build
DOXYGEN     ?= doxygen

# Languages other than English; each is built into build/html/<lang>/.
TRANSLATIONS ?= zh_CN

# Live-reload settings (used by `make live`)
PORT       ?= 8080

# Build up optional watch flags only if the paths exist
WATCH_DIRS :=
ifneq ("$(wildcard $(SOURCEDIR)/_static)","")
WATCH_DIRS += --watch $(SOURCEDIR)/_static
endif
ifneq ("$(wildcard $(SOURCEDIR)/_templates)","")
WATCH_DIRS += --watch $(SOURCEDIR)/_templates
endif

LIVE_FLAGS ?= --open-browser --port $(PORT) \
              --ignore "$(BUILDDIR)/*" \
              $(WATCH_DIRS)

# Default target
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Build HTML. Each language is a separate Sphinx build (DOCS_LANG is read by
# conf.py) with its own doctree cache; site-root adds the files that live
# above the language directories.
html:
	@DOCS_LANG=en $(SPHINXBUILD) -b html -d "$(BUILDDIR)/doctrees/en" \
		"$(SOURCEDIR)" "$(BUILDDIR)/html/en" $(SPHINXOPTS) $(O)
	@$(MAKE) --no-print-directory site-root

html-zh:
	@DOCS_LANG=zh_CN $(SPHINXBUILD) -b html -d "$(BUILDDIR)/doctrees/zh_CN" \
		"$(SOURCEDIR)" "$(BUILDDIR)/html/zh_CN" $(SPHINXOPTS) $(O)

html-all: html html-zh

# The site root redirects to the English pages. GitHub Pages serves
# /404.html for any missing URL, so the English 404 page is copied up.
site-root:
	@printf '%s\n' '<!DOCTYPE html>' \
		'<meta charset="utf-8">' \
		'<title>CNETA documentation</title>' \
		'<meta http-equiv="refresh" content="0; url=en/index.html">' \
		'<link rel="canonical" href="en/index.html">' \
		'<a href="en/index.html">CNETA documentation</a>' \
		> "$(BUILDDIR)/html/index.html"
	@cp "$(BUILDDIR)/html/en/404.html" "$(BUILDDIR)/html/404.html"
	@find "$(BUILDDIR)/html/en" -type f -name '*.html' \
		! -path "$(BUILDDIR)/html/en/index.html" \
		! -path "$(BUILDDIR)/html/en/404.html" -print0 | \
		while IFS= read -r -d '' page; do \
			relative=$${page#$(BUILDDIR)/html/en/}; \
			target="en/$$relative"; \
			directory="$$relative"; \
			while [[ "$$directory" == */* ]]; do \
				target="../$$target"; \
				directory=$${directory%/*}; \
			done; \
			mkdir -p "$(BUILDDIR)/html/$$(dirname "$$relative")"; \
			printf '%s\n' '<!DOCTYPE html>' '<meta charset="utf-8">' \
				'<title>Page moved</title>' \
				"<script>window.location.replace(\"$$target\" + window.location.search + window.location.hash);</script>" \
				"<a href=\"$$target\">This page has moved.</a>" \
				> "$(BUILDDIR)/html/$$relative"; \
		done

# Live-reload while editing (requires sphinx-autobuild). It serves one
# language only, so the language switch links lead nowhere here; use
# `make html-all serve` to check them.
live:
	@DOCS_LANG=en python -m sphinx_autobuild "$(SOURCEDIR)" "$(BUILDDIR)/html/en" \
		-d "$(BUILDDIR)/doctrees/en" $(LIVE_FLAGS) $(SPHINXOPTS)

live-zh:
	@DOCS_LANG=zh_CN python -m sphinx_autobuild "$(SOURCEDIR)" "$(BUILDDIR)/html/zh_CN" \
		-d "$(BUILDDIR)/doctrees/zh_CN" $(LIVE_FLAGS) $(SPHINXOPTS)

# Extract translatable strings, then add new/changed ones to each .po file.
# Changed strings are marked "fuzzy" and shown in English until re-translated.
gettext:
	@$(SPHINXBUILD) -b gettext "$(SOURCEDIR)" "$(BUILDDIR)/gettext" $(SPHINXOPTS) $(O)

intl-update: gettext
	@sphinx-intl update -p "$(BUILDDIR)/gettext" -d "$(SOURCEDIR)/locale" \
		$(addprefix -l ,$(TRANSLATIONS))

# Serve the already-built site (no rebuilds)
serve:
	@python -m http.server $(PORT) -d "$(BUILDDIR)/html"

# Strict build: nitpicky + warnings-as-errors
strict:
	@$(MAKE) html-all SPHINXOPTS="-n -W $(SPHINXOPTS)"

# Regenerate the C++ API XML. conf.py picks up build/doxygen/xml automatically.
# FORCE_DOXYGEN makes conf.py refresh even when XML already exists.
# The mkdir is required: on a fresh checkout Doxygen's OUTPUT_DIRECTORY only
# creates one missing path segment, not a nested one, so it fails to create
# "build/doxygen" when "build/" itself doesn't exist yet.
api:
	@mkdir -p "$(BUILDDIR)/doxygen"
	@FORCE_DOXYGEN=1 $(DOXYGEN) Doxyfile

# Full rebuild including the C++ API
all: clean api html-all

# Link checker (run from time to time - on campus or using VPN)
linkcheck:
	@$(SPHINXBUILD) -M linkcheck "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Clean build artifacts
clean:
	@rm -rf "$(BUILDDIR)"

.PHONY: help html html-zh html-all site-root live live-zh gettext intl-update serve strict api all linkcheck clean Makefile

# Catch-all: route unknown targets to Sphinx "make mode"
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
