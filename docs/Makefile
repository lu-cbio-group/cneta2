# Makefile for the CNETA documentation site
# -----------------------------------------
# Common targets:
#   make html        → build HTML into build/html
#   make live        → live-reload server at http://localhost:8080
#   make serve       → serve build/html without rebuilding
#   make api         → regenerate the C++ API XML via Doxygen
#   make linkcheck   → check external links
#   make strict      → build with -n -W (nitpicky, warnings-as-errors)
#   make all         → clean + api + html
#   make clean       → remove build/
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

# Build HTML (delegates to sphinx "make mode")
html:
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Live-reload while editing (requires sphinx-autobuild)
live:
	@python -m sphinx_autobuild "$(SOURCEDIR)" "$(BUILDDIR)/html" \
		$(LIVE_FLAGS) $(SPHINXOPTS)

# Serve the already-built site (no rebuilds)
serve:
	@python -m http.server $(PORT) -d "$(BUILDDIR)/html"

# Strict build: nitpicky + warnings-as-errors
strict:
	@$(MAKE) html SPHINXOPTS="-n -W $(SPHINXOPTS)"

# Regenerate the C++ API XML. conf.py picks up build/doxygen/xml automatically.
# FORCE_DOXYGEN makes conf.py refresh even when XML already exists.
# The mkdir is required: on a fresh checkout Doxygen's OUTPUT_DIRECTORY only
# creates one missing path segment, not a nested one, so it fails to create
# "build/doxygen" when "build/" itself doesn't exist yet.
api:
	@mkdir -p "$(BUILDDIR)/doxygen"
	@FORCE_DOXYGEN=1 $(DOXYGEN) Doxyfile

# Full rebuild including the C++ API
all: clean api html

# Link checker (run from time to time - on campus or using VPN)
linkcheck:
	@$(SPHINXBUILD) -M linkcheck "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Clean build artifacts
clean:
	@rm -rf "$(BUILDDIR)"

.PHONY: help html live serve strict api all linkcheck clean Makefile

# Catch-all: route unknown targets to Sphinx "make mode"
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
