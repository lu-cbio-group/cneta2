#!/usr/bin/env bash
# Render the two figures, cneta-overview.tex (one glance) and
# cneta-pipeline.tex (every file and mode), into pdf/ (to look at while
# editing), svg/ (light and dark, committed) and png/ (slides, chat).
#
# Needs: pdflatex (TeX Live "basic" or MacTeX is enough). With pdftocairo
# (poppler) installed it also writes svg/ and png/; without it you still
# get pdf/, which is all you need while editing.
#
#   ./build.sh
#
# svg/ is committed, so nobody needs TeX to *use* the figures. Re-run this
# after editing a .tex file, then commit the changed SVGs together with the
# source. pdf/ and png/ are gitignored and are just regenerated here.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
command -v pdflatex >/dev/null || {
  echo "pdflatex not found: install TeX Live (basic) or MacTeX" >&2; exit 1; }
have_cairo=1
command -v pdftocairo >/dev/null || have_cairo=0

# Build files go in a folder next to the source (gitignored), not in a
# system temporary directory, so this also works where those are locked down.
tmp="$here/.build"
rm -rf "$tmp"; mkdir -p "$tmp"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$here/pdf"
if [ "$have_cairo" -eq 1 ]; then mkdir -p "$here/svg" "$here/png"; fi

for source in cneta-overview cneta-pipeline; do
  for scheme in light dark; do
    name=$source
    [ "$scheme" = dark ] && name=$source-dark
    (cd "$here" && pdflatex -interaction=nonstopmode -halt-on-error \
        -output-directory="$tmp" -jobname="$name" \
        "\def\PipelineScheme{$scheme}\input{$source}" \
        >"$tmp/$name.stdout") || { tail -30 "$tmp/$name.stdout"; exit 1; }
    cp "$tmp/$name.pdf" "$here/pdf/$name.pdf"
    if [ "$have_cairo" -eq 1 ]; then
      # SVG: glyphs become paths, so there is no font dependency.
      pdftocairo -svg "$tmp/$name.pdf" "$here/svg/$name.svg"
      # PNG: transparent background, 2000 px wide.
      pdftocairo -png -singlefile -transp -scale-to-x 2000 -scale-to-y -1 \
        "$tmp/$name.pdf" "$here/png/$name"
    fi
    echo "built $name"
  done
done

if [ "$have_cairo" -eq 0 ]; then
  echo "pdftocairo not found: wrote pdf/ only; install poppler to refresh svg/ and png/" >&2
fi
