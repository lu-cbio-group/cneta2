#!/usr/bin/env bash
# Render the two figures, cneta-overview.tex (one glance) and
# cneta-pipeline.tex (every file and mode), into svg/ (light and dark) and png/.
#
# Needs: pdflatex (TeX Live "basic" is enough) and pdftocairo (poppler).
#
#   ./build.sh
#
# svg/ is committed, so nobody needs TeX to *use* the figure. Re-run this
# after editing cneta-pipeline.tex, then commit the changed SVGs. png/ is for
# slides and chat, is gitignored, and is just regenerated here.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
# Build files go in a folder next to the source (gitignored), not in a
# system temporary directory, so this also works where those are locked down.
tmp="$here/.build"
rm -rf "$tmp"; mkdir -p "$tmp"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$here/svg" "$here/png"

for source in cneta-overview cneta-pipeline; do
  for scheme in light dark; do
    name=$source
    [ "$scheme" = dark ] && name=$source-dark
    (cd "$here" && pdflatex -interaction=nonstopmode -halt-on-error \
        -output-directory="$tmp" -jobname="$name" \
        "\def\PipelineScheme{$scheme}\input{$source}" \
        >"$tmp/$name.stdout") || { tail -30 "$tmp/$name.stdout"; exit 1; }
    # SVG: glyphs become paths, so there is no font dependency.
    pdftocairo -svg "$tmp/$name.pdf" "$here/svg/$name.svg"
    # PNG: transparent background, 2000 px wide.
    pdftocairo -png -singlefile -transp -scale-to-x 2000 -scale-to-y -1 \
      "$tmp/$name.pdf" "$here/png/$name"
    echo "built $name"
  done
done
