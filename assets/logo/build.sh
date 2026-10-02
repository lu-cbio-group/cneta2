#!/usr/bin/env bash
# Render every CNETA logo variant from cneta-logo.tex into svg/ and png/, plus
# favicon.ico.
#
# Needs: pdflatex (TeX Live "basic" is enough) and pdftocairo (poppler).
# favicon.ico additionally needs ImageMagick (`magick`); without it that one
# file is skipped.
#
#   ./build.sh
#
# svg/ and favicon.ico are committed, so nobody needs TeX just to *use* the
# logo (the docs build and GitHub render straight from them). Re-run this only
# after editing cneta-logo.tex, then commit the changed files. png/ is for
# one-off needs (avatars, slides), is gitignored, and is just regenerated here.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$here/svg" "$here/png"

# name                  layout   scheme  background  png width (px)
variants=(
  "cneta-logo               full     light   none   1200"
  "cneta-logo-dark          full     dark    none   1200"
  "cneta-logo-tagline       tagline  light   none   1200"
  "cneta-logo-tagline-dark  tagline  dark    none   1200"
  "cneta-icon               icon     light   none   512"
  "cneta-icon-dark          icon     dark    none   512"
  "cneta-icon-tile          icon     light   light  512"
  "cneta-icon-tile-dark     icon     dark    dark   512"
  "cneta-favicon            favicon  dark    dark   512"
)

render() {  # name layout scheme bg -> $tmp/name.pdf
  local name=$1 layout=$2 scheme=$3 bg=$4
  (cd "$here" && pdflatex -interaction=nonstopmode -halt-on-error \
      -output-directory="$tmp" -jobname="$name" \
      "\def\LogoLayout{$layout}\def\LogoScheme{$scheme}\def\LogoBg{$bg}\input{cneta-logo}" \
      >"$tmp/$name.stdout") || { tail -30 "$tmp/$name.stdout"; exit 1; }
}

for v in "${variants[@]}"; do
  read -r name layout scheme bg width <<<"$v"
  render "$name" "$layout" "$scheme" "$bg"
  # SVG: glyphs are converted to paths, so there is no font dependency.
  pdftocairo -svg "$tmp/$name.pdf" "$here/svg/$name.svg"
  # PNG: transparent wherever the variant has no background tile.
  pdftocairo -png -singlefile -transp -scale-to-x "$width" -scale-to-y -1 \
    "$tmp/$name.pdf" "$here/png/$name"
  echo "built $name"
done

# favicon.ico: rasterised straight from the vector favicon at each size
# (downscaling the 512 px PNG would blur the 16 px entry).
if command -v magick >/dev/null; then
  for s in 16 32 48; do
    pdftocairo -png -singlefile -transp -scale-to "$s" \
      "$tmp/cneta-favicon.pdf" "$tmp/favicon-$s"
  done
  magick "$tmp/favicon-16.png" "$tmp/favicon-32.png" "$tmp/favicon-48.png" \
    "$here/favicon.ico"
  echo "built favicon.ico"
else
  echo "magick not found: skipped favicon.ico" >&2
fi
