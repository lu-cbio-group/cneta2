# CNETA logo

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="svg/cneta-logo-dark.svg">
    <img src="svg/cneta-logo.svg" alt="CNETA" width="320">
  </picture>
</p>

The mark combines a rooted phylogeny with copy-number profiles drawn as
heatmap rows, following the usual representation of copy-number trees.

| State | Colour |
|---|---|
| loss | blue `#3B82C4` |
| neutral | grey `#CBD5DE` (dark background: `#5B6E7E`) |
| gain | orange `#F08A3C` |
| amplification | red `#D1382E` |

The tree and wordmark are slate `#1D2B36` on light backgrounds and white on
dark backgrounds. The wordmark uses Latin Modern Sans Bold, converted to
outlines in the exported files, so no font installation is required.

## Choosing a file

Variant names describe the **background they are intended for**:

- `cneta-logo*` — light backgrounds
- `cneta-logo-dark*` — dark backgrounds

Use files from `svg/` whenever possible. PNG versions can be generated locally
with `./build.sh` for software that does not accept SVG.

All variants have transparent backgrounds except the tile versions and
`favicon.ico`.

| Use | File |
|---|---|
| Light background | `svg/cneta-logo.svg` |
| Dark background | `svg/cneta-logo-dark.svg` |
| Logo with tagline | `svg/cneta-logo-tagline.svg`, `svg/cneta-logo-tagline-dark.svg` |
| Mark only | `svg/cneta-icon.svg`, `svg/cneta-icon-dark.svg` |
| Avatar or uncontrolled background | `svg/cneta-icon-tile.svg`, `svg/cneta-icon-tile-dark.svg` |
| Browser favicon | `favicon.ico`, `svg/cneta-favicon.svg` |

The tile variants place the mark on a rounded background and are also
available as generated PNG files.

The favicon uses a simplified three-tip mark with thicker strokes so that it
remains legible at small sizes.

## Repository usage

- **`README.md`** — uses a `<picture>` element to switch between
  `cneta-logo.svg` and `cneta-logo-dark.svg` with the GitHub theme.
- **Documentation sidebar** — `docs/source/conf.py` uses
  `cneta-logo-dark.svg`. The sidebar is always dark, so no theme switching is
  needed.
- **Documentation favicon** — `conf.py` uses `favicon.ico`.

The documentation header colour is set to the logo's slate colour because the
theme's default blue gives insufficient contrast with the blue **loss** blocks.

`conf.py` references the files here directly, so each export has a single
source location. The documentation workflow also runs when files under
`assets/logo/**` change.

## Rebuilding

Rebuild the exports only after changing `cneta-logo.tex`:

```bash
./build.sh
```

The committed `svg/` files and `favicon.ico` mean that building the
documentation does **not** require TeX.

### Requirements

- `pdflatex` — a basic TeX Live installation with TikZ and Latin Modern
- `pdftocairo` — provided by Poppler
- `magick` — optional; used only to generate `favicon.ico`

If ImageMagick is unavailable, favicon generation is skipped.

After changing the source, commit the regenerated `svg/` files and
`favicon.ico` together with `cneta-logo.tex`. Generated files under `png/` are
gitignored.

`cneta-logo.tex` is the single source for all variants. Its header documents
the rendering options (`\LogoLayout`, `\LogoScheme`, and `\LogoBg`), and the
mark geometry is defined directly as TikZ coordinates in a 10 × 10 mm box.