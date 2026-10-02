# CNETA logo

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="svg/cneta-logo-dark.svg">
    <img src="svg/cneta-logo.svg" alt="CNETA" width="320">
  </picture>
</p>

The mark is a rooted phylogeny whose tips run into copy-number profiles drawn
as heatmap rows, the way copy-number trees are usually plotted. Blocks are
coloured by copy-number state:

| State | Colour |
|---|---|
| loss | blue `#3B82C4` |
| neutral | grey `#CBD5DE` (on dark: `#5B6E7E`) |
| gain | orange `#F08A3C` |
| amplification | red `#D1382E` |

Tree and wordmark are slate `#1D2B36` on light backgrounds and white on dark
ones. The wordmark is Latin Modern Sans Bold, converted to outlines in the
exports, so no font is needed to use them.

## Which file to use

The variant name says what the logo is **for**, not its ink colour:
`cneta-logo` goes on light backgrounds, `cneta-logo-dark` on dark ones. All
variants except the tiles and the favicon have a transparent background.
Prefer `svg/`. `png/` (built locally by `./build.sh`, not committed) is for
places that reject SVG.

| Use | File |
|---|---|
| Light background | `svg/cneta-logo.svg` |
| Dark background | `svg/cneta-logo-dark.svg` |
| Same, with the "copy number evolutionary tree analysis" tagline | `svg/cneta-logo-tagline.svg`, `svg/cneta-logo-tagline-dark.svg` |
| Mark only (no wordmark) | `svg/cneta-icon.svg`, `svg/cneta-icon-dark.svg` |
| Background you do not control (avatars, social previews) | `svg/cneta-icon-tile.svg`, `svg/cneta-icon-tile-dark.svg` (or the `png/` equivalents): mark on a rounded tile |
| Browser tab icon | `favicon.ico`, `svg/cneta-favicon.svg` |

The favicon is a simplified mark (three tips, chunkier strokes) rather than a
shrunken logo, because the full mark turns to mush at 16 px.

### Where it is used

- **`README.md`** (repository root): a `<picture>` element that shows
  `cneta-logo.svg` or `cneta-logo-dark.svg` depending on the reader's GitHub
  theme. This is the only place that needs both.
- **Docs sidebar** (`docs/source/conf.py`, `html_logo`): `cneta-logo-dark.svg`.
  The theme's sidebar header is always dark, whatever the reader's setting, so
  there is nothing to switch. `conf.py` also sets the header colour to the
  logo's slate; on the theme's default blue the logo's blue "loss" blocks
  disappear.
- **Docs favicon** (`conf.py`, `html_favicon`): `favicon.ico`.

`conf.py` points straight at the files here, so there is a single copy of each
export. The docs workflow also runs when `assets/logo/**` changes.

## Rebuilding

Only needed after editing `cneta-logo.tex`. `svg/` and `favicon.ico` are
committed, so building the docs or reading the README never requires TeX.

```bash
./build.sh
```

Requires `pdflatex` (TeX Live "basic" is enough: TikZ and Latin Modern, no
extra packages) and `pdftocairo` (poppler). ImageMagick (`magick`) is
optional and only used for `favicon.ico`; without it that file is skipped.

Commit the regenerated `svg/` files and `favicon.ico` together with the source
change. `png/` is gitignored.

`cneta-logo.tex` is the only source. Its header documents the options
(`\LogoLayout`, `\LogoScheme`, `\LogoBg`) if you want to render a single
variant by hand, and the mark's geometry is plain TikZ coordinates in a
10 x 10 mm box.
