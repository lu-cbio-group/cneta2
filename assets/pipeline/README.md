# CNETA pipeline figures

Two figures, drawn in TikZ the same way as the logo in `../logo/`, so they
live under version control and come in a light and a dark variant. **First
drafts by the RSE; the research team is finishing them on this branch.**
The documentation does not use them yet.

## At a glance

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="svg/cneta-overview-dark.svg">
    <img src="svg/cneta-overview.svg" alt="CNETA at a glance" width="900">
  </picture>
</p>

What the tools do, on one line: input, the one evolutionary model, maximum
likelihood inference, the result; simulation from the same model; the
Bayesian sampler as a draft; and what comes next. Source:
`cneta-overview.tex`.

## The pipeline in detail

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="svg/cneta-pipeline-dark.svg">
    <img src="svg/cneta-pipeline.svg" alt="CNETA pipeline" width="900">
  </picture>
</p>

How the programs and the `util/` scripts chain together, and which files
cross the boundaries between them. Source: `cneta-pipeline.tex`.

| Use | File |
|---|---|
| Light background | `svg/cneta-overview.svg`, `svg/cneta-pipeline.svg` |
| Dark background | `svg/cneta-overview-dark.svg`, `svg/cneta-pipeline-dark.svg` |
| Slides, chat | `png/` (generated locally by `./build.sh`, not committed) |
| A look while editing | `pdf/` (generated locally by `./build.sh`, not committed) |

## Reading it

- **Blue** boxes are the programs (`cnets`, `cnetml` in its modes), **purple**
  the scripts in `util/`, **grey** the files that pass between them, **green**
  the results people look at, **amber** the settings to choose, **red
  dashed** what is not released yet (`cnetmcmc`).
- The grey files are the contract between steps: a Nextflow process or a
  Galaxy tool consumes and produces exactly these.
- Every program writes its version and seed into its outputs, so a result can
  be traced to the code and the random draw that produced it.

## What still needs the team

The drawings follow the code and the scripts as they are in the repository;
the words are the RSE's. Please check, and change directly in the `.tex`
files:

- the wording of every box: is each step, mode and script described the way
  you would describe it to a new student?
- the grey file boxes: are the file names right, is anything missing that a
  pipeline would need, is anything listed that nobody uses?
- which `util/` scripts belong in the figure, and under which name;
- the "Next" box in the overview: does it say what the group has decided?
- whether `cnetmcmc` should stay marked as a draft.

## Editing the figures

Everything is in the two `.tex` files; no drawing program is involved.

- **A box** is one line: `\node[style, text width=36mm] (name) at (x,y) {text};`.
  The style word sets the colour (`prog`, `file`, `script`, `result`,
  `settings`, `draft`); `(name)` is how arrows refer to the box; `(x,y)` is
  its centre in millimetres, x to the right and y upwards; the braces hold
  the text. Inside the text, `\\` starts a new line, `\textbf{...}` is bold
  and `\f{...}` sets a file name or option in monospace.
- **An arrow** is `\draw[arr] (from) -- (to);` using the box names. `arr` is
  an input to a step, `arrf` a file being written, `arro` something that
  feeds a result, `arrd` an alternative input.
- **To reword** a box, edit its braces. **To move** it, change `(x,y)`.
  **To add** one, copy a node line, give it a new name and position, and add
  its arrows. **To remove** one, delete the node line and the arrows that
  name it.
- **Colours** live in the palette block near the top of each file (one set
  for a light and one for a dark background); change them only there.
- **To look at the result**, run `./build.sh` and open `pdf/cneta-pipeline.pdf`
  or `pdf/cneta-overview.pdf` (the `-dark` PDFs have light text on a
  transparent background, so they look wrong in a PDF viewer; that is
  expected). `pdflatex cneta-pipeline.tex` on its own also works and gives
  the light version. The files compile unchanged on Overleaf too: upload
  the `.tex`, edit there, and copy the text back into the repository.
- **When you are done**, run `./build.sh` with poppler installed (it
  regenerates `svg/`) and commit the `.tex` and `svg/` files together. If
  you do not have poppler, commit the `.tex` only and say so in the pull
  request; the SVGs can be rebuilt by anyone who has it.

## Rebuilding

```bash
./build.sh
```

Needs `pdflatex` (a basic TeX Live or MacTeX, with TikZ and Latin Modern).
With `pdftocairo` (poppler: `brew install poppler`, `apt install
poppler-utils`) it also writes `svg/` and `png/`; without it you still get
`pdf/`. Commit the regenerated `svg/` files together with the source.
