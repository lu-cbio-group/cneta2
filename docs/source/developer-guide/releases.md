# Versions and releases

Every result CNETA produces should be traceable to the code that produced
it. That needs one version number, recorded wherever results are, and a
simple routine for making releases. This page describes both.

## The version number

CNETA uses `MAJOR.MINOR.PATCH`, starting from `2.0.0`. What each part means
here:

| Change | Example | When |
|---|---|---|
| **PATCH** | `2.1.0` to `2.1.1` | Fixes that change no result: crashes, error messages, documentation, build and test changes. The same input and seed give the same output. |
| **MINOR** | `2.1.0` to `2.2.0` | New features that leave existing command lines working (new options, models, output files), **and any change that alters results** for the same input and seed: a corrected likelihood, a different random-number order, an optimiser change. |
| **MAJOR** | `2.x` to `3.0.0` | Changes that break existing use: an option removed or renamed, an input or output format changed, a program removed; later also an incompatible change to the library or Python interface. |

Changes that alter results go into MINOR rather than PATCH on purpose. Then
the rule is simple: **any `2.1.x` reproduces a result made with any other
`2.1.x`**, given the same input and seed. A paper can cite `2.1` and remain
reproducible across later bug-fix releases.

## Where it lives

- **`VERSION`** at the repository root holds the number, nothing else. It is
  the only place to change it. The CMake build and the legacy `makefile`
  pass it to the programs, the documentation reads it, and a test checks its
  form.
- **The programs** print it with `--version`, print a `version:` line at the
  start of every run, and record it in their outputs next to the seed:
  `VERSION:` at the top of `cnets`' `<prefix>-info.txt`, `version` in
  `cnetml`'s `.summary.txt`.
- **The commit** is added when the build knows it, for example
  `2.1.0 (git v2.1.0-5-gabc1234)`: five commits after the `v2.1.0` tag.
  `-dirty` means the build had uncommitted changes; don't publish results
  from such a build. The commit is read when CMake configures; `build.sh`
  configures every time.
- **`CITATION.cff`** at the repository root holds the citation data: the
  paper to cite for now and the software's authors. It carries no version
  number, so a release does not touch it; see "Keeping the citation data
  up to date" below.
- **Git tags** `vX.Y.Z` mark the commit of each release, and a **GitHub
  Release** made from the tag carries the list of changes. There is no
  changelog file to maintain: GitHub writes the list from the titles of
  the pull requests merged since the previous release ("Generate release
  notes"), so a PR title should say what changed, and a PR that changes
  results should say so in its description.

Between releases `VERSION` keeps the last released number; the commit part
of the version string tells later builds apart.

## Making a release

1. Decide the number from the pull requests merged since the last tag
   (`git log --oneline vX.Y.Z..main`, or the preview that "Generate release
   notes" shows): a new feature or any change that alters results means at
   least MINOR; only fixes that change no result mean PATCH.
2. In one PR, set the number in `VERSION`.
3. Merge, and wait for CI to pass on `main`.
4. Tag the merge commit and push the tag:
   ```bash
   git tag -a vX.Y.Z -m "CNETA X.Y.Z"
   git push origin vX.Y.Z
   ```
   (or create the tag in GitHub's "Draft a new release" form, on `main`).
5. Create a GitHub Release from the tag, click "Generate release notes",
   and add a line for anything that changes results for the same input
   and seed.

## Results for papers

- Produce them with a released build: `--version` shows no
  `-N-g<commit>` part and no `-dirty`.
- Report the CNETA version and the seeds. Both are in the output files and
  the logs; the run scripts also keep `seed-history.txt`.
- Cite the CNETML paper (README, "Citation"). Once the repository is
  connected to Zenodo, each GitHub Release gets a DOI, so the software
  release can be cited as well.

## Keeping the citation data up to date

`CITATION.cff` has two parts that change at different times:

- **`preferred-citation`** is the paper users are asked to cite; GitHub
  shows it under "Cite this repository". When a new paper is published,
  replace this block (title, authors, journal, volume, pages, year, DOI)
  and the reference in the README's "Citation" section.
- **`authors`** are the people who wrote the software. Zenodo credits them
  when a release is archived; it does not read the paper's authors. When
  someone joins the development, add a line with their name (and ORCID,
  if they have one).
- **Once Zenodo is connected**, add the "all versions" DOI as `doi:` in
  `CITATION.cff` and in the README's "Citation" section, so the software
  itself can be cited; each release's own DOI comes from Zenodo.

## Recommended next steps (not set up yet)

- **Zenodo**: connect the GitHub repository, so that every release is
  archived and gets a DOI, plus one DOI for "all versions"; then follow
  "Keeping the citation data up to date".
- **A CI check** that a pushed `vX.Y.Z` tag matches `VERSION`.
