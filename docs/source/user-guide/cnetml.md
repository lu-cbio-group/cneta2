# `cnetml`

`cnetml` builds phylogenetic trees from copy-number profiles by maximum
likelihood, jointly inferring tree topology, node ages, and mutation rates
from (optionally longitudinal) samples of a single patient.

:::{note}
This page is migrated from the top-level `README.md` and `run-cnetml.sh`.
Please review for accuracy and completeness.
:::

## Status

Mature and published — see
[Lu et al., Genome Biology 2023](https://doi.org/10.1186/s13059-023-02983-0).
Model 3 (independent Markov chain) is newer and described in the original
README as "in development, at test branch".

## Usage

```bash
code/cnetml [options]
```

`run-cnetml.sh` is the maintained example driver, parameterised by a
`mode` variable at the top of the script.

## Modes (`mode`)

`0` (default)
: Infer the maximum-likelihood tree from input copy numbers. Add `-b 1` to
  bootstrap.

`1`
: Simple comprehensive test on a simulated tree.

`2`
: Compute the likelihood of a given tree and parameters.

`3`
: Optimize branch lengths (and rates) for a given topology. Add `-b 1` to
  bootstrap.

`4`
: Infer marginal and joint ancestral states for a given tree.

Modes 1-3 are mainly used to validate the likelihood computation, per the
original README.

## Options

Grouped by purpose; see `run-cnetml.sh` for the full set with example
values.

Time constraints
: `--constrained`/`cons` (optimize under patient-age/tip-timing
  constraints), `--estmu` (estimate mutation rates — only reliable with
  informative tip-timing, i.e. `cons=1` and non-zero time differences
  between tips).

Model (`--model`/`-d`)
: `0` Mk, deprecated · `1` bounded total copy number, deprecated ·
  `2` bounded haplotype-specific copy number · `3` independent Markov
  chain, recommended for data with chromosome gain/loss and WGD. Model 3
  differs from the model typically used to *simulate* such data in
  `cnets` (usually model 2). Model 2 can also be used for reconstruction,
  but assumes a strict event order: WGD, then chromosomal gain/loss, then
  duplication/deletion.

Independent Markov chain dimensions (model 3 only)
: `max_wgd` (at most 1 — simulations allow at most one WGD per sample),
  `max_chr_change` (usually 1, higher only if some chromosomes undergo
  more than one gain/loss event), `max_site_change` (often >1, since
  site-level CNAs frequently overlap — e.g. with no WGD/chromosome events
  and a maximum copy number of 5, use `max_site_change = 3`), `m_max`
  (maximum segment copies before chromosome-level events). Setting these
  too low gives an incorrect likelihood; too high slows computation. Use
  `Rscript util/check_site_pattern.R -c sim1-cn.txt -t sim1-patterns.txt`
  to help pick values for real data.

Tree search (`tree_search`)
: `0` genetic algorithm — slow, deprecated · `1` heuristic/hill-climbing —
  for at least 5 samples · `2` exhaustive — efficient for up to 7 samples.

Mutation types (`cn_type`)
: as in `cnets` — `0` segment only, `1` chromosome+WGD, `2` segment+WGD,
  `3` segment+chromosome, `4` all.

## Inputs

Required
: A copy-number file (`*-cn.txt.gz` or `*-haplotype-cn.txt.gz`), space
  separated, **no header row**: `sample_ID chr_ID site_ID CN` for total
  copy number, or with an extra `cnA cnB` for haplotype-specific.
  `sample_ID` must be ordered from 1 to the number of patient samples.

Optional
: A sample-timing file (`*-rel-times.txt`); a directory of initial trees
  for tree search (recommended for large trees, to avoid local optima).

See [File formats](../file-formats/index.md) for full column definitions.

:::{warning}
"Please ensure the input file exists and its name is correct, or else
there may be an error of 'Segmentation fault (core dumped)'" — from the
original README.
:::

## Outputs

`*-tree.txt` / `*-tree.nex` (reconstructed tree, calendar-time branch
lengths), `*-tree.nmut.nex` (branch lengths in mutation counts),
`*-segs.txt` (postprocessed copy-number matrix, written when reading the
input file).

Ancestral-state mode (`mode=4`) additionally writes `*.mrca.cn` /
`*.joint.cn` and `*.mrca.state` / `*.joint.state` — see the original
README (or [File formats](../file-formats/index.md), once migrated) for
column definitions.

## Examples

```bash
./run-cnetml.sh
```

See [Quick start](../quickstart/index.md) for how this chains onto
`cnets`' output.
