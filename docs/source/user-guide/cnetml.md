# `cnetml`

`cnetml` builds phylogenetic trees from copy-number profiles by maximum
likelihood, jointly inferring tree topology, node ages, and mutation rates
from (optionally longitudinal) samples of a single patient.

:::{thumbnail} /_static/briefing.jpeg
:alt: Overview of the cnetml2 workflow, from longitudinal sampling to copy number profiles to evolutionary inference
:align: center
:class: centered-figure
:show_caption: true
:title: Longitudinal tumour sampling, copy-number profiling, and evolutionary inference with cnetml2.
:::


## Status
Mature and published — see
[Lu et al., Genome Biology 2023](https://doi.org/10.1186/s13059-023-02983-0).

## Usage
```bash
bin/cnetml [options]
```

`bin/run-cnetml.sh` is the maintained example driver, parameterised by a
`mode` variable at the top of the script.

## Modes (`mode`)
- `0` (default) — Infer the maximum-likelihood tree from input copy
  numbers. Add `-b 1` to bootstrap.
- `1` — Simple comprehensive test on a simulated tree.
- `2` — Compute the likelihood of a given tree and parameters.
- `3` — Optimize branch lengths (and rates) for a given topology. Add
  `-b 1` to bootstrap.
- `4` — Infer marginal and joint ancestral states for a given tree.

Modes 1-3 are mainly used to validate the likelihood computation, per the
original README.

## Options
Grouped by purpose; see `run-cnetml.sh` for the full set with example
values.

Time constraints
: - `--constrained`/`cons` — optimize under patient-age/tip-timing
    constraints.
  - `--estmu` — estimate mutation rates. Only reliable with informative
    tip-timing, i.e. `cons=1` and non-zero time differences between tips.

Model (`--model`/`-d`)
: - `0` — Mk, deprecated.
  - `1` — bounded total copy number, deprecated.
  - `2` — bounded haplotype-specific copy number.
  - `3` — independent Markov chain, recommended for data with
    chromosome gain/loss and WGD. Differs from the model typically used
    to *simulate* such data in `cnets` (usually model 2). Model 2 can
    also be used for reconstruction, but assumes a strict event order:
    WGD, then chromosomal gain/loss, then duplication/deletion.

Independent Markov chain dimensions (model 3 only)
: - `max_wgd` — at most 1; simulations allow at most one WGD per
    sample.
  - `max_chr_change` — usually 1, higher only if some chromosomes
    undergo more than one gain/loss event.
  - `max_site_change` — often >1, since site-level CNAs frequently
    overlap, e.g. with no WGD/chromosome events and a maximum copy
    number of 5, use `max_site_change = 3`.
  - `m_max` — maximum segment copies before chromosome-level events.

  Setting these too low gives an incorrect likelihood; too high slows
  computation. Use
  `Rscript util/check_site_pattern.R -c sim1-cn.txt -t sim1-patterns.txt`
  to help pick values for real data.

Shift-edge search (`--bsr_mode 3`, with `mode` 3)
: - `rlc_search_method` — `0` exhaustive (every subset of the candidate
    branches, small trees only), `1` stepwise (default), `2` genetic
    algorithm.
  - `rlc_criterion` — `0` AIC (default) or `1` BIC decides which set of
    shift edges is kept.
  - `estimate_bsr0_first` — `0` (default) fits the trunk's own rates
    jointly with the multipliers; `1` calibrates them first and holds
    them fixed.

  AIC and BIC count the parameters the shift edges add, K times the
  number of rate types; the node ages and reference rates fitted in
  every candidate add the same constant to each and are left out. The
  AICc in the log counts every fitted parameter and is for information
  only. BIC's sample size is the number of sites in the likelihood. The
  invariant sites are always in it: with `incl_all 1` (the default) each
  is a site of its own; with `incl_all 0` the likelihood of one all-normal
  site is multiplied by their number (the bias correction). `correct_bias`
  follows `incl_all` (0 with `incl_all 1`, 1 with `incl_all 0`); a
  `--correct_bias 0` given with `incl_all 0` is overridden with a note,
  since it would drop the invariant sites. The search prints these counts
  at its start.

Tree search (`tree_search`)
: - `0` — genetic algorithm. Slow, deprecated.
  - `1` — heuristic/hill-climbing. For at least 5 samples.
  - `2` — exhaustive. Efficient for up to 7 samples.

Mutation types (`cn_type`)
: - `0` — segment only.
  - `1` — chromosome+WGD.
  - `2` — segment+WGD.
  - `3` — segment+chromosome.
  - `4` — all.

  Which rate categories to estimate. This flag is specific to
  `cnetml`/`cnetmcmc` — `cnets` has no `cn_type` selector; it turns
  event types on/off by zeroing the corresponding rate flag instead
  (see [`cnets` options](cnets.md#options)).

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
:::{warning}
The original README (and earlier versions of this page) described these
as `*-tree.txt` / `*-tree.nex` / `*-tree.nmut.nex`. That's not what the
code does: there is no hardcoded `-tree` suffix. All of the files below
are `<ofile>` plus a fixed extension, where `<ofile>` is whatever you
pass to `-o`/`--ofile` (default `maxL-tree.txt`; `run-cnetml.sh` sets it
to `MaxL-<suffix>.txt`). See [File formats](../file-formats/index.md)
for the verified naming.
:::

`<ofile>` / `<ofile>.nex` (reconstructed tree, calendar-time branch
lengths), `<ofile>.nmut.nex` (branch lengths in mutation counts),
`<ofile>.summary.txt` (run metadata: mode, model, `bsr_mode`, likelihood,
mutation-rate estimates), `<ofile>.edge_rates.txt` (per-branch rate and
expected/observed mutation-count breakdown — only meaningful with
`bsr_mode > 0`), `*-segs.txt` (postprocessed copy-number matrix, written
when reading the input file, named by `--seg_file`).

Ancestral-state mode (`mode=4`) additionally writes `<ofile>.mrca.cn` /
`<ofile>.joint.cn` and `<ofile>.mrca.state` / `<ofile>.joint.state` for
models 0-2. Model 3 (independent Markov chains) splits the MRCA state
file by event type instead — `<ofile>.mrca.seg.state` /
`.mrca.chr.state` / `.mrca.wgd.state` — while the joint state file stays
a single `<ofile>.joint.state`. See
[File formats](../file-formats/index.md) for column definitions.

## Examples
```bash
bin/run-cnetml.sh
```

See [Quick start](../quickstart/index.md) for how this chains onto
`cnets`' output.
