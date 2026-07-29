# `cnets`

`cnets` simulates copy-number evolution along a phylogenetic (coalescent)
tree of tumour samples from a single patient, and writes the simulated
copy-number profiles, tree, mutations, and timing information needed as
input to `cnetml` / `cnetmcmc`.

:::{note}
This page is migrated from the top-level `README.md` and `run-cnets.sh`.
Please review for accuracy and completeness.
:::

## Status

Mature and published; see the [citation](../index.md) on the front page.
Model 1 (bounded total copy number) is deprecated in favour of model 2
(haplotype-specific, the default) and model 3 (infinite sites).

## Usage

```bash
code/cnets [options]
```

`run-cnets.sh` is the maintained example driver — copy and edit it rather
than calling `cnets` directly. It sets the seed, verbosity, tree-generation
parameters, mutation-model parameters, and output directory in one place.

## Options

Grouped by purpose; see `run-cnets.sh` for the full set with example
values.

Tree generation
: `--epop` (effective population size, scales branch lengths to years),
  `--tdiff` (spreads tip sampling times by random multiples of this
  value), `--cons` (constrain tree height by patient age), `--age`
  (patient age at first sample).

Mutation model
: `--model` (1: bounded total copy number, deprecated; 2:
  haplotype-specific, default; 3: infinite sites), `--cn_max` (maximum
  copy number — limited by available heap space), `--method` (0:
  simulate waiting times along each branch, the default, supports all
  event types; 1: simulate sequences directly at branch ends, supports
  duplication/deletion only).

Mutation types (`--cn_type`)
: `0` segment-level only (rates r1/r2) · `1` chromosome-level + WGD only
  (r3/r4/r5) · `2` segment-level + WGD (r1/r2/r5) · `3` segment +
  chromosome level (r1-r4) · `4` all of the above (r1-r5).

Branch-specific rates
: `--bsr_mode` (0: constant rate on all branches, the default; 1: one
  shared multiplier per branch; 2: independent per-branch rate per event
  type; 3: random local clock — rates are inherited top-down, changing
  with probability `--bsr_p` at each node), `--bsr_dist` (0: log-normal,
  1: Gamma — both parameterised by `--bsr_variance`).

## Inputs

A tree file (`--tree_file`) is optional; without one, a random coalescent
tree is generated with exponential growth. See
[File formats](../file-formats/index.md) for the tree and timing file
formats.

## Outputs

Required
: `*-cn.txt.gz` (total copy number per site/sample), `*-tree.txt` /
  `*-tree.nex` (the simulated tree, calendar-time branch lengths), and
  `*-tree-nmut.nex` (branch lengths in expected mutations per site),
  `*-info.txt` / `*-mut.txt` (per-branch mutation counts and mutation
  list).

Optional
: `*-rel-times.txt` (tip sampling times), `*-haplotype-cn.txt.gz`
  (haplotype-specific copy number), `*-rcn.txt.gz` /
  `*-haplotype-rcn.txt.gz` (relative copy number), `*-inode-cn.txt.gz`
  (internal-node copy number).

## Examples

```bash
./run-cnets.sh
```

Edit the parameter block at the top of the script (tree size `Ns`,
`cn_max`, mutation rates `r1`-`r5`, `bsr_mode`, etc.) rather than passing
flags by hand. See [Quick start](../quickstart/index.md) for a walkthrough.
