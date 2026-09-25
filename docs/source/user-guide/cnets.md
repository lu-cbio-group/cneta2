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
  value), `--constrained` (constrain tree height by patient age —
  `run-cnets.sh` calls this `cons`, but the actual flag is
  `--constrained`; `cnetmcmc` is the one tool where `--cons` is real,
  see its [Options](cnetmcmc.md#options)), `--age` (patient age at
  first sample).

Mutation model
: `--model` (1: bounded total copy number, deprecated; 2:
  haplotype-specific, default; 3: infinite sites), `--cn_max` (maximum
  copy number — limited by available heap space), `--method` (0:
  simulate waiting times along each branch, the default, supports all
  event types; 1: simulate sequences directly at branch ends, supports
  duplication/deletion only).

Event rates
: `--dup_rate` (r1), `--del_rate` (r2) — site-level duplication/deletion;
  `--chr_gain` (r3), `--chr_loss` (r4) — chromosome-level gain/loss;
  `--wgd` (r5) — whole-genome doubling. Set any of these to `0` to
  exclude that event type from the simulation.

:::{note}
`cnets` has no `--cn_type` selector — that flag (and the r1-r5 grouping
by `cn_type` value) belongs to `cnetml`/`cnetmcmc`, which use it to pick
which rate categories to *estimate* during tree building. An earlier
version of this page borrowed that description for `cnets` by mistake;
`code/cnets.cpp` has no `cn_type` anywhere. For `cnets`, which event
types get simulated is controlled directly by which of the five rate
flags above are non-zero.
:::

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
  list), `*-edge_rates.txt` (true per-branch rates used in the
  simulation — not in the original README; see
  [File formats](../file-formats/index.md)).

Optional
: `*-rel-times.txt` (tip sampling times), `*-haplotype-cn.txt.gz`
  (haplotype-specific copy number), `*-rcn.txt.gz` /
  `*-haplotype-rcn.txt.gz` (relative copy number), `*-inodes-cn.txt.gz` /
  `*-inodes-haplotype-cn.txt.gz` (internal-node copy number, total and
  haplotype-specific).

:::{note}
The README spells this `*-inode-cn.txt.gz` (singular); the code writes
`*-inodes-cn.txt.gz` (plural). The haplotype-specific variant isn't
mentioned in the README at all.
:::

## Examples

```bash
./run-cnets.sh
```

Edit the parameter block at the top of the script (tree size `Ns`,
`cn_max`, mutation rates `r1`-`r5`, `bsr_mode`, etc.) rather than passing
flags by hand. See [Quick start](../quickstart/index.md) for a walkthrough.
