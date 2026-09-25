# File formats

:::{note}
Content below is migrated from the top-level `README.md`, `mcmc.cfg`,
and example output in `example/`. Please review for accuracy — some
sections flag open questions that need confirming against the current
code.
:::

## Copy-number input

A file containing integer absolute/relative copy numbers for all the patient samples and/or the normal sample (*-cn.txt.gz or *-haplotype-cn.txt.gz).

Either compressed file or uncompressed file is fine. There need to be at least four columns, separated by space, in this file: sample_ID, chr_ID, site_ID, CN. Each column is an integer. Note that there should be __no header names__ in this file.

The sample_ID has to be __ordered from 1 to the number of patient samples__.

The chr_ID and site_ID together determine a unique site along the genome of a sample, ordering from 1 to the largest number (1, 2, 3, ...).

The site_ID can be __consecutive numbers__ from 1 to the total number of sites along the genome, or __consecutive numbers__ from 1 to the total number of sites along each chromosome of the genome.

For __haplotype-specific__ CN, there need to be at least five columns, with the last two being cnA, cnB.

If the total CN is larger than the specified maximum CN allowed by the program, the total CN will be automatically decreased to the maximum CN when the input is total CN and the program will exit when the input is haplotype-specific CN.

When the input copy numbers are relative with normal copy being 0 as those output by CGHcall, please specify it with option "--is_rcn 1".

When the input copy numbers are haplotype-specific which have been scaled relative to ploidy or not, please specify it with option "--is_total 0 --is_rcn 0".

When the input copy numbers are in bins of fixed size, please specify it with option "--bin 0" to use original data. By default "--bin 1" is used to get segment-level data by merging consecutive bins with the same copy number in a sample with change points aligned across all the samples.

## Sample timing files

`*-rel-times.txt` (optional input to `cnetml`/`cnetmcmc`; also written
by `cnets` as simulation output) records the sampling time of each tip
node.

Three space-separated columns, no header:

| column      | meaning                                                  |
| ----------- | --------------------------------------------------------- |
| `sample_ID` | ordered from 1 to n (number of patient samples)          |
| `time`      | time relative to the 1st sample, in years (float)        |
| `age`       | patient age at the time of sampling, in years (integer)  |

Providing this file lets `cnetml`/`cnetmcmc` use tip-time information
when estimating divergence times and mutation rates (see `--constrained`
in `cnetml`, `--cons` in `cnetmcmc`, and `--tdiff` in `cnets` — see
[Configuration](../configuration/index.md)).

## Tree formats

:::{warning}
Filenames below were checked against `code/cnets.cpp`, `code/cnetml.cpp`,
and each tool's `run-*.sh` — not just taken from `README.md`, which gets
`cnetml`'s naming wrong (see the `cnetml` subsection).
:::

### `cnets` (tab-delimited `*-tree.txt`, NEWICK `*-tree.nex` / `*-tree-nmut.nex`)

Filenames are hardcoded as `<prefix>-tree.txt`, `<prefix>-tree.nex`,
`<prefix>-tree-nmut.nex` — not user-configurable beyond the run's
`prefix`.

The tab-delimited file has one row per branch (edge) of the tree,
columns `start`, `end`, `length`, `eid`, `nmut`:

| column   | meaning                                                |
| -------- | ------------------------------------------------------- |
| `start`  | parent node ID                                         |
| `end`    | child node ID                                          |
| `length` | branch length, in years                                |
| `eid`    | edge ID (matches the `edge_ID` column in `*-mut.txt`)  |
| `nmut`   | number of mutations simulated on this branch           |

:::{warning}
The top-level `README.md` states node IDs run 0..n-1 for samples, n
for the normal sample, n+1 for the root, and n+2..2n+1 for internal
nodes. The example output in `example/sim-data-1-tree.txt` instead
uses 1-based sample IDs (1..4, with 5 as root). Please confirm which
convention the current code actually uses before relying on this.
:::

The two NEWICK files are standard NEXUS-wrapped NEWICK, one tree per
file:

- `*-tree.nex`: branch length = calendar time (years)
- `*-tree-nmut.nex`: branch length = expected number of CNAs per site

Example (from `example/sim-data-1-tree.nex`):

```text
#nexus
begin trees;
tree 1 = (((2:13.5854296603,3:15.5854296603)6:3.8658646437,1:17.4512943040)7:0.8285947280,4:0.0000000000)5;
end;
```

### `cnetml` (`<ofile>`, `<ofile>.nex`, `<ofile>.nmut.nex`)

:::{warning}
`README.md` (and earlier versions of the `cnetml` User guide page)
describe these as `*-tree.txt` / `*-tree.nex` / `*-tree.nmut.nex`. That's
wrong: `code/cnetml.cpp` has no hardcoded `-tree` suffix anywhere.
:::

`cnetml` writes its reconstructed tree to whatever file you pass with
`-o`/`--ofile` (C++ default: `maxL-tree.txt`), then appends fixed
extensions for the other outputs. `run-cnetml.sh` sets `--ofile` to
`MaxL-<suffix>.txt`, so in practice you'll see e.g. `MaxL-sim1.txt`,
`MaxL-sim1.txt.nex`, `MaxL-sim1.txt.nmut.nex` — not files with `-tree`
in the name at all, unless you choose an `--ofile` value that includes
it.

- `<ofile>`: the tree in the same tab-delimited format as `cnets`
  (calendar-time branch lengths)
- `<ofile>.nex`: NEWICK, calendar-time branch lengths
- `<ofile>.nmut.nex`: NEWICK, branch length = number of mutations

### `cnetml` run reports (`<ofile>.summary.txt`, `<ofile>.edge_rates.txt`)

Written alongside the tree files above; not mentioned in `README.md` at
all.

`<ofile>.summary.txt`
: `key<TAB>value` lines: `mode`, `model`, `bsr_mode`, `cn_type`,
  `constrained`, `estmu`, `estimate_bsr0_first`, `rlc_criterion` (for
  `bsr_mode 3`), `raw_logL`, `penalized_score`, `K_shift_edges`,
  `shift_eids`, per-event-type reference mutation rates
  (`dup_rate_reference`, ...), length-weighted effective rates
  (`dup_rate_final_global_weighted`, ...), `n_edges`, `n_leaves`.

`<ofile>.edge_rates.txt`
: One row per edge. Columns: `eid`, `start`, `end`, `length`,
  `nmut_expected`, `dup_expected`, `del_expected`, `chr_gain_expected`,
  `chr_loss_expected`, `wgd_expected`, `nmut_observed`, `dup_observed`,
  `del_observed`, `chr_gain_observed`, `chr_loss_observed`,
  `wgd_observed`, `dup_rate`, `del_rate`, `chr_gain_rate`,
  `chr_loss_rate`, `wgd_rate`, `total_rate`, `m_shared`, `m_dup`,
  `m_del`, `m_chr_gain`, `m_chr_loss`, `m_wgd`, `is_shift_edge`,
  `local_clock_id`. The last several columns are only meaningful with
  `bsr_mode > 0` (`NA` otherwise).

### MCMC sampled trees (`*.t`)

Written by `cnetmcmc` in a format similar to MrBayes output — a block
of trees sampled during the chain. Analyze with
[TreeAnnotator](https://beast.community/treeannotator) to get a
maximum-clade-credibility summary tree.

## Mutation outputs

### Simulated mutation list (`*-mut.txt`, from `cnets`)

One row per simulated mutation event. Columns, in the order the code
writes them (`code/cnets.cpp`, the `*-mut.txt` writer): `sample_ID`,
`edge_ID`, `muttype_ID`, `mut_btime`, `mut_etime`, `chr_haplotype`,
`chr`, `seg_ID`.

- `edge_ID` matches the `eid` column in `*-tree.txt` (see
  [Tree formats](#tree-formats)), so mutations can be mapped onto
  branches.
- For a chromosome gain/loss event, `seg_ID` is `-1`.
- For a whole-genome-doubling event, `chr` is `0` and `seg_ID` is `-1`.

:::{note}
`muttype_ID` isn't enumerated anywhere in the code near the writer;
treat its values as opaque until confirmed (likely the same event-type
IDs as `--cn_type`/`r1`-`r5`).

`chr_haplotype` and `chr` **are** resolvable, and don't mean quite what
the column names suggest. Internally each chromosome is stored twice —
once per haplotype copy — as a single 0-based index 0..43
(`genome.hpp`/`genome.cpp`: `hap_index = chromosome(0..21) + 22 *
haplotype(0 or 1)`, since `NUM_CHR = 22`, autosomes only). The mutation
struct keeps only this combined index (`mutation.chr` — there's no
separate haplotype field). The writer then outputs:

- `chr_haplotype` = `mutation.chr + 1` (range 1..44) — chromosome *and*
  haplotype together, not haplotype alone despite the column name.
- `chr` = `(mutation.chr + 1) % 22` — recovers the 1-based chromosome
  number, **except** chromosome 22 itself wraps to `0` (there's no
  explicit re-mapping back to 22).

To recover haplotype separately: `mutation.chr >= 22` (i.e.
`chr_haplotype > 22`) is the second haplotype copy.
:::

### Per-branch mutation counts (`*-info.txt`, from `cnets`)

The time of each node and the total number of mutations simulated on
each branch, grouped by the lineages of tip nodes.

### Simulation edge-rate table (`*-edge_rates.txt`, from `cnets`)

Same column layout as `cnetml`'s
[`<ofile>.edge_rates.txt`](#tree-formats) — the true, simulation-time
per-branch rates and (for `bsr_mode > 0`) which branches got a rate
shift, rather than `cnetml`'s reconstructed estimates. Not mentioned in
`README.md`.

### Ancestral-state reconstruction outputs (`cnetml` mode 4)

Named `<ofile>.<suffix>`, same `<ofile>` as in
[Tree formats](#tree-formats).

For models 0-2:

| file             | contents                                                                                                                                                                                      |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `<ofile>.mrca.cn`    | Reconstructed copy numbers for the most recent common ancestor of all tumour samples. `node_ID chr_ID site_ID CN` (total) or `node_ID chr_ID site_ID cnA cnB` (haplotype-specific).           |
| `<ofile>.joint.cn`   | Same, for all internal nodes.                                                                                                                                                                 |
| `<ofile>.mrca.state` | Posterior probability of each possible copy-number state, for each variant site (a site with at least one atypical copy number across samples), at the MRCA node: `node_ID site_ID probability_stateID`. |
| `<ofile>.joint.state`| Possible copy-number state at each variant site, for all internal nodes: `node_ID site_ID cn_stateID`.                                                                                        |

In `<ofile>.mrca.state`/`<ofile>.joint.state`, `site_ID` is written as
`chromosomeID_siteID`, and the state ID matches the state index of the
Markov-model rate matrix used for reconstruction.

:::{note}
Model 3 (independent Markov chains) writes `<ofile>.mrca.cn` the same
way, but splits the MRCA state file by event type instead of one
`.mrca.state`: `<ofile>.mrca.seg.state` (only if segment-level events
are in play), `<ofile>.mrca.chr.state` (chromosome-level), and
`<ofile>.mrca.wgd.state` (whole-genome doubling) — each written only
when that event type has a non-zero state-space dimension. The joint
file stays a single `<ofile>.joint.state`, same format as above. None of
this is in `README.md`; model 3 is still described there as
"in development, at test branch".
:::

## MCMC configuration and trace files

### Config file (`mcmc.cfg`)

Plain-text `key=value` pairs, one per line, `#` for comments and blank
lines allowed. Passed to `cnetmcmc` with `--config_file` (see
`run-cnetmcmc.sh`). Grouped in the file by purpose:

```text
verbose=0

# Parameters about input files
is_bin=0
incl_all=0

# Model parameters
cn_max=4
model=2
cons=0

# MCMC parameters
n_draws=2000
n_burnin=1000
n_gap=10   # sampling every kth sample

# Parameters for initial coalescence tree
epop=90000
beta=1.563e-3
gtime=0.002739726
```

...followed by blocks for the true/initial mutation rates, per-rate-type
prior/proposal parameters (`sigma_l*` for the log-normal prior,
`sigma_*` for the proposal, one pair per event type), tree-height
parameters, and branch-length parameters. See `mcmc.cfg` in the
repository root for the full set with inline comments, and
[`cnetmcmc` options](../user-guide/cnetmcmc.md#options) for what each
group controls.

### Trace files (output)

MrBayes-compatible:

`*.p`
: Parameter traces. Check convergence with
  [Tracer](https://beast.community/tracer) or `util/check_convergence.R`
  (uses [RWTY](https://github.com/danlwarren/RWTY)).

`*.t`
: Sampled trees — see [Tree formats](#tree-formats).

## Identifiers and conventions

`sample_ID`
: 1-based, ordered from 1 to n (the number of patient samples). Used
  consistently in `*-cn.txt.gz`/`*-haplotype-cn.txt.gz` and
  `*-rel-times.txt`. The normal sample is not counted in n and is
  identified separately (see the tree-node-ID warning under
  [Tree formats](#tree-formats) — the exact convention needs
  confirming).

`chr_ID`, `site_ID`
: Both 1-based. Together they identify a unique site on the genome of a
  sample. `site_ID` can be numbered either consecutively across the
  whole genome, or consecutively within each chromosome — the file
  itself doesn't say which, so the two ends of a pipeline (e.g. a
  preprocessing script and `cnetml`) must agree on the convention out of
  band.

`edge_ID` / `eid`
: 1-based, shared between `*-tree.txt` and `*-mut.txt` to map mutations
  onto tree branches.

0- vs 1-based coordinates
: Input files (`*-cn.txt.gz`, `*-rel-times.txt`) are consistently
  1-based. Tree node IDs are where the documentation and an actual
  example output disagree — see the warning under
  [Tree formats](#tree-formats).
