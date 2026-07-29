# `cnetmcmc`

`cnetmcmc` builds phylogenetic trees from copy-number profiles using a
Bayesian MCMC approach.

:::{note}
This page is migrated from the top-level `README.md`, `run-cnetmcmc.sh`,
and `mcmc.cfg`. Please review for accuracy and completeness.
:::

## Status

Under active development. Only the basic MCMC algorithm is implemented,
and it has not been comprehensively tested (per the original README) — the
least mature of the three tools. No workflow or architecture diagrams
exist for it yet; see [Workflows](../workflows/index.md).

## Usage

```bash
code/cnetmcmc [options]
```

`run-cnetmcmc.sh` runs multiple chains and reads most parameters from
`mcmc.cfg` via `--config_file`.

## Modes

Two modes depending on whether a reference tree is supplied (`--rtree`):
with one, the topology is fixed and only branch lengths/rates are sampled;
without one, topology is sampled too.

## Options

`--rtree`
: Optional reference tree; fixes the topology if given.

`--init_tree`
: `0` random tree · `1` a provided tree (`--file_itree`) · `2` a random
  tree sharing the real tree's topology.

`mcmc.cfg` (via `--config_file`)
: MCMC control (`n_draws`, `n_burnin`, `n_gap`, `sample_prior`,
  `fix_topology`); proposal/prior parameters for mutation rates
  (`sigma_l*` for the log-normal prior, `sigma_*` for the proposal, per
  event type); tree-height and branch-length priors (`sigma_height`,
  `tlen_shape`, `tlen_scale`, `dirichlet_param`/`dirichlet_alpha`); and
  proposal step sizes (`lambda`, `lambda_all`).

## Inputs

The same copy-number and sample-timing files as `cnetml` — see
[File formats](../file-formats/index.md).

## Outputs

MrBayes-compatible trace files:

`*.p`
: Parameter traces. Check convergence with
  [Tracer](https://beast.community/tracer) or `util/check_convergence.R`
  (uses [RWTY](https://github.com/danlwarren/RWTY)).

`*.t`
: Sampled trees. Summarise with
  [TreeAnnotator](https://beast.community/treeannotator) into a maximum
  clade credibility tree.

## Examples

```bash
./run-cnetmcmc.sh
```
