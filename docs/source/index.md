# CNETA

CNETA is a suite of tools for analysing copy-number evolution in tumour
samples. It comprises three command-line programs:

`cnets`
: simulation of copy-number evolution.

`cnetml`
: maximum-likelihood phylogenetic inference.

`cnetmcmc`
: Bayesian (MCMC) phylogenetic inference.

This is version 2 of CNETA
([lu-cbio-group/cneta2](https://github.com/lu-cbio-group/cneta2)); the
first version is at [ucl-cssb/cneta](https://github.com/ucl-cssb/cneta).

:::{note}
This documentation is under active development alongside the software.
Sections marked **TODO** are placeholders awaiting content.
:::

## Where to start

- New to CNETA? Start with [Installation](installation/index.md), then the
  [Quick start](quickstart/index.md).
- Chaining the tools together? See [Workflows](workflows/index.md).
- Contributing code? See the [Developer docs](developer-guide/index.md).

## At a glance

<p class="only-light" align="center"><a href="_static/cneta-overview.svg" target="_blank"><img src="_static/cneta-overview.svg" alt="CNETA at a glance" width="100%"></a></p>
<p class="only-dark" align="center"><a href="_static/cneta-overview-dark.svg" target="_blank"><img src="_static/cneta-overview-dark.svg" alt="CNETA at a glance" width="100%"></a></p>

What CNETA does: input, the two inference programs, the result, and the
event types and rates they share.

## Project status

`cnets` and `cnetml` are the mature, published tools; `cnetmcmc` is under
active development and not yet comprehensively tested. This documentation
tracks CNETA version 2.0.

If you use CNETA, please cite:

> Lu B, Curtius K, Graham TA, Yang Z, Barnes CP. CNETML: maximum likelihood
> inference of phylogeny from copy number profiles of multiple samples.
> *Genome Biol* 24, 144 (2023).
> [doi:10.1186/s13059-023-02983-0](https://doi.org/10.1186/s13059-023-02983-0)

```{toctree}
:maxdepth: 2
:hidden:

Home <self>
users/index
developer-guide/index
```
