# CNETA

CNETA is a suite of tools for analysing copy-number evolution in tumour
samples. It comprises three command-line programs:

`cnets`
: simulation of copy-number evolution.

`cnetml`
: maximum-likelihood phylogenetic inference.

`cnetmcmc`
: Bayesian (MCMC) phylogenetic inference.

:::{note}
This documentation is under active development alongside the software.
Sections marked **TODO** are placeholders awaiting content.
:::

## Where to start

- New to CNETA? Start with [Installation](installation/index.md), then the
  [Quick start](quickstart/index.md).
- Chaining the tools together? See [Workflows](workflows/index.md).
- Contributing code? See the [Developer guide](developer-guide/index.md).

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
:caption: Users
:hidden:

installation/index
quickstart/index
user-guide/index
workflows/index
file-formats/index
configuration/index
```

```{toctree}
:maxdepth: 2
:caption: Developers
:hidden:

developer-guide/index
api/cpp
```
