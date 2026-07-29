# Developer guide

```{toctree}
:maxdepth: 1

architecture
documentation
```

## Repository layout

```text
cneta/
├── code/          C++ sources for cnets, cnetml, cnetmcmc
│   ├── gzstream/  vendored: gzip streams
│   ├── lbfgsb/    vendored: L-BFGS-B optimiser
│   └── matexp/    vendored: matrix exponential
├── ilp/           integer-programming experiments
├── util/          R and Python helper scripts
├── docs/          this documentation site
└── run-*.sh       example driver scripts
```

## Build system

TODO: `code/CMakeLists.txt`, `build.sh`, and the legacy `makefile`.

## Architecture

See [Architecture](architecture.md) for how the main files fit together and
how `cnetml` computes and optimizes the DECOMP likelihood.

TODO: the planned `libcneta` core.

## Testing

TODO: CLI smoke tests, Catch2 unit tests, regression datasets.

## Continuous integration

Workflows live in `.github/workflows/`. Currently:

`docs.yml`
: builds the documentation on every pull request into `main`, and
  publishes to GitHub Pages when a pull request is merged.

TODO: build, test, and container workflows.

## Contributing

TODO: branching model, review expectations, commit conventions.
