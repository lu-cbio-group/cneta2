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

`code/CMakeLists.txt` is the build of record (CMake >= 3.10, C++11). It
`find_package`s Boost (`program_options`, `filesystem`), GSL, zlib, and
optionally OpenMP; builds `lbfgsb` as a CMake subdirectory and vendored
`gzstream` as a small library; links both into a shared `cneta_core`
static library (all the evolutionary-model, tree, and I/O sources); and
links `cneta_core` into the three executables `cnets`, `cnetml`,
`cnetmcmc`.

`code/build.sh` wraps that CMake configure/build into one command:

```bash
./build.sh              # no arguments: shows help
./build.sh local 4      # configure + build in code/build/, 4 cores
./build.sh Eureka2 8    # source code/envs/setup_env_Eureka2.sh first, then build
./build.sh clean        # rm -rf code/build/
```

For a `target` other than `local`, it sources
`code/envs/setup_env_<target>.sh` (module loads for that cluster) before
running `cmake -Wno-dev ..` and `make -j<cores>` in `code/build/`. See
[Installation](../installation/index.md) for the per-platform walkthrough.

`code/makefile` is a legacy, pre-CMake build (hardcoded compiler/library
paths edited by hand at the top of the file, no automatic dependency
discovery). It predates `CMakeLists.txt` and is not the maintained path,
but is still checked in. It also has its own OpenMP switch, separate
from the CMake build's automatic detection: set `omp =` (empty) at the
top of the file to turn OpenMP off.

## Architecture

See [Architecture](architecture.md) for how the main files fit together,
how `cnetml` computes and optimizes the DECOMP likelihood, and what the
planned `libcneta` core would change.

## Testing

TODO: CLI smoke tests, Catch2 unit tests, regression datasets.

:::{note}
`code/build/` currently contains some CTest/Catch2-looking generated
files (`Testing/`, `tests-unit/`), but there is no test source directory
and no test target in `CMakeLists.txt` — they don't appear to come from
an active test suite in this repo. Please confirm with whoever generated
that build directory before assuming any tests exist.
:::

## Continuous integration

Workflows live in `.github/workflows/`. Currently:

`docs.yml`
: builds the documentation on every pull request into `main`, and
  publishes to GitHub Pages when a pull request is merged. Since GitHub
  Pages has no per-PR preview environment, the PR build instead posts a
  downloadable HTML artifact as a stand-in — download it from the
  checks tab and open `index.html` locally to preview the change.

TODO: build, test, and container workflows — none exist yet.

## Contributing

**Workflow**: branch from `main` → make your change → open a pull
request into `main` → get it reviewed if the change warrants it → merge.
For a docs change, `make -C docs strict` locally is the same check CI
runs, so a local pass means CI will pass too — see
[Working on the documentation](documentation.md).

**Branch naming**: short and descriptive, e.g. `fix-cnetml-options-page`
or `add-slurm-tutorial`.

```bash
git checkout -b <branch-name>

# ... make your changes ...

git status                             # see what changed
git add path/to/file                   # add specific files — avoid
                                        # `git add .`, it can sweep in
                                        # docs/build/ output or other
                                        # local files by accident
git diff --staged                      # review what will actually commit
git commit -m "Brief description of your changes"
git push -u origin <branch-name>
```

Then open a pull request into `main` on GitHub (a GUI git client works
just as well for these steps).

There's no issue tracker linked to branches/PRs yet — for now, just
branch and go.
