# Developer guide

```{toctree}
:maxdepth: 1

architecture
documentation
```

## Repository layout

```text
cneta/
├── bin/           built executables (created by the build, not in git)
├── code/          C++ sources for cnets, cnetml, cnetmcmc
│   ├── gzstream/  vendored: gzip streams
│   ├── lbfgsb/    vendored: L-BFGS-B optimiser
│   └── matexp/    vendored: matrix exponential
├── tests/         test suite — see tests/README.md
│   ├── unit/      Catch2 unit tests (C++)
│   ├── e2e/       end-to-end tests (pytest)
│   └── data/      committed fixtures
├── ilp/           integer-programming experiments
├── util/          R and Python helper scripts
├── docs/          this documentation site
└── run-*.sh       example driver scripts
```

## Build system

`code/CMakeLists.txt` is the root of the CMake project. It builds one static
library, `cneta_core`, holding all the shared science, and links each of the
three thin command-line programs against it. `cneta_core` is the de-facto
precursor to the planned `libcneta`.

The three executables are written to `bin/` at the repository root, wherever
the build tree happens to live. Pass
`-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=<path>` to send them somewhere else — CI and
the build-system tests do this so a throwaway build does not overwrite the
binaries in your working tree.

`code/build.sh` is a thin wrapper: it optionally sources an HPC environment
script from `code/envs/`, then configures and builds in `code/build/`. See
[Installation](../installation/index.md).

`code/makefile` is the pre-CMake build. It compiles the sources directly and
puts the executables in `code/` rather than `bin/`. It is kept for reference
but is not maintained — use CMake.

The project builds as C++14. That is the minimum Catch2 v3 requires, and the
sources compile cleanly at that level. C++17 additionally needs the
`std::random_shuffle` calls in `tree_op.cpp` and `genome.cpp` replaced with
`std::shuffle`, which changes the order random numbers are drawn in and so
changes fixed-seed simulation output — a deliberate change to make on its own,
not a side effect.

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

The suite lives in `tests/` and has three layers:

`tests/unit/`
: Catch2 unit tests for the shared library. They link `cneta_core` directly
  and cover the state encoding, the rate and transition matrices, and the
  tree helpers. Built only with `-DCNETA_BUILD_TESTING=ON`.

`tests/e2e/`
: pytest tests that run the compiled programs on a small simulated dataset,
  checking output files, error handling, and reproducibility. Includes smoke
  tests that execute the `run-*.sh` scripts as a user would.

`tests/build_system/`
: pytest tests that configure and compile the project from scratch and check
  `build.sh`'s argument handling.

To run everything:

```bash
cmake -S code -B code/build -DCNETA_BUILD_TESTING=ON
cmake --build code/build -j 4
ctest --test-dir code/build --output-on-failure
pip install -r tests/requirements.txt
OMP_NUM_THREADS=1 pytest tests -m "not build"
```

`OMP_NUM_THREADS=1` matters: `cnetml` parallelises its tree search with
OpenMP, and floating-point sums over a varying number of threads are not
bit-reproducible.

The regression anchor is `cnetml --mode 2`, which scores a supplied tree with
no RNG and no search and so is reproducible to the digit. Reference values
live in `tests/data/expected.json`. A change there means the likelihood
calculation changed — which is exactly what needs watching during the
`libcneta` refactor.

`tests/README.md` has the full detail, including what is *not* covered yet and
a list of known failures the suite records but does not fail on.

## Continuous integration

Workflows live in `.github/workflows/`. Currently:

`docs.yml`
: builds the documentation on every pull request into `main`, and
  publishes to GitHub Pages when a pull request is merged.

`ci.yml`
: builds the project and runs the unit and end-to-end tests on
  `ubuntu-latest` and `macos-latest`.

Both run on pull requests into `main` and on pushes to `main` — not on every
pushed commit, so work-in-progress on a branch with no open pull request
costs nothing. Both can also be triggered by hand from the Actions tab
(`workflow_dispatch`).

TODO: container workflow.

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
