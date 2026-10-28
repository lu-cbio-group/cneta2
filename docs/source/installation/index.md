# Installation

`cneta` is written in C++ and built with CMake; a few R and Python scripts
under `util/` support pre/post-processing but are not required to build
the main tools. Everything needed comes from conda environment files in the
repository, so the same commands work on macOS, Linux, HPC clusters with
conda, and Windows via WSL. See the [Developer docs](../developer-guide/index.md)
for the current architecture and what's still planned.

## Environment files

| File | Environment | What it is for | Who needs it |
|---|---|---|---|
| `code/environment-mac.yml` | `cneta-build` | Build and run `cnets` / `cnetml` / `cnetmcmc` on macOS | Everyone on macOS |
| `code/environment-linux.yml` | `cneta-build` | Build and run `cnets` / `cnetml` / `cnetmcmc` on Linux, HPC with conda, Windows via WSL (Ubuntu) | Everyone on Linux |
| `util/environment.yml` | `cneta-util` | Plotting and other R/Python scripts in `util/` (same file on macOS and Linux) | Anyone who plots or post-processes results |
| `docs/environment.yml` | `cneta-docs` | Build the documentation website | Website maintainers only |

You need a conda installation, for example
[Miniforge](https://github.com/conda-forge/miniforge).

## Which system are you on?

| System | Steps |
|---|---|
| macOS | Install the Xcode Command Line Tools (`xcode-select --install`), then use `code/environment-mac.yml` |
| Linux workstation | Use `code/environment-linux.yml` |
| HPC with conda | Use `code/environment-linux.yml` |
| HPC without conda | `cd code && ./build.sh <cluster> <cores>` with `code/envs/setup_env_<cluster>.sh` (module loads); see the HPC tab below |
| Windows | Install WSL with Ubuntu, then follow the Linux steps inside Ubuntu |

## Building from source

Clone the repository:

```bash
git clone https://github.com/lu-cbio-group/cneta2.git
cd cneta2
```

Create the build environment and compile for your platform:

::::{tab-set}
:::{tab-item} macOS

```bash
xcode-select --install          # once, if not already installed
conda env create -f code/environment-mac.yml
conda activate cneta-build
cd code && ./build.sh local 4   # 4 = number of cores
```

Works on Apple Silicon and Intel.
:::

:::{tab-item} Linux / HPC with conda

```bash
conda env create -f code/environment-linux.yml
conda activate cneta-build
cd code && ./build.sh local 4   # 4 = number of cores
```

Why the compiler and `sysroot` versions are pinned is explained at the top
of `code/environment-linux.yml`: in short, so the build works and the
programs also run on older clusters (e.g. CentOS 7, RHEL 8).
:::

:::{tab-item} Windows (WSL)

Install WSL with Ubuntu (in PowerShell: `wsl --install -d Ubuntu`), open
Ubuntu, install conda there, and follow the **Linux / HPC with conda** tab
inside Ubuntu.
:::

:::{tab-item} HPC without conda

On most High Performance Clusters dependencies are provided as environment modules.
The build script, `./build.sh <cluster> <cores>`, sources `code/envs/setup_env_<cluster>.sh`
and loads the modules for you automatically; no need to source it yourself first.

Once your `code/envs/setup_env_<cluster>.sh` script exists, you can build the code like so:

```bash
cd code && ./build.sh Eureka2 <cores>
```

Eureka2, local CPU cluster at Surrey serves as an example here. The script loads this module set:

```text
GCC/13.3.0
CMake/3.29.3
Boost/1.85.0
GSL/2.8
```

To add another cluster, add `code/envs/setup_env_<name>.sh` following the
same pattern and run `./build.sh <name> <cores>`.

TODO: scheduler-specific notes (queue/partition names, resource-request
conventions for `sbatch`/`srun`) — not yet documented here.
:::
::::

`build.sh` detects the active conda environment and takes Boost, GSL, zlib
and OpenMP from it, so libraries installed elsewhere (e.g. Homebrew) are not
picked up.

`build.sh` configures a CMake build in `code/build/`. The three executables
(`cnets`, `cnetml`, `cnetmcmc`) land in `bin/` at the top of the repository,
next to the `run-*.sh` driver scripts, which is where the scripts and the test
suite look for them. The compiled programs are not tracked by git; the scripts
are. Set `CNETA_BIN` if you want the run scripts to use a different build.

With only `cneta-build` active, the run scripts work with `plot=0` (the
default). For plots, see [Running with plots](#running-with-plots).

Other `build.sh` usage:

```bash
./build.sh              # no arguments: shows help
./build.sh clean        # removes code/build/
```

### OpenMP

`cnetml` parallelises its tree search with OpenMP. The configure step reports
whether it was found:

```text
-- OpenMP: enabled (...)
```

or, if it was not:

```text
-- OpenMP: not found, cnetml will run single-threaded (macOS: brew install libomp)
```

Without OpenMP the build still succeeds and runs the same algorithm, just
without the speed-up.

- **conda (`cneta-build`):** both environment files include OpenMP
  (`llvm-openmp`), so nothing extra is needed.
- **HPC without conda:** GCC includes OpenMP.
- **Thread count:** OpenMP uses all available cores by default. Set
  `OMP_NUM_THREADS` to limit it, e.g. `OMP_NUM_THREADS=4 bin/run-cnetml.sh`.
  Floating-point sums over a varying number of threads are not
  bit-reproducible, so pin `OMP_NUM_THREADS=1` when you need exactly
  repeatable scores (the test suite does).

## Plotting and `util/` scripts (optional)

The R and Python scripts in `util/` (see
[Utility scripts](../user-guide/utility-scripts.md)) have their own
environment, the same file on macOS and Linux:

```bash
conda env create -f util/environment.yml
conda activate cneta-util
Rscript util/install_cran_packages.R   # only needed for compare_trees.R and check_convergence.R
```

`install_cran_packages.R` installs `Quartet` and `rwty` from CRAN (they are
not on conda); it skips packages that are already installed and prints
`OK`/`FAILED` for each.

None of this is needed to build or run `cnets`, `cnetml`, or `cnetmcmc`
themselves.

### Running with plots

With `plot=1` in `bin/run-*.sh`, the run scripts need the programs
(`cneta-build`) and R (`cneta-util`) at the same time. Activate both, the
second one stacked on the first:

```bash
conda activate cneta-build
conda activate --stack cneta-util
bin/run-cnets.sh
```

## Running the tests

`cneta-build` also contains pytest. See [Testing](../developer-guide/testing.md#testing)
for the commands.

## Dependencies

What the environment files install, for reference — for example when
building on a system without conda.

Build (`cneta-build`):

- A C++20-capable compiler (GCC >= 10, Clang >= 10, Apple clang >= 12)
- [CMake](https://cmake.org/) >= 3.16
- [Boost](https://www.boost.org/) (`program_options` and `filesystem`
  components)
- [GSL](https://www.gnu.org/software/gsl/) (GNU Scientific Library)
- zlib (used by the vendored `gzstream` for compressed I/O)
- OpenMP — optional, accelerates `cnetml`'s tree search; the build still
  succeeds without it, but `cnetml` then runs single-threaded. See
  [OpenMP](#openmp) above.
- Python 3.9+ with pytest, for the tests. Catch2 is fetched automatically at
  configure time if it is not already installed.

`code/CMakeLists.txt` discovers these via `find_package`; there is nothing
to configure by hand.

`util/` scripts (`cneta-util`):

- R: `optparse`, `copynumber`, `tools`, `tidyverse`, `ggtree`, `treeio`,
  `tidytree`, `ggrepel`, `ggpubr`, `ape`, `phangorn`, `xml2`, `seqinr`,
  `MCMCpack`, `rstudioapi`, plus `Quartet` and `rwty` from CRAN
- Python: `networkx` (used by `newick2elist.py`)

## Containers

Planned, not yet available. In the meantime, the conda environment files
above (see [Environment files](#environment-files)) give a reproducible
installation.
