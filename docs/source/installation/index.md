# Installation

`cneta` is written in C++ and built with CMake; a few R and Python scripts
under `util/` support pre/post-processing but are not required to build
the main tools. This page covers what has actually been set up and
verified so far — see the [Developer guide](../developer-guide/index.md)
for the current architecture and what's still planned.

## Dependencies

Required to build `cnets`, `cnetml`, and `cnetmcmc`:

- A C++11-capable compiler (GCC or Clang)
- [CMake](https://cmake.org/) >= 3.10
- [Boost](https://www.boost.org/) (`program_options` and `filesystem`
  components)
- [GSL](https://www.gnu.org/software/gsl/) (GNU Scientific Library)
- zlib (used by the vendored `gzstream` for compressed I/O)
- OpenMP — optional, accelerates `cnetml`'s tree search; the build still
  succeeds without it

`code/CMakeLists.txt` discovers all of these automatically via
`find_package` — there is nothing to configure by hand.

R and Python are only needed for the postprocessing scripts in `util/` —
see [Utility scripts](../user-guide/utility-scripts.md).

## Building from source

```bash
git clone https://github.com/lu-cbio-group/cneta2.git
cd cneta2/code
```

Install dependencies and compile for your platform:

::::{tab-set}
:::{tab-item} macOS / Apple Silicon

```bash
brew install cmake boost gsl
./build.sh local        # or: ./build.sh local <cores>, e.g. ./build.sh local 4
```

Builds are confirmed working on Apple Silicon.
:::

:::{tab-item} Linux (local)

Requires sudo — for shared HPC clusters where you typically don't have
that, see the HPC tab instead. Package names below are for Debian/Ubuntu;
adjust for other distributions.

```bash
sudo apt install cmake libboost-all-dev libgsl-dev zlib1g-dev
./build.sh local        # or: ./build.sh local <cores>, e.g. ./build.sh local 4
```

Not yet verified against a from-scratch Linux desktop install, only
against Eureka2's HPC module environment (HPC tab).
:::

:::{tab-item} HPC

No sudo needed — dependencies are provided as environment modules.
`./build.sh <cluster> <cores>` sources `code/envs/setup_env_<cluster>.sh`
for you automatically; no need to source it yourself first.

```bash
./build.sh Eureka2 <cores>
```

Eureka2 is the currently-supported example; it loads this module set:

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

`build.sh` configures an out-of-source CMake build in `code/build/`. The
three executables (`cnets`, `cnetml`, `cnetmcmc`) land directly in
`code/build/`.

Other `build.sh` usage:

```bash
./build.sh              # no arguments: show help
./build.sh clean         # remove code/build/
```

## Containers

Planned, not yet available.
