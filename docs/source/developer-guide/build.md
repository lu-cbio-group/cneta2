# Build system

## Repository layout

```text
cneta/
├── bin/           what you run
│   ├── run-*.sh   example driver scripts (in git)
│   └── cnets, cnetml, cnetmcmc   built executables (written by the build, not in git)
├── config/        configuration: common.conf (shared by the scripts), one .cfg per tool
├── code/          C++ sources for cnets, cnetml, cnetmcmc
│   ├── gzstream/  vendored: gzip streams
│   ├── lbfgsb/    vendored: L-BFGS-B optimiser
│   └── matexp/    vendored: matrix exponential
├── tests/         unit, end-to-end and build-system tests
├── util/          R and Python helper scripts
├── assets/        logo artwork
└── docs/          this documentation site
```

`bin/` mixes tracked scripts with untracked build output. `.gitignore`
excludes the three programs by name rather than ignoring `bin/` as a whole,
so a new script added there is picked up by git as usual;
`tests/build_system/test_repo_layout.py` checks this.

## Building

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
