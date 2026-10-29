# Environment files

Everything that sets up an environment to build or run CNETA, in one place.
Each file says at its top what it is for and how to use it.

| File | Environment | What it is for |
|---|---|---|
| `environment-mac.yml` | `cneta-build` | Build and run the programs on macOS (conda) |
| `environment-linux.yml` | `cneta-build` | Build and run the programs on Linux, HPC with conda, Windows via WSL (conda) |
| `environment-util.yml` | `cneta-util` | The R and Python scripts in `util/` (conda) |
| `setup_env_<cluster>.sh` | modules | HPC without conda: the module loads for one cluster, sourced by `code/build.sh <cluster> <cores>` |

The documentation website has its own file, `docs/environment.yml`, which
stays next to the `docs/requirements.txt` it installs.

Steps per platform: the
[installation guide](https://lu-cbio-group.github.io/cneta2/en/installation/index.html).
