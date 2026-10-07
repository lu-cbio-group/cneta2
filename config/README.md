# config/

Configuration read by the run scripts in `bin/` and by the programs.

| File | Read by | Syntax |
|---|---|---|
| `common.conf` | all three `bin/run-*.sh` scripts | shell |
| `cnet_mcmc.cfg` | `cnetmcmc`, via `--config_file` (passed by `bin/run-cnetmcmc.sh`) | program options |

The two syntaxes look alike but are not interchangeable:

- **`.conf` — shell.** Sourced by the scripts, so no spaces around `=`.
  Write settings as `VAR="${VAR:-value}"` so the environment can override them.
- **`.cfg` — program options.** One `key=value` per line, where `key` is the
  program's long option name without the leading `--`. `#` starts a comment.
  An unknown key is an error.

A value given on the program's command line overrides the same key in its
`.cfg` file.

The scripts look for this directory relative to the repository. To use
configuration kept elsewhere, set `CNETA_CONFIG_DIR`; to use a different MCMC
config file, set `CNETA_CONFIG`.

See the [Configuration](../docs/source/configuration/index.md) page of the
documentation for details.
