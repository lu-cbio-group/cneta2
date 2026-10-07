# Configuration

## Command-line options

Each tool's options are documented on its own User guide page, grouped
by purpose rather than as a flat alphabetical list:

- [`cnets` options](../user-guide/cnets.md#options) — tree generation,
  mutation model, event rates, branch-specific rates.
- [`cnetml` options](../user-guide/cnetml.md#options) — time
  constraints, evolutionary model, independent Markov chain dimensions,
  tree search, mutation types.
- [`cnetmcmc` options](../user-guide/cnetmcmc.md#options) — reference
  tree/init tree, and everything read from `config/cnet_mcmc.cfg`.

There is no single global options reference; each binary parses its own
set independently.

## Run scripts and config files

The driver scripts live in `bin/`, next to the compiled programs, and their
configuration in `config/`:

```text
bin/
├── run-cnets.sh
├── run-cnetml.sh
└── run-cnetmcmc.sh
config/
├── common.conf      settings shared by all three scripts (shell syntax)
└── cnet_mcmc.cfg    parameters for cnetmcmc (program-options syntax)
```

### Run scripts

The current pattern is: **edit the script, don't type flags by hand.**
Each `bin/run-*.sh` sets its tool's options as shell variables in a block
at the top of the file (input paths, model parameters, etc.), then builds
the command line from those variables. To run with different settings,
copy the relevant script (or edit it in place) and change the variable
values — see the "Examples" section on each tool's
[User guide](../user-guide/index.md) page.

The scripts find the programs, `config/` and `util/` relative to the
repository, so they can be run from any directory. Output goes to
`CNETA_OUT` (default `./example`), which is relative to the directory you
run them from — usually the top of the repository:

```bash
bin/run-cnets.sh
```

### Shared settings: `config/common.conf`

Sourced by all three scripts before anything else. It holds settings that
belong to a whole pipeline run rather than to one tool:

| Variable | Meaning | Default |
|---|---|---|
| `CNETA_OUT` | Working directory. `cnets` writes its simulated data here; `cnetml` and `cnetmcmc` read their input from here and write their results alongside. | `./example` |
| `CNETA_SEED` | Random seed. Empty means a new random seed for every run; the scripts write the seed they used at the top of their log. Set a number to reproduce a run. | empty |

It is read by the shell, so use shell syntax: no spaces around `=`. Write
each setting as `VAR="${VAR:-value}"`, so that a value set in the
environment still takes precedence over the file.

### Environment overrides

Every setting can be overridden for a single run from the environment,
e.g. `CNETA_SEED=42 bin/run-cnets.sh`. The order of precedence is:
environment, then `config/common.conf`, then the script's own default.

| Variable | Used by | Default |
|---|---|---|
| `CNETA_OUT` | all scripts | `./example` |
| `CNETA_SEED` | all scripts | a new random seed each run |
| `CNETA_BIN` | all scripts | `bin/` in the repository |
| `CNETA_CONFIG_DIR` | all scripts | `config/` in the repository |
| `CNETA_CONFIG` | `run-cnetmcmc.sh` | `$CNETA_CONFIG_DIR/cnet_mcmc.cfg` |

### Tool config files: `config/*.cfg`

`cnetmcmc` reads most of its parameters from `config/cnet_mcmc.cfg`,
passed with `--config_file` — see
[MCMC configuration and trace files](../file-formats/index.md#mcmc-configuration-and-trace-files)
for its format. `bin/run-cnetmcmc.sh` stops with an error if the file
does not exist, because `cnetmcmc` itself would silently fall back to its
built-in defaults.

A value given on the command line overrides the same key in the file. A
parameter the script passes explicitly — for example `--seed` or
`--is_total` — therefore cannot be changed from the config file; change it
in the script.

`cnets` and `cnetml` don't read a config file yet; they take all options
as command-line flags.

:::{warning}
`.conf` and `.cfg` files look alike but are not interchangeable. A `.cfg`
file is read by the program and accepts `dup_rate = 0.001`; a `.conf` file
is read by the shell, where the same line is an error. Keep the extension
that matches the reader.
:::

## Planned

More config files are planned: one per tool (for `cnets` and `cnetml`) and
task-specific configurations. This page will describe them as they land.
