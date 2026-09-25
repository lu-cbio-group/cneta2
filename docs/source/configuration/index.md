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
  tree/init tree, and everything read from `mcmc.cfg`.

There is no single global options reference; each binary parses its own
set independently.

## Config files and run scripts

The repository ships `run-cnets.sh`, `run-cnetml.sh`, and
`run-cnetmcmc.sh` together with `mcmc.cfg`.

The current pattern is: **edit the script, don't type flags by hand.**
Each `run-*.sh` sets its tool's options as shell variables in a block at
the top of the file (input/output paths, model parameters, etc.), then
builds the command line from those variables. To run with different
settings, copy the relevant script (or edit it in place) and change the
variable values — see the "Examples" section on each tool's
[User guide](../user-guide/index.md) page.

`cnetmcmc` additionally reads most of its parameters from a separate
config file, `mcmc.cfg`, passed with `--config_file` — see
[MCMC configuration and trace files](../file-formats/index.md#mcmc-configuration-and-trace-files)
for its format. `cnets` and `cnetml` take all options as command-line
flags; they don't read a config file.

## Planned structured configuration

TODO: schema, validation, presets.

:::{note}
No roadmap for this exists in `README.md` or elsewhere in the repo yet —
left as a placeholder rather than invented. Fill in once there's an
actual plan to describe.
:::
