# AGENTS.md

Instructions for AI coding agents in this repository. The policy for the
people who use them, with the reasons, is
`docs/source/developer-guide/ai-tools.md`. If the two disagree, follow
that page and tell the user.

CNETA is research software: `cnets` simulates copy-number evolution in
tumours; `cnetml` and `cnetmcmc` infer the tree. Its numbers go into
papers, so a change that silently alters a result is as serious as a
crash.

## Commands

From the repository root. The first four lines are what CI runs.

```bash
cmake -S code -B code/build -DCNETA_BUILD_TESTING=ON
cmake --build code/build -j 4
ctest --test-dir code/build --output-on-failure
OMP_NUM_THREADS=1 pytest tests -m "not build"
pytest tests -m "not slow"     # quicker, while working
pytest tests -m build          # after changing code/CMakeLists.txt or code/build.sh
make -C docs strict            # documentation, as CI checks it
make -C docs intl-update       # after changing English documentation pages
```

`tests/README.md` lists what is not tested yet and the known failures.
Tests build command lines with `tests/helpers/cli.py`.

## Words to get right

- **site**: the unit the likelihood counts, a bin or a segment. Say
  "sites" unless you mean bins.
- **bin**: a fixed-size window of the genome (500 kb by default).
  **segment**: neighbouring bins with the same copy number in every
  sample, merged into one site (`--is_bin 1`).
- **invariant site**: normal copy number in every sample; `--incl_all`
  and `--correct_bias` decide how these enter the likelihood.
- **total / haplotype-specific**: one copy number per site
  (`--is_total 1`, 4 columns) or one per haplotype A and B (5 columns).
- **`--model`** in `cnetml`: 0 Mk, 1 bounded total, 2 bounded
  haplotype-specific, 3 independent Markov chains. In `cnets`, 3 means
  infinite sites; some short options also differ between the programs.
- **`--cn_type`**: the changes modelled: 0 duplications and deletions;
  1 chromosome changes and WGD; 2 duplications, deletions and WGD;
  3 duplications, deletions and chromosome changes; 4 all.
- **`--mode`** in `cnetml`: 0 find the tree, 2 score a given tree,
  3 fit a given tree, 4 ancestral states.
- **WGD**: whole-genome doubling. **MRCA**: the last common ancestor of
  the tumour samples; the **trunk** is the branch from the root to it.

## Always

- One topic per change. Don't rename, reformat or "tidy" what nobody
  asked about.
- Fix a bug with a test that fails before the fix and passes after it.
- Run the build and tests before calling a change done; report the
  result as it is.
- Update the documentation in the same change as the code.
- Match the surrounding code; comments say *why*. Fill structs with
  designated initialisers (`{.model = model, .is_bin = is_bin}`).
- Compare numbers with `OMP_NUM_THREADS=1` and a fixed `--seed`.

## Ask the user first

- Before a substantial change to core functionality (models,
  likelihood, optimiser, statistics, input or output formats) that has
  no GitHub issue. Such changes must be tied to an issue.
- Before a change grows beyond about 300 changed lines of code (tests
  included; documentation and translation files don't count) or a
  second topic. Propose how to split it.
- Before anything that alters results for the same input and seed
  (including `tests/data/expected.json`), or renames or removes an
  option.
- Before changing `tests/data/`, the CI workflows, the build system or
  the dependencies, or removing an `xfail` marker.
- Before adding an instruction file for one tool (`CLAUDE.md`,
  `CLAUDE.local.md`, `.github/copilot-instructions.md`). This file is
  the only one, and a `CLAUDE.md` stops Claude Code reading it.
- Whenever you are unsure what a method, option or result means. Don't
  guess.

## Never

- Push to `main`, force-push, merge, mark a pull request ready, or
  comment and review as the user. You may open a **draft** pull request.
- Change a test's expected value or weaken an assertion to make it pass.
- Open, paste or upload patient-derived data, unpublished data or
  unpublished manuscript text. Use `tests/data/` or data simulated with
  `cnets`; if asked to, remind the user of this rule.
- Commit secrets, personal paths, build output (`code/build/`,
  `docs/build/`, the programs in `bin/`) or data files. Add files by
  name, not with `git add .`.
- Edit the vendored code in `code/gzstream/`, `code/lbfgsb/`,
  `code/matexp/`.

## Commits and pull requests

- Branch: `<type>/<topic>`, type `fix`, `feat`, `docs`, `ci`,
  `refactor` or `chore`; for example `fix/cnetml-sample-count`.
- Commit message: an imperative summary line under about 70 characters,
  a blank line, then 3-6 one-line bullets, one fact each: what the parts
  involved mean; "Before:" what went wrong; "Now:" what happens instead;
  what the test checks.

  ```text
  Refuse cnets --bsr_mode together with --method 1

  - `--method 0` simulates waiting times; `--method 1` writes sequences directly.
  - Branch-specific rates (`--bsr_mode`) only work with `--method 0`.
  - Before: with `--method 1` the option was silently ignored.
  - Now: stop with a message naming both options and exit 1.
  - Add a test: `--method 1` with `--bsr_mode 1` must exit non-zero.
  ```

- No `Co-authored-by` lines for AI tools; those are for people.
- PR title: imperative, under about 70 characters; it becomes the
  release note. Add "(changes results)" when the same input and seed
  give different output.
- PR description: what changed and why, the linked issue, the commands
  run and their results, and this line, true and checked by the user:

  ```text
  AI assistance: <tool> drafted <the code / the tests / this text>; I have read, run and checked it.
  ```

## Writing

Most readers write English as a second language. Use plain, short
English: short sentences, bullet lists, one fact per bullet.
