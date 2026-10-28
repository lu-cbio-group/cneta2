# Contributing

**Workflow**: branch from `main` → make your change → open a pull
request into `main` → get it reviewed if the change warrants it → merge.
For a docs change, `make -C docs strict` locally is the same check CI
runs, so a local pass means CI will pass too — see
[Working on the documentation](documentation.md).

**Branches**: one branch per pull request, made from an up-to-date `main`
and deleted after the merge. Name it `<type>/<topic>`, lower case with
hyphens, two to four words:

| Type | For | Example |
|---|---|---|
| `fix/` | a bug, a crash, a wrong message | `fix/cnetml-sample-count` |
| `feat/` | a new option, model or output | `feat/version-file` |
| `docs/` | documentation only | `docs/conda-environments` |
| `ci/` | workflows and tests | `ci/plot-validation` |
| `refactor/` | code moved or reshaped, results unchanged | `refactor/search-into-library` |
| `chore/` | housekeeping: `.gitignore`, dependencies | `chore/drop-r-gitignore-rule` |

Exploratory work that is not yet meant for `main` uses your initials
instead of a type, e.g. `rw/rlc-aicc`.

**Pull request titles**: in the imperative and under about 70 characters,
e.g. `Reject invalid copy-number files with clear error messages`. The
titles become the release notes (see [Versions and releases](releases.md)).

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
