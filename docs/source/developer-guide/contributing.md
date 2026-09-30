# Contributing

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
