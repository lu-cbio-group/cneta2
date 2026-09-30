# Testing and CI

## Testing

The suite lives in `tests/` and has three layers:

`tests/unit/`
: Catch2 unit tests for the shared library. They link `cneta_core` directly
  and cover the state encoding, the rate and transition matrices, and the
  tree helpers. Built only with `-DCNETA_BUILD_TESTING=ON`.

`tests/e2e/`
: pytest tests that run the compiled programs on a small simulated dataset,
  checking output files, error handling, and reproducibility. Includes smoke
  tests that execute the `run-*.sh` scripts as a user would.

`tests/build_system/`
: pytest tests that configure and compile the project from scratch and check
  `build.sh`'s argument handling.

To run everything:

```bash
cmake -S code -B code/build -DCNETA_BUILD_TESTING=ON
cmake --build code/build -j 4
ctest --test-dir code/build --output-on-failure
pip install -r tests/requirements.txt
OMP_NUM_THREADS=1 pytest tests -m "not build"
```

`OMP_NUM_THREADS=1` matters: `cnetml` parallelises its tree search with
OpenMP, and floating-point sums over a varying number of threads are not
bit-reproducible.

The regression anchor is `cnetml --mode 2`, which scores a supplied tree with
no RNG and no search and so is reproducible to the digit. Reference values
live in `tests/data/expected.json`. A change there means the likelihood
calculation changed — which is exactly what needs watching during the
`libcneta` refactor.

[Writing a unit test](writing-tests.md) is the walkthrough for adding one of
your own.

`tests/README.md` has the full detail, including what is *not* covered yet and
a list of known failures the suite records but does not fail on.

## Continuous integration

Workflows live in `.github/workflows/`. Currently:

`docs.yml`
: builds the documentation on every pull request into `main`, and
  publishes to GitHub Pages when a pull request is merged.

`ci.yml`
: builds the project and runs the unit and end-to-end tests on
  `ubuntu-latest` and `macos-latest`.

Both run on pull requests into `main` and on pushes to `main` — not on every
pushed commit, so work-in-progress on a branch with no open pull request
costs nothing. Both can also be triggered by hand from the Actions tab
(`workflow_dispatch`).

TODO: container workflow.
