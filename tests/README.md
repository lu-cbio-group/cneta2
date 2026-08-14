# cneta tests

Three layers, each answering a different question.

| Layer | Directory | Question it answers | Runtime |
|---|---|---|---|
| **Unit** (Catch2, C++) | `unit/` | Do the shared library functions still compute the right thing? | < 1 s |
| **End-to-end** (pytest) | `e2e/` | Do the three programs still run, produce the right files, and give the same numbers? | ~40 s |
| **Build system** (pytest) | `build_system/` | Does the project still configure and compile from scratch? | ~20 s |

The suite is deliberately small. It is a starting point meant to grow, not a
finished safety net — see [What is not covered yet](#what-is-not-covered-yet).

## Before you open a pull request

CI runs on pull requests, not on every pushed commit, so this is the check to
run locally first. From the repository root:

```bash
# 1. Build, with the unit tests enabled
cmake -S code -B code/build -DCNETA_BUILD_TESTING=ON
cmake --build code/build -j 4

# 2. Unit tests
ctest --test-dir code/build --output-on-failure

# 3. End-to-end tests
pip install -r tests/requirements.txt      # once
OMP_NUM_THREADS=1 pytest tests -m "not build"
```

That is exactly what CI does. If all three pass locally, CI should be green.

**One-time setup on macOS** (the whole team is on Apple Silicon):

```bash
brew install cmake boost gsl
```

pytest can come from Homebrew (`brew install pytest`) or from a virtual
environment (`python3 -m venv .venv && . .venv/bin/activate && pip install -r
tests/requirements.txt`). Either works; the suite has no dependencies beyond
pytest itself.

## Running less, or more

```bash
pytest tests                       # everything except the build-system tests
pytest tests -m "not slow"         # ~12 s: skips the run-script smoke tests
pytest tests -m build              # build-system tests only
pytest tests/e2e/test_cnetml.py    # one file
pytest tests -k reproducib -v      # one topic
ctest --test-dir code/build -R model   # one group of unit tests
```

`pytest` skips the whole end-to-end suite with a build hint if the binaries
are missing, rather than reporting a wall of failures.

## How it fits together

**Binaries.** The build writes `cnets`, `cnetml` and `cnetmcmc` to `bin/` at
the repository root. Tests look there; set `CNETA_BIN` to point at a
different build.

**Determinism.** Everything runs with `OMP_NUM_THREADS=1` and a fixed seed.
`cnetml`'s tree search is parallelised with OpenMP, and floating-point sums
over a varying number of threads are not bit-reproducible, so scores would
otherwise drift between runs.

**The regression anchor.** `cnetml --mode 2` scores a supplied tree with no
RNG and no search, so its output is reproducible to the digit.
`tests/data/expected.json` holds the reference values, and
`test_cnetml.py::test_scoring_a_tree_matches_the_reference` compares against
them. **If that test fails after a refactor, the likelihood calculation
changed.** That may be intentional — but it should never happen by accident.
This is the guard the project plan asks for ahead of the `libcneta` work.

**Fixtures.** `tests/data/` holds a tiny committed dataset: 3 tumour regions
plus a normal, 50 sites, simulated with seed 12345. Committed rather than
generated so the regression values do not depend on `cnets` continuing to
simulate the same thing. Regenerate with `tests/regenerate_fixtures.sh`, and
commit any change on its own with an explanation.

**Unit tests link `cneta_core`**, the static library in `code/CMakeLists.txt`
that already holds all the shared science. When that becomes `libcneta`, only
the `target_link_libraries` line in `unit/CMakeLists.txt` needs to change.

**Catch2** is found on the system if installed, otherwise fetched (pinned to
v3.7.1) at configure time. It is only needed with `-DCNETA_BUILD_TESTING=ON`,
which is off by default — a plain `./build.sh local` needs no network access.
Offline, point `-DFETCHCONTENT_SOURCE_DIR_CATCH2=` at a local checkout.

## Known failures recorded in the suite

These run and report, but do not fail the build. Each marks a real problem in
the code, not in the test. Remove the marker when the underlying issue is
fixed — the assertion should then pass as written.

| Where | Issue |
|---|---|
| `unit/test_tree_op.cpp` (`[!mayfail]`) | `order_tree_string_uniq` sorts before dropping the trailing empty field, so it discards the lexicographically largest node instead. Every canonical tree string is missing one node, which could make tree search treat two distinct topologies as the same. |
| `e2e/test_cnets.py` (`xfail`) | `cnets --help` and `cnetml --help` exit 1. Asking for help is not an error, and a non-zero status breaks `cmd --help` in scripts and Makefiles. |
| `e2e/test_run_scripts.py` (`xfail`) | `run-cnetmcmc.sh` passes total copy numbers (`is_total=1`, 4 columns) while `mcmc.cfg` sets `model=2`, which needs the 5-column haplotype file. `cnetmcmc` correctly exits 1, but the script never checks the status inside its chain loop, so it reports success and exits 0 with no traces written. |

## Adding a test

- **A pure function in the C++ library** → `unit/`, next to the closest
  existing case. Add new `.cpp` files to the `add_executable` list in
  `unit/CMakeLists.txt`.
- **Behaviour of a whole program** → `e2e/`. Build the command line with a
  helper from `helpers/cli.py` rather than by hand, overriding only the
  options your test cares about:

  ```python
  run([cnetml, *cli.cnetml_search_args(sim.cn, sim.times, out,
                                       **{"--cn_max": 6})]).assert_ok()
  ```

- **Reading an output file** → add a reader to `helpers/outputs.py` so the
  parsing lives in one place. Standard library only, please.
- **Anything over a few seconds** → mark it `@pytest.mark.slow`.

Use the `tiny_sim` fixture for input; it simulates once per session and is
shared. Only run `cnets` yourself if you need to compare two simulations.

## What is not covered yet

Worth knowing before trusting a green run:

- Only one model configuration is exercised (`model=2`, `cn_type=0`, no
  chromosome-level events and no WGD). The other models and `--cn_type`
  settings are untested.
- `cnetml` modes 1, 3, 4 and 5 have no coverage; only tree search (0) and
  scoring (2) are tested.
- Branch-specific rates (`--bsr_mode`) and the random-local-clock search are
  untested.
- Bootstrapping, ancestral state reconstruction, and the `util/` R scripts
  have no tests.
- The regression anchor covers one dataset. It will catch a change in the
  likelihood; it will not localise it.
