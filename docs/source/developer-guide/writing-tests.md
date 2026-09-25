# Writing a unit test

Unit tests call the C++ functions directly and check what they return. They
are the cheapest tests to write and the fastest to run — the whole suite
finishes in about a second — so they are the right place for anything that can
be checked without running a whole program.

This page is the walkthrough for adding one. See
[Testing](index.md#testing) for what the three test layers are and how to run
them all.

## What can be tested this way

The unit tests link `cneta_core`, the static library built from the shared
sources (`model.cpp`, `tree_op.cpp`, `parse_cn.cpp`, `stats.cpp`,
`likelihood.cpp` and the rest). Any function declared in one of their headers
can be called from a test.

Functions defined in `cnets.cpp`, `cnetml.cpp` or `cnetmcmc.cpp` cannot: those
files are compiled into the three programs, not into the library, so a test
that calls them will fail to link. Behaviour that only exists inside one of
the programs belongs in `tests/e2e/` instead.

## Adding a test file

There are two steps, and the second is easy to forget.

1. **Create `tests/unit/test_<area>.cpp`** — or, if one of the existing files
   already covers the area, add your `TEST_CASE` there instead. The four files
   are `test_model.cpp`, `test_tree_op.cpp`, `test_parse_cn.cpp` and
   `test_stats.cpp`.

2. **List the new file** in the `add_executable` block of
   `tests/unit/CMakeLists.txt`:

   ```cmake
   add_executable(cneta_unit_tests
       test_stats.cpp
       test_model.cpp
       test_parse_cn.cpp
       test_tree_op.cpp
       test_mine.cpp)      # <- your new file
   ```

:::{warning}
A file that is not listed there is never compiled. There is no error message —
your tests simply do not exist. If a test you just wrote does not appear in the
output, check this first.
:::

## Building and running

From the repository root:

```bash
cmake -S code -B code/build -DCNETA_BUILD_TESTING=ON
cmake --build code/build -j 4
ctest --test-dir code/build --output-on-failure
```

The first command is only needed once, and again whenever you edit
`CMakeLists.txt`. After that, rebuild and re-run with the last two commands.

Note that configuring alone builds nothing — if you skip the `cmake --build`
step, `ctest` runs the previous executable and your new test will not be there.

To run only your own tests while working on them:

```bash
# by test name (a regular expression)
ctest --test-dir code/build -R "largest node time" --output-on-failure

# by tag, calling the test executable directly
code/build/tests/cneta_unit_tests "[tree_op]"

# see everything that exists, with its tags
code/build/tests/cneta_unit_tests --list-tests
```

Running the executable directly also shows anything the code prints to
`stdout`, which `ctest` hides unless a test fails.

## The shape of a test

```cpp
TEST_CASE("get_tree_height returns the largest node time", "[tree_op]"){
    std::vector<double> node_times{0.0, 3.5, 1.25, 7.0, 2.0};

    REQUIRE_THAT(get_tree_height(node_times),
                 Catch::Matchers::WithinAbs(7.0, 1e-12));
}
```

- The **name** is what `ctest` prints, so write it as a statement of what must
  be true rather than as a label. "returns the largest node time" tells a
  future reader what broke; "test1" does not.
- The **tag** in square brackets groups tests so they can be run together. Use
  the one already used by the file you are working in: `[model]`, `[tree_op]`,
  `[parse_cn]` or `[stats]`.
- **`REQUIRE`** stops the test case at the first failure. Use `CHECK` instead
  if you want the remaining assertions in that test to run anyway.

## Comparing numbers

`REQUIRE(a == b)` is right for integers, states and sizes. For `double` values
it is almost always wrong: arithmetic that is mathematically equivalent can
differ in the last bits, so an exact comparison fails for no useful reason.

Compare doubles with a tolerance instead, using one of the two matchers from
`<catch2/matchers/catch_matchers_floating_point.hpp>`:

| Matcher | Use it when | Example |
|---|---|---|
| `WithinAbs(expected, tol)` | the expected value is zero or close to it | a rate-matrix row must sum to `0.0` |
| `WithinRel(expected, tol)` | everything else — the tolerance scales with the size of the number | a likelihood, a branch length, a rate |

```cpp
REQUIRE_THAT(row_sum(q.data(), n, i), Catch::Matchers::WithinAbs(0.0, 1e-12));
REQUIRE_THAT(rate, Catch::Matchers::WithinRel(0.03, 1e-12));
```

`WithinRel` cannot be used against an expected value of zero: a relative
tolerance around zero is zero, so the comparison becomes exact again.

## Matrices are stored column-major

Every rate and transition matrix in the codebase is a flat `double*` in
column-major order, so:

```cpp
m[i + j * n]   // row i, column j, in an n x n matrix
```

Writing `m[j + i * n]` by mistake reads the transpose. That is worth being
careful about, because many of these matrices are close to symmetric, so a
transposed test can still pass and hide the error. `test_model.cpp` has a small
`row_sum` helper that shows the indexing in use.

## Checking the same property for several inputs

When a property should hold for many values, loop over them and use `CAPTURE`
so that a failure reports which values caused it:

```cpp
for(const auto& rates : kRatePairs){
    const double dup_rate = rates.first;
    const double del_rate = rates.second;
    CAPTURE(dup_rate, del_rate);

    // ... assertions using those rates
}
```

Without `CAPTURE`, a failure tells you the assertion broke but not which
iteration broke it.

Choose the values deliberately. Inputs that are equal to each other, or zero,
can make a wrong result look identical to a right one — if two rates are both
`0.01`, then code that uses them in the wrong places still produces the correct
matrix and the test cannot detect the mistake. Varying the values, and the
ratios between them, is what gives the test its power.

## Using a fixture file

`tests/data/` holds a small committed dataset. Tests reach it through
`CNETA_TEST_DATA_DIR`, a path that `tests/unit/CMakeLists.txt` defines at
compile time so that tests work regardless of the directory `ctest` runs them
from:

```cpp
std::string data_path(const std::string& name){
    return std::string(CNETA_TEST_DATA_DIR) + "/" + name;
}
```

`test_tree_op.cpp` and `test_parse_cn.cpp` both use this. Prefer an existing
fixture over adding a new one; if you do add a file, commit it on its own with
an explanation of what generated it.

## Recording a bug you cannot fix yet

If a test documents a genuine bug that is not going to be fixed in the same
change, tag it `[!mayfail]`. It still runs and still reports, but it does not
fail the build:

```cpp
TEST_CASE("order_tree_string_uniq is idempotent", "[tree_op][!mayfail]"){
```

Write the assertion the way it *should* pass, add a row to the known-failures
table in `tests/README.md` explaining the underlying problem, and remove the
tag when the bug is fixed. This is much better than deleting the test or
weakening it until it passes.

## Worked examples in the repository

Rather than copying a template, read the test that already does what you need.
All of these are compiled on every build, so they cannot fall out of date:

| What you want to do | Look at |
|---|---|
| The simplest possible test | `test_model.cpp` — *"compute_haplotype_change_dim follows the parity rule"* |
| Compare doubles; group variations with `SECTION` | `test_model.cpp` — *"get_rate_matrix_bounded is a valid generator"* |
| Run the same check over several inputs with `CAPTURE` | `test_model.cpp` — *"hand-filled site-level haplotype rate matrix matches the generated one"* |
| Read a fixture file | `test_tree_op.cpp` — *"read_tree_info loads the fixture tree"* |
| Check two implementations agree | `test_model.cpp` — the two *"hand-filled ... matches the generated one"* cases |
| Record a known bug | `test_tree_op.cpp` — *"order_tree_string_uniq is idempotent"* |

## Unit test or end-to-end test?

A function that can be called directly and returns something checkable belongs
here. Anything that needs a program to run — command-line handling, output
files, the `run-*.sh` scripts, reproducibility of a whole simulation — belongs
in `tests/e2e/` as a pytest test.

`tests/README.md` covers that side, along with what the suite does not cover
yet and the current list of known failures.
