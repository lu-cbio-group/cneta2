"""Does the same seed give the same result on every platform?

The other tests check that a seed reproduces a result twice on *this*
machine. These check that it reproduces the *recorded* result, which was
produced on another machine (see "seeded"/"_generated" in expected.json).
CI runs them on x86-64 Linux with GCC and on arm64 macOS with Clang, so a
failure on one runner only means the two platforms disagree. That is the
finding these tests exist to make, not a flaky test: look at the values
before touching the pins.

What is pinned:

- the simulated copy numbers for a fixed seed, as a digest. They are
  integers produced from the random draws alone, so this isolates the
  random-number plumbing (generator, shuffle order, sequence of draws);
- the topology hill climbing returns for a fixed seed on 5 samples, as a
  list of clades, plus its likelihood. The topology is the scientific
  result and is compared exactly; the likelihood is compared loosely
  because the last digits may differ between math libraries.

If a change to the simulator or the search is intended to move these, the
failure message prints the new values ready to paste into expected.json.
Commit them together with the reason.
"""

from __future__ import annotations

import gzip
import hashlib
import json

import pytest

from helpers import cli
from helpers.outputs import clades, parse_search_neg_logl, read_tree

from conftest import Simulation

LOGL_TOLERANCE = 1e-5   # cnetml prints 6 significant figures
SAMPLES = 5             # hill climbing needs at least 5


@pytest.fixture(scope="module")
def expected(data_dir):
    return json.loads((data_dir / "expected.json").read_text())


def _cn_digest(sim: Simulation) -> str:
    """Digest of the decompressed copy numbers (gzip headers carry a time)."""
    return hashlib.sha256(gzip.open(sim.cn, "rb").read()).hexdigest()


def _pins(expected: dict) -> dict:
    """The recorded values; empty if none, so a first run prints them all."""
    return expected.get("seeded", {})


def _explain(key: str, actual) -> str:
    return (
        f"\"{key}\" differs from the recorded value. If this platform "
        "should agree with the one in expected.json, something changed "
        "the random-number plumbing or the search. If the change was "
        "intended, record the new value:\n"
        f"    \"{key}\": {json.dumps(actual)}"
    )


def test_simulated_data_matches_the_recorded_digest(tiny_sim: Simulation,
                                                    expected):
    key = f"cnets_ns{cli.NS}_seed{cli.SEED}_cn_sha256"
    actual = _cn_digest(tiny_sim)
    assert _pins(expected).get(key) == actual, _explain(key, actual)


@pytest.mark.slow
def test_hill_climbing_matches_the_recorded_tree(tmp_path, cnets, cnetml, run,
                                                 expected):
    simdir = tmp_path / "sim"
    simdir.mkdir()
    run([cnets, *cli.cnets_args(simdir, **{"-r": SAMPLES})]).assert_ok()
    sim = Simulation.at(simdir)

    out_file = tmp_path / "MaxL.txt"
    result = run([
        cnetml,
        *cli.cnetml_search_args(sim.cn, sim.times, out_file,
                                **{"--tree_search": 1, "-s": SAMPLES}),
    ]).assert_ok()

    # Check everything before failing, so one run reports every value.
    problems = []

    key = f"cnets_ns{SAMPLES}_seed{cli.SEED}_cn_sha256"
    digest = _cn_digest(sim)
    if _pins(expected).get(key) != digest:
        problems.append(_explain(key, digest))

    key = f"hill_climbing_ns{SAMPLES}_seed{cli.SEED}"
    actual = {
        "neg_logl": parse_search_neg_logl(result.stdout),
        "clades": clades(read_tree(out_file)),
    }
    recorded = _pins(expected).get(key) or {}
    if actual["clades"] != recorded.get("clades"):
        problems.append(_explain(key, actual))
    elif actual["neg_logl"] != pytest.approx(recorded.get("neg_logl"),
                                             rel=LOGL_TOLERANCE):
        problems.append(_explain(key, actual))

    assert not problems, "\n\n".join(problems)
