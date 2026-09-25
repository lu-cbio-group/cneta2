"""End-to-end tests for cnetml, the maximum-likelihood inference tool.

Two very different kinds of test live here.

`--mode 2` scores a supplied tree: no RNG, no tree search, so its output is
reproducible to the digit and can be compared against a committed reference.
That is the regression anchor for the libcneta refactor.

`--mode 0` searches for a tree. It is stochastic and parallel, so the tests
check structure and ranking instead of exact numbers.
"""

from __future__ import annotations

import json

import pytest

from helpers import cli
from helpers.outputs import (
    parse_scored_logl,
    parse_search_neg_logl,
    read_summary,
    read_tree,
)

from conftest import Simulation

EXPECTED_LEAVES = cli.NS + 1

# cnetml prints 6 significant figures, which sets the floor on how tightly
# these can be compared.
LOGL_TOLERANCE = 1e-5


@pytest.fixture(scope="module")
def expected(data_dir):
    return json.loads((data_dir / "expected.json").read_text())


# --------------------------------------------------------------------------
# mode 2: scoring a given tree
# --------------------------------------------------------------------------

@pytest.mark.parametrize("tree_name", ["tiny-tree.txt", "tiny-tree-alt.txt"])
def test_scoring_a_tree_matches_the_reference(tree_name, data_dir, cnetml, run,
                                              expected):
    """The regression test: a fixed tree and dataset must keep scoring the same.

    If this fails after a refactor, the likelihood calculation changed. That
    might be intentional -- but it must never happen by accident.
    """
    result = run([
        cnetml,
        *cli.cnetml_score_args(
            cn_file=data_dir / "tiny-cn.txt.gz",
            times_file=data_dir / "tiny-rel-times.txt",
            tree_file=data_dir / tree_name,
        ),
    ]).assert_ok()

    logl = parse_scored_logl(result.stdout)
    reference = expected["mode2_logl"][tree_name]

    assert logl == pytest.approx(reference, rel=LOGL_TOLERANCE), (
        f"log-likelihood for {tree_name} moved from {reference} to {logl}.\n"
        "If this was intentional, update tests/data/expected.json and say why "
        "in the commit message."
    )


def test_scoring_is_reproducible(data_dir, cnetml, run):
    """Same inputs, same answer -- twice in the same session."""
    scores = []
    for _ in range(2):
        result = run([
            cnetml,
            *cli.cnetml_score_args(
                cn_file=data_dir / "tiny-cn.txt.gz",
                times_file=data_dir / "tiny-rel-times.txt",
                tree_file=data_dir / "tiny-tree.txt",
            ),
        ]).assert_ok()
        scores.append(parse_scored_logl(result.stdout))

    assert scores[0] == scores[1]


def test_true_tree_beats_a_wrong_topology(data_dir, cnetml, run):
    """A cheap sanity check on the science, not just the plumbing.

    tiny-tree.txt is the topology the data was simulated under;
    tiny-tree-alt.txt is a different topology on the same leaves. The true
    one must score higher, or the likelihood is not measuring what it claims.
    """
    scores = {}
    for tree_name in ("tiny-tree.txt", "tiny-tree-alt.txt"):
        result = run([
            cnetml,
            *cli.cnetml_score_args(
                cn_file=data_dir / "tiny-cn.txt.gz",
                times_file=data_dir / "tiny-rel-times.txt",
                tree_file=data_dir / tree_name,
            ),
        ]).assert_ok()
        scores[tree_name] = parse_scored_logl(result.stdout)

    assert scores["tiny-tree.txt"] > scores["tiny-tree-alt.txt"], (
        f"the true topology scored worse than a wrong one: {scores}"
    )


# --------------------------------------------------------------------------
# mode 0: searching for a tree
# --------------------------------------------------------------------------

def test_tree_search_writes_a_usable_tree(tmp_path, tiny_sim: Simulation,
                                          cnetml, run):
    out_file = tmp_path / "MaxL.txt"

    result = run([
        cnetml,
        *cli.cnetml_search_args(tiny_sim.cn, tiny_sim.times, out_file),
    ]).assert_ok()

    assert out_file.is_file(), str(result)

    tree = read_tree(out_file)
    assert tree.n_leaves == EXPECTED_LEAVES
    assert tree.n_edges == 2 * EXPECTED_LEAVES - 2
    assert all(edge.length >= 0 for edge in tree.edges), (
        "inferred tree has a negative branch length"
    )

    # cnetml reports the negated log-likelihood, so this is positive.
    neg_logl = parse_search_neg_logl(result.stdout)
    assert neg_logl > 0, f"expected a positive -logL, got {neg_logl}"

    # stdout and the summary sidecar must describe the same run.
    summary = read_summary(out_file.with_name(out_file.name + ".summary.txt"))
    assert float(summary["raw_logL"]) == pytest.approx(-neg_logl, rel=LOGL_TOLERANCE)


def test_tree_search_writes_a_consistent_summary(tmp_path, tiny_sim: Simulation,
                                                 cnetml, run):
    """The .summary.txt sidecar is the machine-readable record of a run."""
    out_file = tmp_path / "MaxL.txt"

    run([
        cnetml,
        *cli.cnetml_search_args(tiny_sim.cn, tiny_sim.times, out_file),
    ]).assert_ok()

    summary_file = out_file.with_name(out_file.name + ".summary.txt")
    assert summary_file.is_file()

    summary = read_summary(summary_file)

    assert summary["mode"] == "0"
    assert summary["model"] == "2"
    assert int(summary["n_leaves"]) == EXPECTED_LEAVES
    assert int(summary["n_edges"]) == 2 * EXPECTED_LEAVES - 2

    # The summary's score must agree with the tree file it describes.
    assert float(summary["raw_logL"]) < 0

    # Estimated rates must be usable numbers, not NaN escaping the optimiser.
    for key in ("dup_rate_reference", "del_rate_reference"):
        rate = float(summary[key])
        assert rate == rate, f"{key} is NaN"
        assert rate >= 0, f"{key} is negative: {rate}"


def test_tree_search_writes_edge_rates(tmp_path, tiny_sim: Simulation,
                                       cnetml, run):
    out_file = tmp_path / "MaxL.txt"

    run([
        cnetml,
        *cli.cnetml_search_args(tiny_sim.cn, tiny_sim.times, out_file),
    ]).assert_ok()

    edge_rates = out_file.with_name(out_file.name + ".edge_rates.txt")
    assert edge_rates.is_file()

    lines = [line for line in edge_rates.read_text().splitlines() if line.strip()]
    header, rows = lines[0].split("\t"), lines[1:]

    assert header[0] == "eid"
    # One row per edge.
    assert len(rows) == 2 * EXPECTED_LEAVES - 2


def test_seg_file_is_written_not_read(tmp_path, tiny_sim: Simulation,
                                      cnetml, run):
    """--seg_file names an *output*, despite reading like an input.

    run-cnetml.sh points it at `<prefix>-segs.txt`, a file cnets never
    creates, which looks like a bug until you notice cnetml writes it. This
    test pins the direction down so nobody else has to work it out.
    """
    out_file = tmp_path / "MaxL.txt"
    seg_file = tmp_path / "does-not-exist-yet-segs.txt"

    assert not seg_file.exists()

    run([
        cnetml,
        *cli.cnetml_search_args(tiny_sim.cn, tiny_sim.times, out_file,
                                **{"--seg_file": seg_file}),
    ]).assert_ok()

    assert seg_file.is_file(), "cnetml did not write the --seg_file output"
    assert seg_file.stat().st_size > 0


# --------------------------------------------------------------------------
# Error handling
# --------------------------------------------------------------------------

def test_missing_input_file_is_rejected(tmp_path, data_dir, cnetml, run):
    result = run([
        cnetml,
        *cli.cnetml_score_args(
            cn_file=tmp_path / "no-such-file-cn.txt.gz",
            times_file=data_dir / "tiny-rel-times.txt",
            tree_file=data_dir / "tiny-tree.txt",
        ),
    ])

    assert result.returncode != 0, (
        "cnetml accepted a missing copy-number file\n" + str(result)
    )
