"""End-to-end tests for cnetml, the maximum-likelihood inference tool.

Two very different kinds of test live here.

`--mode 2` scores a supplied tree: no RNG, no tree search, so its output is
reproducible to the digit and can be compared against a committed reference.
That is the regression anchor for the libcneta refactor.

`--mode 0` searches for a tree. It is stochastic and parallel, so the tests
check structure and ranking instead of exact numbers.
"""

from __future__ import annotations

import gzip
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


def test_cn_type_is_required(data_dir, cnetml, run):
    """`--cn_type` decides which event types the model considers, and there is
    no right default. It used to have none and no check either: without it the
    variable was never set, and mode 3 stopped with "Unknown CN model type!"
    after reading the uninitialised value. It is now a required option, so
    leaving it out gives a clean error that names it."""
    result = run([
        cnetml,
        *cli.cnetml_score_args(
            cn_file=data_dir / "tiny-cn.txt.gz",
            times_file=data_dir / "tiny-rel-times.txt",
            tree_file=data_dir / "tiny-tree.txt",
            **{"--cn_type": None},
        ),
    ])

    assert_clean_rejection(result, mentions="cn_type")


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
                                                 cnetml, run, cneta_version):
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
    # The seed and the version travel with the result, whatever launched
    # the program.
    assert summary["seed"] == str(cli.SEED)
    assert summary["version"].split()[0] == cneta_version
    assert int(summary["n_leaves"]) == EXPECTED_LEAVES
    assert int(summary["n_edges"]) == 2 * EXPECTED_LEAVES - 2

    # The summary's score must agree with the tree file it describes.
    assert float(summary["raw_logL"]) < 0

    # Estimated rates must be usable numbers, not NaN escaping the optimiser.
    for key in ("dup_rate_reference", "del_rate_reference"):
        rate = float(summary[key])
        assert rate == rate, f"{key} is NaN"
        assert rate >= 0, f"{key} is negative: {rate}"


def test_fixed_tree_with_branch_specific_rates_runs(tmp_path, data_dir, cnetml,
                                                     run):
    """Mode 3 with bsr_mode 3 (random local clock) on the tiny fixture: the
    path the researchers use for rate shifts, and the only test that runs
    the shift-edge search. Checks that the run completes and that the
    summary records the settings. About 0.04 s on the fixture."""
    out_file = tmp_path / "fixed.txt"
    run([
        cnetml,
        *cli.cnetml_score_args(
            cn_file=data_dir / "tiny-cn.txt.gz",
            times_file=data_dir / "tiny-rel-times.txt",
            tree_file=data_dir / "tiny-tree.txt",
            **{"--mode": 3, "-o": out_file, "--optim": 1, "--estmu": 0,
               "--bsr_mode": 3},
        ),
    ]).assert_ok()

    summary = read_summary(tmp_path / "fixed.txt.summary.txt")
    assert summary["mode"] == "3"
    assert summary["bsr_mode"] == "3"
    assert summary["rlc_criterion"] == "AIC"
    assert float(summary["raw_logL"]) < 0
    assert (tmp_path / "fixed.txt.edge_rates.txt").is_file(), (
        "mode 3 with bsr_mode 3 wrote no edge-rate table"
    )


@pytest.mark.parametrize("incl_all", [1, 0])
def test_shift_edge_search_states_its_parameter_count(incl_all, tmp_path, data_dir,
                                                      cnetml, run):
    """The shift-edge search prints what its information criteria rest on:
    the number of sites behind BIC's sample size, the shift parameters per
    edge, and the shared parameters that only AICc includes. On the tiny
    fixture (3 regions, age-constrained tree, segment-level rates) that is
    3 node-age parameters plus 2 reference rates, so 5 shared ones.

    The sample size is the number of sites in the likelihood, whatever the
    input layout: all 50 sites with incl_all 1, or the 9 variable sites
    plus the 41 invariant ones that the bias correction adds back with
    incl_all 0. Both must report n=50."""
    out_file = tmp_path / "fixed.txt"
    result = run([
        cnetml,
        *cli.cnetml_score_args(
            cn_file=data_dir / "tiny-cn.txt.gz",
            times_file=data_dir / "tiny-rel-times.txt",
            tree_file=data_dir / "tiny-tree.txt",
            **{"--mode": 3, "-o": out_file, "--optim": 1, "--estmu": 0,
               "--bsr_mode": 3, "-d": 3, "--incl_all": incl_all},
        ),
    ])
    result.assert_ok()

    assert "information criteria: n=50 sites in the likelihood" in result.stdout
    assert "AIC/BIC count K x 2 shift parameters" in result.stdout
    assert "AICc also counts the 5 shared ones" in result.stdout


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


def _search_outputs(out_dir, sim: Simulation, cnetml, run, threads: str,
                    **overrides) -> tuple[bytes, bytes]:
    """Run one seeded tree search with the given thread count; return the
    tree file and its summary sidecar as bytes."""
    out_dir.mkdir()
    out_file = out_dir / "MaxL.txt"
    run([
        cnetml,
        *cli.cnetml_search_args(sim.cn, sim.times, out_file, **overrides),
    ], env={"OMP_NUM_THREADS": threads}).assert_ok()

    summary = out_file.with_name(out_file.name + ".summary.txt")
    return out_file.read_bytes(), summary.read_bytes()


def test_same_seed_gives_the_same_tree_whatever_the_thread_count(
        tmp_path, tiny_sim: Simulation, cnetml, run):
    """The search is parallelised with OpenMP. Each thread optimises its own
    copy of a tree and no random number is drawn inside a parallel region,
    so the thread count must not change the answer -- the seed alone does.
    """
    outputs = [_search_outputs(tmp_path / f"threads-{threads}", tiny_sim,
                               cnetml, run, threads)
               for threads in ("1", "4")]

    assert outputs[0] == outputs[1], "the thread count changed the result"


@pytest.mark.slow
def test_hill_climbing_is_reproducible_across_thread_counts(tmp_path, cnets,
                                                            cnetml, run):
    """Hill climbing draws random numbers between its parallel sections
    (random start trees, then random perturbations of the best ones), so it
    is the search most likely to pick up a thread-order dependence. It needs
    at least 5 samples, hence its own simulation."""
    samples = 5
    simdir = tmp_path / "sim"
    simdir.mkdir()
    run([cnets, *cli.cnets_args(simdir, **{"-r": samples})]).assert_ok()
    sim = Simulation.at(simdir)

    outputs = [_search_outputs(tmp_path / f"threads-{threads}", sim, cnetml,
                               run, threads,
                               **{"--tree_search": 1, "-s": samples})
               for threads in ("1", "4")]

    assert outputs[0] == outputs[1], "the thread count changed the result"


# --------------------------------------------------------------------------
# Error handling
# --------------------------------------------------------------------------

def _write_gz(path, text: str):
    with gzip.open(path, "wt") as handle:
        handle.write(text)
    return path


def assert_clean_rejection(result, mentions: str):
    """The program refused the input with a message, rather than crashing.

    A crash also exits non-zero, so `returncode != 0` alone cannot tell the
    two apart: killed by a signal, subprocess reports a negative code
    (-11 for a segfault). This test once passed while cnetml segfaulted.
    """
    assert result.returncode > 0, (
        "expected a clean non-zero exit, got "
        + ("a crash (signal %d)" % -result.returncode if result.returncode < 0
           else "success")
        + "\n" + str(result)
    )
    assert mentions in result.stdout + result.stderr, (
        f"the error message does not mention {mentions!r}\n{result}"
    )


@pytest.mark.parametrize("problem", ["missing", "empty", "short-line"])
def test_unusable_input_file_is_rejected_with_a_message(tmp_path, data_dir,
                                                        cnetml, run, problem):
    """A missing file, an empty one, and a line with too few columns all used
    to reach code that indexes into the data (two of them segfaulted)."""
    cn_file = tmp_path / f"{problem}-cn.txt.gz"
    if problem == "empty":
        _write_gz(cn_file, "")
    elif problem == "short-line":
        _write_gz(cn_file, "1\t1\n")

    result = run([
        cnetml,
        *cli.cnetml_score_args(
            cn_file=cn_file,
            times_file=data_dir / "tiny-rel-times.txt",
            tree_file=data_dir / "tiny-tree.txt",
        ),
    ])

    assert_clean_rejection(result, mentions=str(cn_file))


def test_sample_count_larger_than_the_file_is_rejected(data_dir, cnetml, run):
    """`-s` larger than the number of samples in the file (here: the default
    `-s 5` on the 4-sample fixture) left samples without data, and cnetml
    segfaulted after "Reading data and group regions by chromosome". It
    must stop with a message that names the file and the missing sample.

    No sample-times file and no age constraint: the times reader already
    checks its own entry count against `-s`, which would mask this."""
    cn_file = data_dir / "tiny-cn.txt.gz"
    result = run([
        cnetml,
        *cli.cnetml_score_args(
            cn_file=cn_file,
            times_file=None,
            tree_file=data_dir / "tiny-tree.txt",
            **{"-s": cli.NS + 2, "--constrained": 0},
        ),
    ])

    assert_clean_rejection(result, mentions=str(cn_file))
    assert f"sample {cli.NS + 2} has no rows" in result.stdout + result.stderr, str(result)
