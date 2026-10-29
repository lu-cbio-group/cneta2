"""Does the optimiser reach the optimum?

Two checks that need no second machine and no reference numbers.

A model with one more shift edge contains the smaller one: a multiplier of
1 on the extra edge gives exactly the smaller model, so the bigger model's
maximised likelihood can never be worse. In an exhaustive shift-edge search
every such nested pair is fitted, so a pair where the bigger model scores
worse means a fit stopped short of its optimum. The researchers' Mac and
Linux runs disagreed because their fits stopped at different points; this
is the same symptom, seen on one machine.

The second check refits the same model from different starting points; at
a true optimum the start does not matter.
"""

from __future__ import annotations

import re

import pytest

from helpers import cli
from helpers.outputs import read_summary

from conftest import Simulation

CANDIDATE = re.compile(
    r"\[RLC\]\[exhaustive\] candidate: K=(\d+) shift_eids=\[([ \d]*)\] nlnl=([-\d.eE+]+)"
)

# The strict stopping settings: the projected-gradient test at 1e-3 instead
# of the default 1e-2, and the relative-reduction test at 1e5 instead of
# 1e7 so that it does not fire first. Fits then agree to the print
# precision; the defaults leave up to about 0.05 on a flat likelihood
# surface, at half the run time.
STRICT = {"--tolerance": 1e-3, "--factr": 1e5}

# cnetml prints six significant figures (0.001 at a -logL of 100), so two
# subsets can only be told apart beyond one unit of the last printed digit;
# the margin covers that plus the optimiser's own tolerance.
TOLERANCE = 5e-3
# What the default settings guarantee on the same data.
TOLERANCE_DEFAULTS = 0.1


def candidates(stdout: str) -> dict[tuple[int, ...], float]:
    """Minimised -logL per shift-edge subset, from the --verbose 1 log."""
    found = {}
    for match in CANDIDATE.finditer(stdout):
        subset = tuple(sorted(int(s) for s in match.group(2).split()))
        found[subset] = float(match.group(3))
    return found


def nested_pair_violations(subsets: dict[tuple[int, ...], float],
                           tolerance: float) -> list[str]:
    """Pairs (S, S + one edge) where the bigger model scored worse."""
    edges = sorted({edge for subset in subsets for edge in subset})
    bad = []
    for subset, nlnl in subsets.items():
        for edge in edges:
            if edge in subset:
                continue
            bigger = tuple(sorted(subset + (edge,)))
            if bigger in subsets and subsets[bigger] > nlnl + tolerance:
                bad.append(f"{list(subset)} -> {list(bigger)}: "
                           f"-logL {nlnl} -> {subsets[bigger]}")
    return bad


def _exhaustive_search(cnetml, run, out_file, cn_file, times_file, tree_file,
                       threads="1", **overrides):
    return run([
        cnetml,
        *cli.cnetml_score_args(
            cn_file=cn_file, times_file=times_file, tree_file=tree_file,
            **{"--mode": 3, "-o": out_file, "--optim": 1, "--estmu": 0,
               "--bsr_mode": 3, "-d": 3, "--rlc_search_method": 0,
               "--verbose": 1, **overrides},
        ),
    ], env={"OMP_NUM_THREADS": threads}).assert_ok()


def _assert_no_violations(subsets, tolerance, expected_subsets, result):
    assert len(subsets) == expected_subsets, (
        f"expected {expected_subsets} subsets in the log\n{result}"
    )
    violations = nested_pair_violations(subsets, tolerance)
    assert not violations, (
        f"{len(violations)} nested pairs where the bigger model scored worse "
        f"by more than {tolerance}:\n" + "\n".join(violations[:20])
    )


def test_exhaustive_search_never_scores_a_bigger_model_worse(tmp_path, data_dir,
                                                            cnetml, run):
    """The tiny fixture with the strict settings: 4 candidate edges, 16
    subsets, 32 nested pairs."""
    result = _exhaustive_search(
        cnetml, run, tmp_path / "fixed.txt",
        data_dir / "tiny-cn.txt.gz", data_dir / "tiny-rel-times.txt",
        data_dir / "tiny-tree.txt", **STRICT,
    )
    _assert_no_violations(candidates(result.stdout), TOLERANCE, 16, result)


def test_exhaustive_search_at_the_defaults_is_close(tmp_path, data_dir, cnetml,
                                                    run):
    """The same search with the default settings, which the run script
    uses. They stop earlier, so a bigger model may score worse by a few
    hundredths, never by whole log units as before."""
    result = _exhaustive_search(
        cnetml, run, tmp_path / "fixed.txt",
        data_dir / "tiny-cn.txt.gz", data_dir / "tiny-rel-times.txt",
        data_dir / "tiny-tree.txt",
    )
    _assert_no_violations(candidates(result.stdout), TOLERANCE_DEFAULTS, 16, result)


@pytest.mark.slow
def test_exhaustive_search_on_simulated_rate_shifts(tmp_path, cnets, cnetml, run):
    """Five samples simulated with rate shifts (cnets --bsr_mode 3): 8
    candidate edges, 256 subsets, 1024 nested pairs, strict settings. The
    same check on the kind of data the researchers run; the subsets are
    fitted in parallel, about a minute on four threads."""
    samples = 5
    simdir = tmp_path / "sim"
    simdir.mkdir()
    run([cnets, *cli.cnets_args(simdir, **{
        "-r": samples, "--method": 0, "--bsr_mode": 3,
        "--bsr_p": 0.5, "--bsr_variance": 1, "--seg_max": 50,
    })]).assert_ok()
    sim = Simulation.at(simdir)

    result = _exhaustive_search(
        cnetml, run, tmp_path / "fixed.txt", sim.cn, sim.times, sim.tree,
        threads="4", **{"-s": samples, **STRICT},
    )
    _assert_no_violations(candidates(result.stdout), TOLERANCE, 256, result)


def _scaled_tree(source, target, factor):
    """Copy a tree file with every branch length multiplied by factor."""
    lines = source.read_text().splitlines()
    out = [lines[0]]
    for line in lines[1:]:
        fields = line.split("\t")
        fields[2] = repr(float(fields[2]) * factor)
        out.append("\t".join(fields))
    target.write_text("\n".join(out) + "\n")
    return target


@pytest.mark.parametrize("constrained", [1, 0])
def test_refits_from_different_starts_agree(tmp_path, data_dir, cnetml, run,
                                            constrained):
    """Mode 3 with the rates estimated (estmu 1), default settings, from
    three starting points. With the time constraint the branch lengths
    follow the sampling times, so the starts differ in the rates; without
    it they differ in the branch lengths as well."""
    tree = data_dir / "tiny-tree.txt"
    if constrained:
        starts = {
            "as-given": {},
            "rates-x10": {"--dup_rate": 0.01, "--del_rate": 0.01},
            "rates-x0.1": {"--dup_rate": 0.0001, "--del_rate": 0.0001},
        }
    else:
        starts = {
            "as-given": {},
            "branches-x1.5": {"--tree_file": _scaled_tree(tree, tmp_path / "longer.txt", 1.5)},
            "branches-x0.6": {"--tree_file": _scaled_tree(tree, tmp_path / "shorter.txt", 0.6)},
        }

    scores = {}
    for name, overrides in starts.items():
        out_file = tmp_path / f"{name}.txt"
        run([
            cnetml,
            *cli.cnetml_score_args(
                cn_file=data_dir / "tiny-cn.txt.gz",
                times_file=data_dir / "tiny-rel-times.txt",
                tree_file=tree,
                **{"--mode": 3, "-o": out_file, "--optim": 1, "--estmu": 1,
                   "--constrained": constrained, **overrides},
            ),
        ]).assert_ok()
        scores[name] = float(read_summary(tmp_path / f"{name}.txt.summary.txt")["raw_logL"])

    spread = max(scores.values()) - min(scores.values())
    assert spread <= TOLERANCE, (
        f"the maximised log-likelihood depends on the starting point: {scores}"
    )
