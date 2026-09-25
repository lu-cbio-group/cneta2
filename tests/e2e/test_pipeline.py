"""The full chain: simulate, infer, sample.

This is the project's acceptance criterion in test form -- "executes
successfully end-to-end on synthetic data". It is also the test that will be
replaced by, or wrapped in, the Nextflow pipeline later on.

Unlike the per-tool tests, this one runs everything in a single directory so
that the file names each stage produces really are the ones the next stage
looks for.
"""

from __future__ import annotations

from helpers import cli
from helpers.outputs import (
    is_balanced_newick,
    read_copy_numbers,
    read_mcmc_trees,
    read_summary,
    read_tree,
)

from conftest import Simulation

EXPECTED_LEAVES = cli.NS + 1


def test_simulate_infer_sample(tmp_path, data_dir, cnets, cnetml, cnetmcmc, run):
    workdir = tmp_path / "pipeline"
    workdir.mkdir()

    # ---- 1. Simulate -----------------------------------------------------
    run([cnets, *cli.cnets_args(workdir)]).assert_ok()

    sim = Simulation.at(workdir)
    assert sim.cn.is_file()
    assert sim.tree.is_file()
    assert sim.times.is_file()

    truth = read_tree(sim.tree)
    assert truth.n_leaves == EXPECTED_LEAVES

    observed = read_copy_numbers(sim.cn)
    assert observed.samples == set(range(1, EXPECTED_LEAVES + 1))

    # ---- 2. Infer a tree from the simulated copy numbers -----------------
    ml_tree = workdir / "MaxL.txt"
    run([
        cnetml,
        *cli.cnetml_search_args(sim.cn, sim.times, ml_tree),
    ]).assert_ok()

    assert ml_tree.is_file()

    inferred = read_tree(ml_tree)
    # The inferred tree spans the same taxa as the tree the data came from.
    assert inferred.n_leaves == truth.n_leaves
    assert inferred.n_edges == truth.n_edges

    summary = read_summary(ml_tree.with_name(ml_tree.name + ".summary.txt"))
    assert int(summary["n_leaves"]) == EXPECTED_LEAVES
    assert float(summary["raw_logL"]) < 0

    # ---- 3. Sample the posterior from the same simulated data ------------
    trace_param = workdir / "chain.p"
    trace_tree = workdir / "chain.t"
    run([
        cnetmcmc,
        *cli.cnetmcmc_args(
            cn_file=sim.haplotype_cn,
            times_file=sim.times,
            tree_file=sim.tree,
            config_file=data_dir / "mcmc-ci.cfg",
            trace_param=trace_param,
            trace_tree=trace_tree,
        ),
    ]).assert_ok()

    trees = read_mcmc_trees(trace_tree)
    assert trees, "MCMC produced no sampled trees"
    assert all(is_balanced_newick(newick) for newick in trees)


def test_inferred_tree_can_be_rescored(tmp_path, cnets, cnetml, run):
    """cnetml must be able to read back the tree it just wrote.

    Round-tripping through the tree file format is the interface every
    downstream tool depends on, and the one most likely to drift during the
    libcneta refactor.
    """
    workdir = tmp_path / "roundtrip"
    workdir.mkdir()

    run([cnets, *cli.cnets_args(workdir)]).assert_ok()
    sim = Simulation.at(workdir)

    ml_tree = workdir / "MaxL.txt"
    run([
        cnetml,
        *cli.cnetml_search_args(sim.cn, sim.times, ml_tree),
    ]).assert_ok()

    # Feed the inferred tree straight back in as an input.
    rescored = run([
        cnetml,
        *cli.cnetml_score_args(sim.cn, sim.times, ml_tree),
    ]).assert_ok()

    from helpers.outputs import parse_scored_logl

    assert parse_scored_logl(rescored.stdout) < 0
