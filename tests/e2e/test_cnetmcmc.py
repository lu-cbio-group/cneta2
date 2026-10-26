"""End-to-end tests for cnetmcmc, the Bayesian sampler.

cnetmcmc is the least mature of the three tools, so these tests stay at the
level of "it runs and writes traces of the right shape". They use
tests/data/mcmc-ci.cfg, which is config/cnet_mcmc.cfg with the chain shortened
from 2000 draws to 50 so the suite stays quick.
"""

from __future__ import annotations

import pytest

from helpers import cli
from helpers.outputs import is_balanced_newick, read_mcmc_params, read_mcmc_trees

from conftest import Simulation

# From tests/data/mcmc-ci.cfg: 50 draws, 10 burn-in, sampling every 5th.
N_DRAWS = 50
N_BURNIN = 10
N_GAP = 5

EXPECTED_LEAVES = cli.NS + 1


@pytest.fixture(scope="module")
def mcmc_run(tmp_path_factory, data_dir, tiny_sim: Simulation, cnetmcmc, run):
    """Run one chain and share it across the tests in this module."""
    outdir = tmp_path_factory.mktemp("mcmc")
    trace_param = outdir / "chain.p"
    trace_tree = outdir / "chain.t"

    result = run([
        cnetmcmc,
        *cli.cnetmcmc_args(
            # model=2 in the config means haplotype-specific input.
            cn_file=tiny_sim.haplotype_cn,
            times_file=tiny_sim.times,
            tree_file=tiny_sim.tree,
            config_file=data_dir / "mcmc-ci.cfg",
            trace_param=trace_param,
            trace_tree=trace_tree,
        ),
    ]).assert_ok()

    return result, trace_param, trace_tree


def test_writes_both_trace_files(mcmc_run):
    _, trace_param, trace_tree = mcmc_run

    assert trace_param.is_file(), "no parameter trace written"
    assert trace_tree.is_file(), "no tree trace written"


def test_same_seed_gives_the_same_traces(mcmc_run, tmp_path, data_dir,
                                         tiny_sim: Simulation, cnetmcmc, run):
    """A chain is a random walk, so the seed is the only thing that makes it
    repeatable. Compares the traces byte for byte with the shared run."""
    _, trace_param, trace_tree = mcmc_run
    again_param = tmp_path / "again.p"
    again_tree = tmp_path / "again.t"

    run([
        cnetmcmc,
        *cli.cnetmcmc_args(
            cn_file=tiny_sim.haplotype_cn,
            times_file=tiny_sim.times,
            tree_file=tiny_sim.tree,
            config_file=data_dir / "mcmc-ci.cfg",
            trace_param=again_param,
            trace_tree=again_tree,
        ),
    ]).assert_ok()

    assert again_param.read_bytes() == trace_param.read_bytes(), (
        "same seed, different parameter trace"
    )
    assert again_tree.read_bytes() == trace_tree.read_bytes(), (
        "same seed, different tree trace"
    )


def test_parameter_trace_has_the_expected_number_of_samples(mcmc_run):
    _, trace_param, _ = mcmc_run

    header, samples = read_mcmc_params(trace_param)

    assert header[0] == "state"
    assert "lnl" in header

    # Post-burn-in draws, thinned by n_gap.
    expected = (N_DRAWS - N_BURNIN) // N_GAP
    assert len(samples) == expected, (
        f"expected {expected} post-burn-in samples, got {len(samples)}"
    )


def test_sampled_log_likelihoods_are_finite(mcmc_run):
    _, trace_param, _ = mcmc_run

    header, samples = read_mcmc_params(trace_param)
    lnl_index = header.index("lnl")

    for row in samples:
        value = row[lnl_index]
        assert value == value, "NaN log-likelihood in the trace"
        assert value < 0, f"log-likelihood should be negative, got {value}"


def test_sampled_states_advance_monotonically(mcmc_run):
    _, trace_param, _ = mcmc_run

    _, samples = read_mcmc_params(trace_param)
    states = [row[0] for row in samples]

    assert states == sorted(states), "MCMC state numbers are not increasing"
    assert states[0] >= N_BURNIN, "burn-in samples leaked into the trace"


def test_tree_trace_contains_valid_newick(mcmc_run):
    _, _, trace_tree = mcmc_run

    trees = read_mcmc_trees(trace_tree)
    assert trees, "no trees found in the tree trace"

    for newick in trees:
        assert is_balanced_newick(newick), f"malformed Newick: {newick}"
        # Every leaf appears exactly once, so the commas separating them
        # number one fewer than the leaves.
        assert newick.count(",") == EXPECTED_LEAVES - 1, (
            f"expected {EXPECTED_LEAVES} leaves in {newick}"
        )


def test_shipped_config_loads_and_the_command_line_overrides_it(
    tmp_path, repo_root, tiny_sim: Simulation, cnetmcmc, run
):
    """config/cnet_mcmc.cfg must stay loadable by the current cnetmcmc.

    cnetmcmc rejects unknown keys, so a renamed or removed option would break
    every default run of bin/run-cnetmcmc.sh. The chain is cut short on the
    command line, which takes precedence over the config file -- the sample
    count shows which one won.
    """
    trace_param = tmp_path / "chain.p"
    n_draws, n_burnin, n_gap = 20, 10, 5    # the file says 2000, 1000, 10

    run([
        cnetmcmc,
        *cli.cnetmcmc_args(
            cn_file=tiny_sim.haplotype_cn,
            times_file=tiny_sim.times,
            tree_file=tiny_sim.tree,
            config_file=repo_root / "config" / "cnet_mcmc.cfg",
            trace_param=trace_param,
            trace_tree=tmp_path / "chain.t",
            **{"--n_draws": n_draws, "--n_burnin": n_burnin, "--n_gap": n_gap},
        ),
    ]).assert_ok()

    _, samples = read_mcmc_params(trace_param)
    assert len(samples) == (n_draws - n_burnin) // n_gap, (
        "the command-line chain length did not override the config file"
    )


def test_total_copy_number_input_is_rejected_for_the_haplotype_model(
    tmp_path, data_dir, tiny_sim: Simulation, cnetmcmc, run
):
    """model=2 needs five columns; passing total copy numbers must not be
    silently accepted.

    This is exactly the mismatch that bin/run-cnetmcmc.sh currently ships
    with -- see test_run_scripts.py.
    """
    result = run([
        cnetmcmc,
        *cli.cnetmcmc_args(
            cn_file=tiny_sim.cn,          # total CN: only four columns
            times_file=tiny_sim.times,
            tree_file=tiny_sim.tree,
            config_file=data_dir / "mcmc-ci.cfg",
            trace_param=tmp_path / "chain.p",
            trace_tree=tmp_path / "chain.t",
            **{"--is_total": 1},
        ),
    ])

    combined = result.stdout + result.stderr
    assert "5 columns" in combined, (
        "expected a column-count complaint, got:\n" + str(result)
    )
    assert not (tmp_path / "chain.p").exists(), (
        "cnetmcmc wrote a trace despite rejecting the input"
    )


def test_rejected_input_exits_non_zero(tmp_path, data_dir, tiny_sim: Simulation,
                                       cnetmcmc, run):
    """A fatal input error must be visible to the caller as a non-zero status.

    This is what lets a shell script, CI job, or Nextflow process notice that
    a chain never ran. bin/run-cnetmcmc.sh currently discards it -- see
    test_run_scripts.py.
    """
    result = run([
        cnetmcmc,
        *cli.cnetmcmc_args(
            cn_file=tiny_sim.cn,
            times_file=tiny_sim.times,
            tree_file=tiny_sim.tree,
            config_file=data_dir / "mcmc-ci.cfg",
            trace_param=tmp_path / "chain.p",
            trace_tree=tmp_path / "chain.t",
            **{"--is_total": 1},
        ),
    ])

    assert result.returncode != 0


def test_missing_input_file_is_rejected_with_a_message(tmp_path, data_dir,
                                                       tiny_sim: Simulation,
                                                       cnetmcmc, run):
    """A missing copy-number file used to segfault after printing "0 bins".
    A crash is also non-zero, so require a positive exit code (a signal
    gives a negative one) and a message naming the file."""
    missing = tmp_path / "no-such-file-haplotype-cn.txt.gz"
    result = run([
        cnetmcmc,
        *cli.cnetmcmc_args(
            cn_file=missing,
            times_file=tiny_sim.times,
            tree_file=tiny_sim.tree,
            config_file=data_dir / "mcmc-ci.cfg",
            trace_param=tmp_path / "chain.p",
            trace_tree=tmp_path / "chain.t",
        ),
    ])

    assert result.returncode > 0, (
        f"expected a clean non-zero exit, got {result.returncode}\n{result}"
    )
    assert str(missing) in result.stdout + result.stderr, str(result)


def test_missing_config_file_is_rejected_with_a_message(tmp_path, data_dir,
                                                        tiny_sim: Simulation,
                                                        cnetmcmc, run):
    """A config file that does not exist used to be parsed as empty, so a
    wrong path silently ran the chain with the built-in defaults."""
    missing = tmp_path / "no-such.cfg"
    result = run([
        cnetmcmc,
        *cli.cnetmcmc_args(
            cn_file=tiny_sim.haplotype_cn,
            times_file=tiny_sim.times,
            tree_file=tiny_sim.tree,
            config_file=missing,
            trace_param=tmp_path / "chain.p",
            trace_tree=tmp_path / "chain.t",
        ),
    ])

    assert result.returncode > 0, str(result)
    assert str(missing) in result.stdout + result.stderr, str(result)
    assert not (tmp_path / "chain.p").exists(), (
        "cnetmcmc ran a chain despite the missing config file"
    )
