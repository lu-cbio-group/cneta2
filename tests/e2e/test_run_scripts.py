"""Smoke tests for the run-*.sh driver scripts.

The other e2e tests build their own command lines, which is fast and flexible
but leaves the scripts -- the path the documentation actually tells users to
follow -- untested. These run them as a user would.

Each script honours three environment overrides added for exactly this
purpose, all defaulting to the previous behaviour when unset:

    CNETA_BIN     where the executables live       (default <repo>/bin)
    CNETA_OUT     the working/output directory     (default ./example)
    CNETA_SEED    the random seed                  (default $RANDOM)

run-cnetmcmc.sh also takes CNETA_CONFIG for the configuration file.

These are marked `slow`: the scripts run at their real defaults (seg_max=1000
rather than the 50 the other tests use), so the cnetml one takes ~30 s.
"""

from __future__ import annotations

import pytest

from helpers.outputs import read_copy_numbers, read_tree

pytestmark = pytest.mark.slow

EXPECTED_LEAVES = 4      # Ns=3 in the scripts, plus the normal sample
SEED = "12345"


@pytest.fixture
def script_env(bin_dir, tmp_path):
    """Environment that redirects a run script into a temporary directory."""
    outdir = tmp_path / "example"
    outdir.mkdir()
    return outdir, {
        "CNETA_BIN": str(bin_dir),
        "CNETA_OUT": str(outdir),
        "CNETA_SEED": SEED,
    }


def test_run_cnets(repo_root, script_env, run):
    outdir, env = script_env

    result = run([repo_root / "run-cnets.sh"], cwd=repo_root, env=env)
    assert result.returncode == 0, str(result)

    cn_file = outdir / "sim-data-1-cn.txt.gz"
    tree_file = outdir / "sim-data-1-tree.txt"

    assert cn_file.is_file(), f"run-cnets.sh wrote no copy numbers\n{result}"
    assert tree_file.is_file(), f"run-cnets.sh wrote no tree\n{result}"

    assert read_tree(tree_file).n_leaves == EXPECTED_LEAVES
    assert read_copy_numbers(cn_file).samples == set(range(1, EXPECTED_LEAVES + 1))


def test_run_cnetml_after_cnets(repo_root, script_env, run):
    """run-cnetml.sh consumes what run-cnets.sh produces, so run both."""
    outdir, env = script_env

    run([repo_root / "run-cnets.sh"], cwd=repo_root, env=env).assert_ok()

    result = run([repo_root / "run-cnetml.sh"], cwd=repo_root, env=env,
                 timeout=900)
    assert result.returncode == 0, str(result)

    # The script prints its own verdict; make sure it agrees.
    assert "cnetml main run SUCCEEDED" in result.stdout, str(result)

    # cnetml writes sidecars alongside the tree (.summary.txt, .edge_rates.txt,
    # .nex), and those also end in .txt -- match only the tree itself.
    ml_trees = [path for path in outdir.glob("MaxL-*.txt")
                if path.name.count(".") == 1]
    assert ml_trees, f"run-cnetml.sh wrote no ML tree\n{result}"

    tree = read_tree(ml_trees[0])
    assert tree.n_leaves == EXPECTED_LEAVES
    assert all(edge.length >= 0 for edge in tree.edges)


@pytest.mark.xfail(
    reason="run-cnetmcmc.sh ships with an inconsistent configuration: it "
           "passes is_total=1 with <prefix>-cn.txt.gz (4 columns), but "
           "mcmc.cfg sets model=2, which needs the 5-column "
           "<prefix>-haplotype-cn.txt.gz. cnetmcmc correctly prints 'There "
           "should be 5 columns...' and exits 1, but the script never checks "
           "the status inside its chain loop, so it still prints 'Finish "
           "running cnetmcmc' and exits 0 with no traces written. Two fixes "
           "needed: point input at the haplotype file with is_total=0 (or set "
           "model=1 in mcmc.cfg), and check the exit status per chain the way "
           "run-cnetml.sh already does.",
    strict=False,
)
def test_run_cnetmcmc_after_cnets(repo_root, script_env, data_dir, run):
    outdir, env = script_env
    env = {**env, "CNETA_CONFIG": str(data_dir / "mcmc-ci.cfg")}

    run([repo_root / "run-cnets.sh"], cwd=repo_root, env=env).assert_ok()

    result = run([repo_root / "run-cnetmcmc.sh"], cwd=repo_root, env=env,
                 timeout=900)
    assert result.returncode == 0, str(result)

    traces = list((outdir / "mcmc").glob("mcmc_*.p"))
    assert traces, (
        "run-cnetmcmc.sh produced no parameter trace -- the run failed "
        f"silently\n{result}"
    )


def test_scripts_respect_the_output_override(repo_root, script_env, run):
    """CNETA_OUT must fully redirect output, or a test run would scribble
    into the developer's ./example directory."""
    outdir, env = script_env

    # A developer may well have run the scripts by hand already, so compare
    # before and after rather than requiring ./example to be absent.
    default_dir = repo_root / "example"
    before = set(default_dir.iterdir()) if default_dir.is_dir() else set()

    run([repo_root / "run-cnets.sh"], cwd=repo_root, env=env).assert_ok()

    assert list(outdir.iterdir()), "nothing was written to CNETA_OUT"

    after = set(default_dir.iterdir()) if default_dir.is_dir() else set()
    assert after == before, (
        f"run-cnets.sh wrote to ./example despite CNETA_OUT being set: "
        f"{sorted(path.name for path in after - before)}"
    )
