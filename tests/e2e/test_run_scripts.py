"""Smoke tests for the bin/run-*.sh driver scripts.

The other e2e tests build their own command lines, which is fast and flexible
but leaves the scripts -- the path the documentation actually tells users to
follow -- untested. These run them as a user would.

Each script honours these environment overrides, all defaulting to the
normal behaviour when unset:

    CNETA_BIN         where the executables live    (default <repo>/bin)
    CNETA_CONFIG_DIR  where the config files live   (default <repo>/config)
    CNETA_OUT         the working/output directory  (default ./example)
    CNETA_SEED        the random seed               (default $RANDOM)

run-cnetmcmc.sh also takes CNETA_CONFIG for its configuration file.

The scripts source config/common.conf, which may set CNETA_OUT and CNETA_SEED;
a value in the environment wins over the file.

These are marked `slow`: the scripts run at their real defaults (seg_max=1000
rather than the 50 the other tests use), so the cnetml one takes ~30 s.
"""

from __future__ import annotations

import pytest

from helpers.outputs import read_copy_numbers, read_tree

pytestmark = pytest.mark.slow

EXPECTED_LEAVES = 4      # Ns=3 in the scripts, plus the normal sample
SEED = "12345"


@pytest.fixture(scope="session")
def scripts_dir(repo_root):
    """Where the run scripts live -- always the repository's own bin/, even
    when CNETA_BIN points the scripts at a different build."""
    return repo_root / "bin"


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


def test_run_cnets(repo_root, scripts_dir, script_env, run):
    outdir, env = script_env

    result = run([scripts_dir / "run-cnets.sh"], cwd=repo_root, env=env)
    assert result.returncode == 0, str(result)

    cn_file = outdir / "sim-data-1-cn.txt.gz"
    tree_file = outdir / "sim-data-1-tree.txt"

    assert cn_file.is_file(), f"run-cnets.sh wrote no copy numbers\n{result}"
    assert tree_file.is_file(), f"run-cnets.sh wrote no tree\n{result}"

    assert read_tree(tree_file).n_leaves == EXPECTED_LEAVES
    assert read_copy_numbers(cn_file).samples == set(range(1, EXPECTED_LEAVES + 1))


def test_run_cnetml_after_cnets(repo_root, scripts_dir, script_env, run):
    """run-cnetml.sh consumes what run-cnets.sh produces, so run both."""
    outdir, env = script_env

    run([scripts_dir / "run-cnets.sh"], cwd=repo_root, env=env).assert_ok()

    result = run([scripts_dir / "run-cnetml.sh"], cwd=repo_root, env=env,
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
           "config/cnet_mcmc.cfg sets model=2, which needs the 5-column "
           "<prefix>-haplotype-cn.txt.gz. cnetmcmc correctly prints 'There "
           "should be 5 columns...' and exits 1, but the script never checks "
           "the status inside its chain loop, so it still prints 'Finish "
           "running cnetmcmc' and exits 0 with no traces written. Two fixes "
           "needed: point input at the haplotype file with is_total=0 (or set "
           "model=1 in the config), and check the exit status per chain the way "
           "run-cnetml.sh already does.",
    strict=False,
)
def test_run_cnetmcmc_after_cnets(repo_root, scripts_dir, script_env, data_dir,
                                  run):
    outdir, env = script_env
    env = {**env, "CNETA_CONFIG": str(data_dir / "mcmc-ci.cfg")}

    run([scripts_dir / "run-cnets.sh"], cwd=repo_root, env=env).assert_ok()

    result = run([scripts_dir / "run-cnetmcmc.sh"], cwd=repo_root, env=env,
                 timeout=900)
    assert result.returncode == 0, str(result)

    traces = list((outdir / "mcmc").glob("mcmc_*.p"))
    assert traces, (
        "run-cnetmcmc.sh produced no parameter trace -- the run failed "
        f"silently\n{result}"
    )


def test_scripts_respect_the_output_override(repo_root, scripts_dir, script_env,
                                              run):
    """CNETA_OUT must fully redirect output, or a test run would scribble
    into the developer's ./example directory."""
    outdir, env = script_env

    # A developer may well have run the scripts by hand already, so compare
    # before and after rather than requiring ./example to be absent.
    default_dir = repo_root / "example"
    before = set(default_dir.iterdir()) if default_dir.is_dir() else set()

    run([scripts_dir / "run-cnets.sh"], cwd=repo_root, env=env).assert_ok()

    assert list(outdir.iterdir()), "nothing was written to CNETA_OUT"

    after = set(default_dir.iterdir()) if default_dir.is_dir() else set()
    assert after == before, (
        f"run-cnets.sh wrote to ./example despite CNETA_OUT being set: "
        f"{sorted(path.name for path in after - before)}"
    )


# --------------------------------------------------------------------------
# Locations and configuration
# --------------------------------------------------------------------------

def _seed_in_log(outdir):
    """The seed the script wrote at the top of cnets' log."""
    (log,) = outdir.glob("std_cnets_*")
    return log.read_text().splitlines()[0]


@pytest.fixture
def config_dir(tmp_path):
    """A config directory whose common.conf points somewhere recognisable."""
    configured_out = tmp_path / "from-common-conf"
    directory = tmp_path / "config"
    directory.mkdir()
    (directory / "common.conf").write_text(
        f'CNETA_OUT="${{CNETA_OUT:-{configured_out}}}"\n'
        'CNETA_SEED="${CNETA_SEED:-777}"\n'
    )
    return directory, configured_out


def test_common_conf_supplies_the_defaults(scripts_dir, bin_dir, config_dir,
                                           tmp_path, run):
    directory, configured_out = config_dir
    # Empty counts as unset, so this also neutralises any value exported in
    # the developer's own shell.
    env = {"CNETA_BIN": str(bin_dir), "CNETA_CONFIG_DIR": str(directory),
           "CNETA_OUT": "", "CNETA_SEED": ""}

    run([scripts_dir / "run-cnets.sh"], cwd=tmp_path, env=env).assert_ok()

    assert (configured_out / "sim-data-1-cn.txt.gz").is_file(), (
        "run-cnets.sh ignored CNETA_OUT from common.conf"
    )
    assert _seed_in_log(configured_out) == "seed 777"


def test_the_environment_overrides_common_conf(scripts_dir, script_env,
                                               config_dir, tmp_path, run):
    outdir, env = script_env
    directory, configured_out = config_dir
    env = {**env, "CNETA_CONFIG_DIR": str(directory)}

    run([scripts_dir / "run-cnets.sh"], cwd=tmp_path, env=env).assert_ok()

    assert (outdir / "sim-data-1-cn.txt.gz").is_file()
    assert _seed_in_log(outdir) == f"seed {SEED}"
    assert not configured_out.exists(), (
        "common.conf overrode CNETA_OUT from the environment"
    )


def test_scripts_run_from_any_directory(scripts_dir, bin_dir, tmp_path, run):
    """The scripts find the programs and config/ relative to the repository;
    a relative CNETA_OUT stays relative to where they are run from."""
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    env = {"CNETA_BIN": str(bin_dir), "CNETA_OUT": "out", "CNETA_SEED": SEED}

    run([scripts_dir / "run-cnets.sh"], cwd=elsewhere, env=env).assert_ok()

    assert (elsewhere / "out" / "sim-data-1-cn.txt.gz").is_file()


def test_run_cnetmcmc_stops_when_its_config_is_missing(scripts_dir, script_env,
                                                       tmp_path, run):
    """cnetmcmc would silently fall back to its built-in defaults, so the
    script must refuse to start instead."""
    _, env = script_env
    missing = tmp_path / "no-such.cfg"
    env = {**env, "CNETA_CONFIG": str(missing)}

    result = run([scripts_dir / "run-cnetmcmc.sh"], cwd=tmp_path, env=env)

    assert result.returncode != 0, str(result)
    assert str(missing) in result.stderr, str(result)
    assert "Start running cnetmcmc" not in result.stdout
