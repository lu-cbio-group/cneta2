"""Command lines for the three cneta programs.

The parameter values mirror the defaults in the repository's ``run-*.sh``
scripts, shrunk to a size that finishes in seconds. Keeping them here as
dictionaries rather than inline in the tests means a test can override a
single option without restating thirty of them, and means there is one place
to update when the CLI changes.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

# Number of tumour regions. Output has NS + 1 samples, including the normal.
NS = 3

# Fixed seed everywhere, so a failure can always be reproduced by hand.
SEED = 12345


def _flatten(options: Mapping[str, Any]) -> list[str]:
    """Turn ``{"--foo": 1}`` into ``["--foo", "1"]``, dropping ``None``."""
    argv: list[str] = []
    for flag, value in options.items():
        if value is None:
            continue
        argv.extend([flag, str(value)])
    return argv


def cnets_args(outdir: Path, **overrides: Any) -> list[str]:
    """Simulate a small copy-number dataset.

    ``outdir`` needs its trailing separator: cnets concatenates the output
    directory and the file prefix without inserting one.
    """
    options: dict[str, Any] = {
        "-o": f"{outdir}{'' if str(outdir).endswith('/') else '/'}",
        "-p": "",              # empty prefix -> files are named sim-data-1-*
        "-r": NS,
        "-n": 1,               # one simulated patient
        "--mode": 1,           # segments of random size
        "--method": 1,         # simulate sequences directly
        "--fix_nseg": 1,
        "--seg_max": 50,       # run-cnets.sh uses 1000; 50 keeps CI quick
        "--cn_max": 4,
        "--model": 2,          # haplotype-specific bounded model
        "--dup_rate": 0.001,
        "--del_rate": 0.001,
        "--chr_gain": 0,
        "--chr_loss": 0,
        "--wgd": 0,
        "--dup_size": 5,
        "--del_size": 5,
        "-e": 9000000,         # effective population size
        "-b": 1.563e-3,        # exponential growth rate
        "--gtime": 0.002739726,
        "-t": 2,               # time step between sampling times
        "--age": 60,
        "--constrained": 1,
        "--print_relative": 1,
        "--verbose": 0,
        "--seed": SEED,
    }
    options.update(overrides)
    return _flatten(options)


def cnetml_search_args(cn_file: Path, times_file: Path, out_file: Path,
                       **overrides: Any) -> list[str]:
    """Search for the maximum-likelihood tree (``--mode 0``).

    Writes ``out_file`` plus ``out_file.summary.txt``,
    ``out_file.edge_rates.txt`` and two NEXUS files alongside it.
    """
    options: dict[str, Any] = {
        "-c": cn_file,
        "-t": times_file,
        "-o": out_file,
        "-s": NS,
        "--is_total": 1,
        "--is_bin": 0,
        "--incl_all": 1,
        "--m_max": 1,
        "-d": 2,               # model
        "--cn_max": 4,
        "--cn_type": 0,        # segment-level changes only
        "--tree_search": 2,    # exhaustive: tiny and deterministic at NS=3
        "--init_tree": 0,
        "-p": 100,             # population size for the genetic algorithm
        "-g": 20,
        "-e": 5,
        "-r": 1e-2,            # convergence tolerance
        "--epop": 90000,
        "--beta": 1.563e-3,
        "--gtime": 0.002739726,
        "--optim": 1,          # L-BFGS-B
        "--constrained": 1,
        "--estmu": 1,
        "--correct_bias": 1,
        "-x": 0,
        "--dup_rate": 0.001,
        "--del_rate": 0.001,
        "--chr_gain_rate": 0,
        "--chr_loss_rate": 0,
        "--wgd_rate": 0,
        "--speed_nni": 1,
        "--mode": 0,
        "--verbose": 0,
        "--seed": SEED,
    }
    options.update(overrides)
    return _flatten(options)


def cnetml_score_args(cn_file: Path, times_file: Path, tree_file: Path,
                      **overrides: Any) -> list[str]:
    """Score one supplied tree (``--mode 2``).

    No RNG and no tree search, so this is the one invocation whose output is
    reproducible to the digit. It prints the result rather than writing a
    file: "The log likelihood of the input tree is <x>".
    """
    options: dict[str, Any] = {
        "-c": cn_file,
        "-t": times_file,
        "--tree_file": tree_file,
        "-s": NS,
        "--is_total": 1,
        "--is_bin": 0,
        "--incl_all": 1,
        "--m_max": 1,
        "-d": 2,
        "--cn_max": 4,
        "--cn_type": 0,
        "--constrained": 1,
        "--correct_bias": 1,
        "-x": 0,
        "--dup_rate": 0.001,
        "--del_rate": 0.001,
        "--chr_gain_rate": 0,
        "--chr_loss_rate": 0,
        "--wgd_rate": 0,
        "--mode": 2,
        "--verbose": 0,
        "--seed": SEED,
    }
    options.update(overrides)
    return _flatten(options)


def cnetmcmc_args(cn_file: Path, times_file: Path, tree_file: Path,
                  config_file: Path, trace_param: Path, trace_tree: Path,
                  **overrides: Any) -> list[str]:
    """Run the MCMC sampler.

    ``--is_total 0`` with a haplotype-specific input file, because the
    configuration uses ``model=2``. Passing total copy numbers with that model
    makes cnetmcmc stop with a column-count error -- see test_run_scripts.py.
    """
    options: dict[str, Any] = {
        "-s": NS,
        "--is_total": 0,
        "-c": cn_file,
        "-t": times_file,
        "--rtree": "",
        "--trace_param_file": trace_param,
        "--trace_tree_file": trace_tree,
        "--config_file": config_file,
        "--init_tree": 0,
        "--file_itree": tree_file,
        "--seed": SEED,
    }
    options.update(overrides)
    return _flatten(options)
