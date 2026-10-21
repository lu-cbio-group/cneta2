"""One version number, from the VERSION file, everywhere it appears.

The version lives in VERSION at the repository root (MAJOR.MINOR.PATCH).
The build passes it to the programs, which print it with --version and
record it in their outputs next to the seed, so that any result can be
traced to the code that produced it. See
docs/source/developer-guide/releases.md.
"""

from __future__ import annotations

import re

import pytest

from helpers import cli

from conftest import Simulation


def test_version_file_holds_major_minor_patch(cneta_version):
    assert re.fullmatch(r"\d+\.\d+\.\d+", cneta_version), (
        f"VERSION must hold MAJOR.MINOR.PATCH, e.g. 2.1.0; found {cneta_version!r}"
    )


@pytest.mark.parametrize("program", ["cnets", "cnetml", "cnetmcmc"])
def test_programs_report_the_version_file(program, request, run,
                                          cneta_version):
    """--version prints the number in VERSION (plus the git commit, when the
    build knows it). Only the text is checked: --version currently exits 1,
    a known failure covered in test_cnets.py."""
    result = run([request.getfixturevalue(program), "--version"])

    assert f"[version {cneta_version}" in result.stdout, str(result)


def test_every_run_logs_the_version(cnets, run, tmp_path, cneta_version):
    """Every run prints the version once, just before the seed, so a captured
    log is enough to tell which code and which seed produced a result."""
    outdir = tmp_path / "out"
    outdir.mkdir()
    result = run([cnets, *cli.cnets_args(outdir)]).assert_ok()

    lines = result.stdout.splitlines()
    version_lines = [i for i, line in enumerate(lines)
                     if line.startswith("version:\t")]
    seed_lines = [i for i, line in enumerate(lines) if line.startswith("seed:\t")]
    assert len(version_lines) == 1, f"expected one version line\n{result}"
    assert lines[version_lines[0]].startswith(f"version:\t{cneta_version}")
    assert version_lines[0] < seed_lines[0], "the version should precede the seed"


def test_the_simulation_records_the_version(tiny_sim: Simulation,
                                            cneta_version):
    """The record travels with the data: <prefix>-info.txt starts with the
    version, followed by the seeds."""
    info = tiny_sim.info.read_text()

    assert info.startswith(f"VERSION: {cneta_version}"), info[:200]
    assert "\nSEED: " in info[:400], info[:400]
