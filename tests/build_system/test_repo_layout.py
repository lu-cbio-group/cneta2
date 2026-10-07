"""Checks on the repository layout that are easy to break without noticing.

bin/ holds both the tracked run scripts and the untracked programs the build
writes there. Ignoring bin/ as a whole would keep existing scripts (git keeps
tracking them) but silently leave out any new one, so check the patterns
directly. Paths need not exist for `git check-ignore`.

Fast and needs no build, so -- unlike test_build.py -- this runs in CI.
"""

from __future__ import annotations

import shutil

import pytest

PROGRAMS = ("cnets", "cnetml", "cnetmcmc")


@pytest.fixture(scope="module")
def ignored(repo_root, run):
    git = shutil.which("git")
    if git is None:
        pytest.skip("git is not installed")
    if run([git, "rev-parse", "--is-inside-work-tree"], cwd=repo_root).returncode:
        pytest.skip("not a git checkout")

    def _ignored(path: str) -> bool:
        return run([git, "check-ignore", "-q", path], cwd=repo_root).returncode == 0

    return _ignored


@pytest.mark.parametrize("program", PROGRAMS)
def test_built_programs_are_ignored(ignored, program):
    assert ignored(f"bin/{program}")


@pytest.mark.parametrize("path", [
    "bin/run-new-tool.sh",       # a script added later, not yet tracked
    "config/new-tool.cfg",
    "config/common.conf",
])
def test_new_scripts_and_configs_are_not_ignored(ignored, path):
    assert not ignored(path), f".gitignore would silently leave out {path}"
