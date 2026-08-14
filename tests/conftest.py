"""Shared fixtures for the cneta end-to-end tests.

Every test here drives a compiled binary, so the fixtures are mostly about
locating those binaries, running them safely, and simulating one small
dataset that the rest of the suite can share.
"""

from __future__ import annotations

import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pytest

# Make `helpers` importable from the test modules without turning the whole
# tree into a package.
sys.path.insert(0, str(Path(__file__).parent))

from helpers import cli  # noqa: E402

PROGRAMS = ("cnets", "cnetml", "cnetmcmc")

BUILD_HINT = (
    "cneta binaries not found in {bin_dir}. Build them first:\n"
    "    cd code && ./build.sh local 4\n"
    "or point CNETA_BIN at an existing build directory."
)


# --------------------------------------------------------------------------
# Locations
# --------------------------------------------------------------------------

@pytest.fixture(scope="session")
def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


@pytest.fixture(scope="session")
def data_dir() -> Path:
    """Committed fixtures, small enough to live in git."""
    return Path(__file__).resolve().parent / "data"


@pytest.fixture(scope="session")
def bin_dir(repo_root: Path) -> Path:
    """Where the three executables live.

    The CMake build writes them to <repo>/bin; CNETA_BIN overrides that, which
    is how CI points the tests at an out-of-tree build.
    """
    override = os.environ.get("CNETA_BIN")
    return Path(override).resolve() if override else repo_root / "bin"


@pytest.fixture(scope="session")
def binaries(bin_dir: Path) -> dict[str, Path]:
    """Skip the whole suite politely when the project has not been built."""
    missing = [name for name in PROGRAMS if not (bin_dir / name).is_file()]
    if missing:
        pytest.skip(BUILD_HINT.format(bin_dir=bin_dir), allow_module_level=True)
    return {name: bin_dir / name for name in PROGRAMS}


@pytest.fixture(scope="session")
def cnets(binaries: dict[str, Path]) -> Path:
    return binaries["cnets"]


@pytest.fixture(scope="session")
def cnetml(binaries: dict[str, Path]) -> Path:
    return binaries["cnetml"]


@pytest.fixture(scope="session")
def cnetmcmc(binaries: dict[str, Path]) -> Path:
    return binaries["cnetmcmc"]


# --------------------------------------------------------------------------
# Running things
# --------------------------------------------------------------------------

@dataclass
class Result:
    argv: list[str]
    returncode: int
    stdout: str
    stderr: str

    def __str__(self) -> str:
        return (
            f"command: {' '.join(self.argv)}\n"
            f"exit code: {self.returncode}\n"
            f"--- stdout (tail) ---\n{self.stdout[-3000:]}\n"
            f"--- stderr (tail) ---\n{self.stderr[-3000:]}"
        )

    def assert_ok(self) -> "Result":
        assert self.returncode == 0, f"command failed\n{self}"
        return self


@pytest.fixture(scope="session")
def run():
    """Run a command with a reproducible environment.

    OMP_NUM_THREADS=1 matters: cnetml parallelises its tree search with
    OpenMP, and floating-point reductions over a variable number of threads
    are not bit-reproducible. Pinning to one thread makes scores comparable
    between runs and between machines.
    """

    def _run(argv, cwd: Path | None = None, env: dict[str, str] | None = None,
             timeout: int = 600) -> Result:
        argv = [str(a) for a in argv]

        full_env = os.environ.copy()
        full_env["OMP_NUM_THREADS"] = "1"
        if env:
            full_env.update(env)

        completed = subprocess.run(
            argv,
            cwd=str(cwd) if cwd else None,
            env=full_env,
            capture_output=True,
            text=True,
            timeout=timeout,
        )

        return Result(argv, completed.returncode, completed.stdout, completed.stderr)

    return _run


# --------------------------------------------------------------------------
# A shared simulated dataset
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Simulation:
    """Paths to the files one cnets run produces."""

    directory: Path
    cn: Path
    haplotype_cn: Path
    tree: Path
    times: Path
    info: Path

    @classmethod
    def at(cls, directory: Path, prefix: str = "sim-data-1") -> "Simulation":
        return cls(
            directory=directory,
            cn=directory / f"{prefix}-cn.txt.gz",
            haplotype_cn=directory / f"{prefix}-haplotype-cn.txt.gz",
            tree=directory / f"{prefix}-tree.txt",
            times=directory / f"{prefix}-rel-times.txt",
            info=directory / f"{prefix}-info.txt",
        )


@pytest.fixture(scope="session")
def tiny_sim(tmp_path_factory, cnets: Path, run) -> Simulation:
    """One small simulated dataset, shared by every test that needs input.

    Session-scoped because simulating is the cheapest step and repeating it
    per test would triple the suite's runtime for no extra coverage. Tests
    that need to compare two simulations run cnets themselves.
    """
    directory = tmp_path_factory.mktemp("tiny_sim")

    run([cnets, *cli.cnets_args(directory)]).assert_ok()

    sim = Simulation.at(directory)
    assert sim.cn.is_file(), f"cnets did not write {sim.cn}"

    return sim
