"""Tests for the build system itself.

CI builds the project anyway, so these are excluded from the default CI run
(`-m "not build"`). They exist so that a developer changing CMakeLists.txt or
build.sh can check the contract locally, and so the build script's own
argument handling is covered.

The configure-and-compile test builds into a temporary directory and never
touches code/build/ or bin/.
"""

from __future__ import annotations

import shutil

import pytest

pytestmark = [pytest.mark.build, pytest.mark.slow]

PROGRAMS = ("cnets", "cnetml", "cnetmcmc")


@pytest.fixture(scope="module")
def cmake() -> str:
    path = shutil.which("cmake")
    if path is None:
        pytest.skip("cmake is not installed")
    return path


# --------------------------------------------------------------------------
# build.sh argument handling
# --------------------------------------------------------------------------

def test_build_script_shows_help_without_arguments(repo_root, run):
    result = run(["./build.sh"], cwd=repo_root / "code")

    assert result.returncode == 0, str(result)
    assert "cneta Build Script" in result.stdout
    assert "Usage:" in result.stdout


def test_build_script_rejects_an_unknown_target(repo_root, run):
    # Anything that is not 'local' or 'clean' is treated as a cluster name and
    # must have a matching envs/setup_env_<name>.sh.
    result = run(["./build.sh", "NotACluster"], cwd=repo_root / "code")

    assert result.returncode != 0, str(result)
    assert "Cannot find environment script" in result.stdout


def test_build_script_rejects_a_non_numeric_core_count(repo_root, run):
    result = run(["./build.sh", "local", "many"], cwd=repo_root / "code")

    assert result.returncode != 0, str(result)
    assert "must be a positive integer" in result.stdout


# --------------------------------------------------------------------------
# A clean configure-and-compile
# --------------------------------------------------------------------------

def test_configures_and_builds_from_scratch(repo_root, tmp_path, cmake, run):
    build_dir = tmp_path / "build"
    bin_dir = tmp_path / "bin"

    configure = run([
        cmake,
        "-S", repo_root / "code",
        "-B", build_dir,
        f"-DCMAKE_RUNTIME_OUTPUT_DIRECTORY={bin_dir}",
    ], timeout=900)
    assert configure.returncode == 0, str(configure)

    compile_result = run([cmake, "--build", build_dir, "-j", "4"], timeout=1800)
    assert compile_result.returncode == 0, str(compile_result)

    for name in PROGRAMS:
        program = bin_dir / name
        assert program.is_file(), f"{name} was not built"
        assert program.stat().st_mode & 0o111, f"{name} is not executable"


def test_the_output_directory_holds_only_the_user_facing_tools(
    repo_root, tmp_path, cmake, run
):
    """Test binaries and vendored libraries must not land in bin/.

    `ls bin/` should keep meaning "the programs a user runs", even when the
    tests are enabled.
    """
    build_dir = tmp_path / "build"
    bin_dir = tmp_path / "bin"

    run([
        cmake,
        "-S", repo_root / "code",
        "-B", build_dir,
        "-DCNETA_BUILD_TESTING=ON",
        f"-DCMAKE_RUNTIME_OUTPUT_DIRECTORY={bin_dir}",
    ], timeout=900).assert_ok()

    run([cmake, "--build", build_dir, "-j", "4"], timeout=1800).assert_ok()

    assert {entry.name for entry in bin_dir.iterdir()} == set(PROGRAMS)


def test_unit_tests_are_off_by_default(repo_root, tmp_path, cmake, run):
    """A plain `./build.sh local` must not need network access for Catch2."""
    build_dir = tmp_path / "build"

    run([
        cmake,
        "-S", repo_root / "code",
        "-B", build_dir,
        f"-DCMAKE_RUNTIME_OUTPUT_DIRECTORY={tmp_path / 'bin'}",
    ], timeout=900).assert_ok()

    # No test target is generated, so ctest finds nothing to run.
    listing = run([cmake, "--build", build_dir, "--target", "help"], timeout=300)
    assert "cneta_unit_tests" not in listing.stdout
