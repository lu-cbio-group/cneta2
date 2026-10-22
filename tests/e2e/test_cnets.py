"""End-to-end tests for cnets, the simulator.

cnets is the head of the pipeline: everything downstream consumes what it
writes, so these tests care about the shape and reproducibility of its
output rather than about the biology.
"""

from __future__ import annotations

import gzip
import hashlib

import pytest

from helpers import cli
from helpers.outputs import read_copy_numbers, read_tree

from conftest import Simulation

EXPECTED_SITES = 50          # --seg_max in cli.cnets_args
EXPECTED_SAMPLES = cli.NS + 1  # tumour regions plus the normal


def test_writes_the_expected_files(tiny_sim: Simulation):
    for path in (tiny_sim.cn, tiny_sim.haplotype_cn, tiny_sim.tree,
                 tiny_sim.times, tiny_sim.info):
        assert path.is_file(), f"cnets did not write {path.name}"
        assert path.stat().st_size > 0, f"{path.name} is empty"


def test_copy_number_file_is_well_formed(tiny_sim: Simulation):
    cn = read_copy_numbers(tiny_sim.cn)

    # Sample IDs run 1..Ns+1 with no gaps, as the file-format docs require.
    assert cn.samples == set(range(1, EXPECTED_SAMPLES + 1))

    # Every sample is profiled at the same sites.
    assert set(cn.sites_per_sample.values()) == {EXPECTED_SITES}

    # Copy numbers are non-negative and respect --cn_max.
    assert all(0 <= value <= 4 for value in cn.total_cn())


def test_haplotype_file_splits_copy_number_by_allele(tiny_sim: Simulation):
    # Haplotype-specific output is sample, chr, site, cnA, cnB -- five
    # columns, with no separate total. The total is cnA + cnB.
    cn = read_copy_numbers(tiny_sim.haplotype_cn)

    for row in cn.rows:
        assert len(row) == 5, f"expected 5 columns, got {len(row)}: {row}"
        cn_a, cn_b = row[3], row[4]
        assert cn_a >= 0 and cn_b >= 0, f"negative allele copy number in {row}"
        assert cn_a + cn_b <= 4, f"total CN exceeds --cn_max in {row}"


def test_haplotype_and_total_files_agree(tiny_sim: Simulation):
    # The two views of the same simulation must describe the same genome.
    total = read_copy_numbers(tiny_sim.cn)
    haplotype = read_copy_numbers(tiny_sim.haplotype_cn)

    assert len(total.rows) == len(haplotype.rows)

    for total_row, hap_row in zip(total.rows, haplotype.rows):
        assert total_row[:3] == hap_row[:3], "row ordering differs between files"
        assert total_row[3] == hap_row[3] + hap_row[4], (
            f"total CN {total_row[3]} != cnA + cnB for site {total_row[:3]}"
        )


def test_simulated_tree_is_a_rooted_binary_tree(tiny_sim: Simulation):
    tree = read_tree(tiny_sim.tree)

    assert tree.n_leaves == EXPECTED_SAMPLES
    assert tree.n_edges == 2 * EXPECTED_SAMPLES - 2
    assert all(edge.length >= 0 for edge in tree.edges)

    # Node IDs are 1..2n-1 with none skipped.
    assert tree.node_ids == set(range(1, 2 * EXPECTED_SAMPLES))


def test_sample_times_cover_every_tumour_region(tiny_sim: Simulation):
    lines = [line for line in tiny_sim.times.read_text().splitlines() if line.strip()]

    assert len(lines) == cli.NS

    for line in lines:
        fields = line.split()
        assert len(fields) >= 2, f"unexpected timing row: {line!r}"


def test_same_seed_reproduces_the_same_data(tmp_path, cnets, run):
    """The seed is the only reproducibility guarantee the tools offer.

    Compares decompressed bytes: gzip records a modification time in its
    header, so two identical payloads produce different .gz files.
    """
    digests = []
    for name in ("first", "second"):
        outdir = tmp_path / name
        outdir.mkdir()
        run([cnets, *cli.cnets_args(outdir, **{"--seed": 4242})]).assert_ok()

        payload = gzip.open(Simulation.at(outdir).cn, "rb").read()
        digests.append(hashlib.sha256(payload).hexdigest())

    assert digests[0] == digests[1], "same seed produced different copy numbers"


def test_different_seeds_produce_different_data(tmp_path, cnets, run):
    """Guards against the seed being ignored, which would make the test above
    pass for the wrong reason."""
    digests = []
    for seed in (4242, 9999):
        outdir = tmp_path / str(seed)
        outdir.mkdir()
        run([cnets, *cli.cnets_args(outdir, **{"--seed": seed})]).assert_ok()

        payload = gzip.open(Simulation.at(outdir).cn, "rb").read()
        digests.append(hashlib.sha256(payload).hexdigest())

    assert digests[0] != digests[1], "different seeds produced identical output"


def _cn_digest(outdir) -> str:
    """Digest of the decompressed copy numbers (see the gzip note above)."""
    payload = gzip.open(Simulation.at(outdir).cn, "rb").read()
    return hashlib.sha256(payload).hexdigest()


def _printed_seed(result) -> int:
    """The seed the program reports on stdout as `seed:<tab><number>`."""
    lines = [line for line in result.stdout.splitlines()
             if line.startswith("seed:\t")]
    assert len(lines) == 1, f"expected exactly one seed line\n{result}"
    return int(lines[0].split("\t")[1])


def test_the_seed_in_use_is_printed(tmp_path, cnets, run):
    """Every run prints the seed it used, so a log is enough to repeat it."""
    outdir = tmp_path / "out"
    outdir.mkdir()
    result = run([cnets, *cli.cnets_args(outdir, **{"--seed": 4242})]).assert_ok()

    assert _printed_seed(result) == 4242


def test_a_clock_seed_is_printed_and_can_be_replayed(tmp_path, cnets, run):
    """--seed 0, the programs' default, lets cnets pick a seed from the clock.
    That seed must be printed, and passing it back must reproduce the run."""
    first = tmp_path / "first"
    first.mkdir()
    result = run([cnets, *cli.cnets_args(first, **{"--seed": 0})]).assert_ok()
    seed = _printed_seed(result)
    assert seed > 0

    second = tmp_path / "second"
    second.mkdir()
    run([cnets, *cli.cnets_args(second, **{"--seed": seed})]).assert_ok()

    assert _cn_digest(first) == _cn_digest(second), (
        f"replaying the printed seed {seed} gave different data"
    )


def _info_text(outdir) -> str:
    return Simulation.at(outdir).info.read_text()


def test_seeds_are_recorded_in_the_info_file(tmp_path, cnets, run):
    """The record travels with the data: whatever launched cnets, the info
    file says which seeds produced the data set."""
    plain = tmp_path / "plain"
    plain.mkdir()
    run([cnets, *cli.cnets_args(plain, **{"--seed": 4242})]).assert_ok()
    assert "SEED: 4242\n" in _info_text(plain)
    assert "TREE SEED: same stream as SEED\n" in _info_text(plain)

    separate = tmp_path / "separate"
    separate.mkdir()
    run([cnets, *cli.cnets_args(separate, **{"--seed": 4242,
                                             "--tree_seed": 7})]).assert_ok()
    assert "SEED: 4242\n" in _info_text(separate)
    assert "TREE SEED: 7\n" in _info_text(separate)


def _edge_lists(outdir, replicates: int):
    """(start, end, length) per edge, per replicate; mutation counts ignored."""
    return [[(e.start, e.end, e.length)
             for e in read_tree(outdir / f"sim-data-{i}-tree.txt").edges]
            for i in range(1, replicates + 1)]


def test_tree_seed_gives_the_same_trees_across_rate_modes(tmp_path, cnets, run):
    """Branch-specific-rate modes draw rate multipliers from the stream, so
    with one shared stream the second replicate onwards gets different trees
    from a constant-rate run with the same --seed. A separate --tree_seed
    makes the trees independent of whatever is drawn after them."""
    replicates = 3

    def simulate(name, **overrides):
        outdir = tmp_path / name
        outdir.mkdir()
        # --method 0 (waiting times): the only path that applies --bsr_mode.
        run([cnets, *cli.cnets_args(outdir, **{"-n": replicates, "--method": 0,
                                               **overrides})]).assert_ok()
        return _edge_lists(outdir, replicates)

    constant = simulate("constant", **{"--bsr_mode": 0, "--tree_seed": 11})
    local_clock = simulate("local-clock", **{"--bsr_mode": 3, "--tree_seed": 11})
    assert constant == local_clock, "same --tree_seed, different trees"

    constant_shared = simulate("constant-shared", **{"--bsr_mode": 0})
    local_clock_shared = simulate("local-clock-shared", **{"--bsr_mode": 3})
    assert constant_shared[0] == local_clock_shared[0], (
        "the first replicate is drawn before any rate multiplier, so it "
        "should agree even on the shared stream"
    )
    assert constant_shared[1:] != local_clock_shared[1:], (
        "without --tree_seed the later replicates were expected to differ; "
        "if they agree now, the desync this option exists for has gone away"
    )


def test_bsr_mode_is_refused_when_sequences_are_simulated_directly(tmp_path,
                                                                    cnets, run):
    """Branch-specific rates are only drawn on the waiting-times path
    (--method 0). With --method 1, the suite's default, --bsr_mode used to be
    silently ignored and the run came out constant-rate. It must now stop
    with a message that names both options."""
    outdir = tmp_path / "out"
    outdir.mkdir()
    result = run([cnets, *cli.cnets_args(outdir, **{"--bsr_mode": 1})])

    assert result.returncode > 0, str(result)
    combined = result.stdout + result.stderr
    assert "bsr_mode" in combined and "method" in combined, str(result)


def test_help_lists_key_options(cnets, run):
    result = run([cnets, "--help"])

    combined = result.stdout + result.stderr
    for option in ("--cn_max", "--model", "--seed"):
        assert option in combined, f"{option} missing from --help output"


@pytest.mark.parametrize("program", ["cnets", "cnetml", "cnetmcmc"])
@pytest.mark.parametrize("flag", ["--help", "--version"])
def test_help_and_version_exit_successfully(program, flag, request, run):
    """Asking for help or the version is not an error. The three programs
    used to exit 1 after printing, which breaks `cmd --help` in shell
    pipelines, Makefiles and package tests."""
    result = run([request.getfixturevalue(program), flag])

    assert result.returncode == 0, str(result)
    assert result.stdout.strip(), f"{program} {flag} printed nothing"
