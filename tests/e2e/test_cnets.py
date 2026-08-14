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


def test_help_lists_key_options(cnets, run):
    result = run([cnets, "--help"])

    combined = result.stdout + result.stderr
    for option in ("--cn_max", "--model", "--seed"):
        assert option in combined, f"{option} missing from --help output"


@pytest.mark.xfail(
    reason="cnets and cnetml exit 1 after printing --help; asking for help "
           "is not an error, and a non-zero status breaks `cmd --help` in "
           "shell pipelines and Makefiles",
    strict=False,
)
def test_help_exits_successfully(cnets, run):
    assert run([cnets, "--help"]).returncode == 0
