"""End-to-end tests for cnetml ancestral state reconstruction (``--mode 4``)
with branch-specific rates read from ``--edge_rate_file``.

The numerical side (per-edge transition matrices, agreement with the fitting
likelihood) is covered in tests/unit/test_state.cpp. These tests check the
command line: the rate table is read, checked against the tree, and refused
when it cannot be trusted.
"""

from __future__ import annotations

import pytest

from helpers import cli

RATE_COLUMNS = ["dup_rate", "del_rate", "chr_gain_rate", "chr_loss_rate", "wgd_rate"]
GLOBAL_RATES = {"dup_rate": 0.001, "del_rate": 0.001,
                "chr_gain_rate": 0, "chr_loss_rate": 0, "wgd_rate": 0}

OUTPUTS = [".mrca.cn", ".mrca.seg.state", ".joint.cn", ".joint.state"]


def tree_edges(tree_file):
    """(start, end, length) for each edge, as written in a cneta tree file."""
    lines = tree_file.read_text().splitlines()
    header = lines[0].split("\t")
    i_start, i_end, i_len = (header.index(c) for c in ("start", "end", "length"))
    return [(f[i_start], f[i_end], f[i_len]) for f in (l.split("\t") for l in lines[1:] if l)]


def write_rate_table(path, edges, rates_for_edge=lambda k: GLOBAL_RATES):
    """Write a rate table with the columns cnetml reads (extra columns are ignored)."""
    header = ["eid", "start", "end", "length", *RATE_COLUMNS]
    rows = ["\t".join(header)]
    for k, (start, end, length) in enumerate(edges):
        rates = rates_for_edge(k)
        rows.append("\t".join([str(k), start, end, length, *(str(rates[c]) for c in RATE_COLUMNS)]))
    path.write_text("\n".join(rows) + "\n")
    return path


@pytest.fixture
def asr(data_dir, cnetml, run, tmp_path):
    """Run mode 4 on the fixture; returns (result, output prefix)."""
    def _asr(name, **overrides):
        out = tmp_path / name
        result = run([cnetml, *cli.cnetml_asr_args(
            cn_file=data_dir / "tiny-cn.txt.gz",
            times_file=data_dir / "tiny-rel-times.txt",
            tree_file=data_dir / "tiny-tree.txt",
            out_file=out,
            **overrides,
        )])
        return result, out
    return _asr


@pytest.fixture
def edges(data_dir):
    return tree_edges(data_dir / "tiny-tree.txt")


def test_global_rates_on_every_edge_reproduce_the_plain_run(asr, edges, tmp_path):
    """A rate table that repeats the global rates must not change any output."""
    _, plain = asr("plain")
    table = write_rate_table(tmp_path / "rates.txt", edges)
    result, with_table = asr("with_table", **{"--edge_rate_file": table})
    result.assert_ok()
    assert "Reading branch-specific rates" in result.stdout

    for suffix in OUTPUTS:
        assert (plain.parent / (plain.name + suffix)).read_bytes() == \
               (with_table.parent / (with_table.name + suffix)).read_bytes(), suffix


def test_different_edge_rates_change_the_reconstruction(asr, edges, tmp_path):
    _, plain = asr("plain")
    faster = {**GLOBAL_RATES, "dup_rate": 0.02, "del_rate": 0.02}
    table = write_rate_table(tmp_path / "rates.txt", edges,
                             lambda k: faster if k % 2 else GLOBAL_RATES)
    asr("varied", **{"--edge_rate_file": table})[0].assert_ok()

    plain_states = (plain.parent / (plain.name + ".mrca.seg.state")).read_text()
    varied_states = (tmp_path / "varied.mrca.seg.state").read_text()
    assert plain_states != varied_states


@pytest.mark.parametrize("overrides, message", [
    ({"-d": 2}, "requires model 3"),
    ({"--mode": 3}, "only used with mode=4"),
])
def test_edge_rate_file_is_refused_outside_mode_4_model_3(asr, edges, tmp_path,
                                                          overrides, message):
    table = write_rate_table(tmp_path / "rates.txt", edges)
    result, _ = asr("refused", **{"--edge_rate_file": table, **overrides})
    assert result.returncode != 0
    assert message in result.stderr


def _drop_column(text):
    lines = text.splitlines()
    i = lines[0].split("\t").index("wgd_rate")
    return "\n".join("\t".join(f for j, f in enumerate(l.split("\t")) if j != i)
                     for l in lines) + "\n"


def _set_field(text, row, column, value):
    lines = text.splitlines()
    i = lines[0].split("\t").index(column)
    fields = lines[row].split("\t")
    fields[i] = value
    lines[row] = "\t".join(fields)
    return "\n".join(lines) + "\n"


@pytest.mark.parametrize("corrupt, message", [
    (_drop_column, "has no column 'wgd_rate'"),
    (lambda t: _set_field(t, 1, "dup_rate", "-0.1"), "is negative"),
    (lambda t: _set_field(t, 2, "del_rate", "abc"), "is not a number"),
    (lambda t: "\n".join(t.splitlines()[:-1]) + "\n", "is missing from"),
    (lambda t: t + t.splitlines()[1] + "\n", "appears more than once"),
    (lambda t: _set_field(t, 1, "length", "99"), "must come from the same run"),
    (lambda t: _set_field(t, 1, "end", "42"), "is not in the tree"),
    (lambda t: _set_field(t, 1, "start", "6.5"), "is not a node ID"),
])
def test_a_bad_rate_table_is_refused(asr, edges, tmp_path, corrupt, message):
    table = write_rate_table(tmp_path / "rates.txt", edges)
    table.write_text(corrupt(table.read_text()))
    result, _ = asr("bad", **{"--edge_rate_file": table})
    assert result.returncode != 0
    assert message in result.stderr
