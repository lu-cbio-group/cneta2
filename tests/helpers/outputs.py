"""Readers for the file formats the cneta programs produce.

Deliberately dependency-free: the standard library is enough for every format
here, so the test suite installs nothing beyond pytest itself.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path


# --------------------------------------------------------------------------
# Copy-number matrices (*-cn.txt.gz, *-haplotype-cn.txt.gz)
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class CopyNumbers:
    """A parsed copy-number file.

    Columns are sample_ID, chr_ID, site_ID, CN -- and for haplotype-specific
    files, cnA and cnB after that. There is no header.
    """

    rows: list[tuple[int, ...]]

    @property
    def samples(self) -> set[int]:
        return {row[0] for row in self.rows}

    @property
    def sites_per_sample(self) -> dict[int, int]:
        counts: dict[int, int] = {}
        for row in self.rows:
            counts[row[0]] = counts.get(row[0], 0) + 1
        return counts

    def total_cn(self) -> list[int]:
        return [row[3] for row in self.rows]


def read_copy_numbers(path: Path) -> CopyNumbers:
    """Read a copy-number file, gzipped or not, and check it is all integers."""
    opener = gzip.open if str(path).endswith(".gz") else open

    rows: list[tuple[int, ...]] = []
    with opener(path, "rt") as handle:
        for lineno, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) < 4:
                raise AssertionError(
                    f"{path}:{lineno}: expected at least 4 columns, got {len(fields)}: {line!r}"
                )
            try:
                rows.append(tuple(int(f) for f in fields))
            except ValueError as exc:
                raise AssertionError(
                    f"{path}:{lineno}: every column must be an integer: {line!r}"
                ) from exc

    if not rows:
        raise AssertionError(f"{path}: file is empty")

    return CopyNumbers(rows)


# --------------------------------------------------------------------------
# Tree files (*-tree.txt, MaxL-*.txt)
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Edge:
    start: int
    end: int
    length: float


@dataclass(frozen=True)
class Tree:
    """An edge-list tree, as written by cnets and cnetml.

    Tab-separated with a header row: start, end, length, then some mix of
    eid and nmut depending on which program wrote it.
    """

    edges: list[Edge]

    @property
    def n_edges(self) -> int:
        return len(self.edges)

    @property
    def n_leaves(self) -> int:
        # A rooted binary tree on n leaves has 2n-2 edges.
        return (self.n_edges + 2) // 2

    @property
    def node_ids(self) -> set[int]:
        return {e.start for e in self.edges} | {e.end for e in self.edges}


def read_tree(path: Path) -> Tree:
    edges: list[Edge] = []

    with open(path) as handle:
        for lineno, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            fields = line.split()
            if fields[0] == "start":       # header
                continue
            if len(fields) < 3:
                raise AssertionError(
                    f"{path}:{lineno}: expected at least 3 columns, got {line!r}"
                )
            edges.append(Edge(int(fields[0]), int(fields[1]), float(fields[2])))

    if not edges:
        raise AssertionError(f"{path}: no edges found")

    return Tree(edges)


# --------------------------------------------------------------------------
# cnetml run summary (<output>.summary.txt)
# --------------------------------------------------------------------------

def read_summary(path: Path) -> dict[str, str]:
    """Read the tab-separated key/value summary cnetml writes after a search.

    Values stay as strings; "NA" is common and callers decide what it means.
    """
    summary: dict[str, str] = {}

    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                summary[fields[0]] = fields[1]

    if not summary:
        raise AssertionError(f"{path}: no key/value pairs found")

    return summary


# --------------------------------------------------------------------------
# MCMC traces
# --------------------------------------------------------------------------

def read_mcmc_params(path: Path) -> tuple[list[str], list[list[float]]]:
    """Read a ``.p`` trace: a ``# Parameters`` banner, a header, then samples."""
    header: list[str] = []
    samples: list[list[float]] = []

    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if not header:
                header = fields
                continue
            samples.append([float(f) for f in fields])

    return header, samples


NEWICK_TREE_RE = re.compile(r"^\s*tree\s+\S+\s*=\s*(?P<newick>.+;)\s*$")


def read_mcmc_trees(path: Path) -> list[str]:
    """Pull the Newick strings out of a ``.t`` trace file."""
    trees: list[str] = []

    with open(path) as handle:
        for line in handle:
            match = NEWICK_TREE_RE.match(line)
            if match:
                trees.append(match.group("newick"))

    return trees


def is_balanced_newick(newick: str) -> bool:
    """Cheap structural check: brackets balance and the string is terminated."""
    if not newick.endswith(";"):
        return False

    depth = 0
    for char in newick:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth < 0:
                return False

    return depth == 0


# --------------------------------------------------------------------------
# Program output
# --------------------------------------------------------------------------

LOGL_RE = re.compile(r"log likelihood of the input tree is\s+(?P<value>-?[\d.eE+]+)")
MIN_NLOGL_RE = re.compile(r"MIN -ve logL\s*=\s*(?P<value>-?[\d.eE+]+)")


def parse_scored_logl(stdout: str) -> float:
    """Extract the log-likelihood printed by ``cnetml --mode 2``."""
    match = LOGL_RE.search(stdout)
    if not match:
        raise AssertionError(
            "no log-likelihood line in cnetml output:\n" + stdout[-2000:]
        )
    return float(match.group("value"))


def parse_search_neg_logl(stdout: str) -> float:
    """Extract the best score printed at the end of a ``--mode 0`` search.

    Note the sign: cnetml prints "MIN -ve logL", the *negated* log-likelihood,
    so a well-fitting tree gives a positive number here. The ``raw_logL`` in
    the summary file is the same quantity with the opposite sign.
    """
    match = MIN_NLOGL_RE.search(stdout)
    if not match:
        raise AssertionError(
            "no 'MIN -ve logL' line in cnetml output:\n" + stdout[-2000:]
        )
    return float(match.group("value"))
