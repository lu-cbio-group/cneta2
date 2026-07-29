# Quick start

A minimal end-to-end run: simulate data with `cnets`, then infer a tree
from it with `cnetml`.

## 1. Build

```bash
cd code && ./build.sh local
```

See [Installation](../installation/index.md) for HPC targets.

## 2. Simulate a small dataset

```bash
./run-cnets.sh
```

Edit the parameter block at the top of the script rather than passing
flags by hand. The defaults simulate `Ns=3` tumour regions plus the normal
sample, using the haplotype-specific model (`model=2`), `cn_max=4`, with
only duplication/deletion events (`r3=r4=r5=0`, i.e. no chromosome
gain/loss or WGD). Output goes to `./example/`.

## 3. Infer a tree

```bash
./run-cnetml.sh
```

By default this reads `./example/sim-data-1-cn.txt.gz` — note the matching
`prefix=sim-data-1` at the top of the script. If you changed the prefix (or
left it empty) in `run-cnets.sh`, update `prefix` here to match. `mode=0`
(the default) runs maximum-likelihood tree search.

## 4. Understand the outputs

`run-cnets.sh` writes `*-cn.txt.gz` (simulated copy numbers), `*-tree.txt`
/ `*-tree.nex` (the true tree), and `*-info.txt` / `*-mut.txt` (mutations
per branch) into `./example/`.

`run-cnetml.sh` writes the reconstructed tree to `MaxL-<suffix>.txt` in the
same directory, plus a `std_cnetml_<suffix>` log file — check it for
`cnetml main run SUCCEEDED` or `FAILED`.

See [File formats](../file-formats/index.md) for the column-level detail of
each file, and the [User guide](../user-guide/index.md) for every option
`cnets` and `cnetml` accept.
