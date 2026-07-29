# Utility scripts

The `util/` directory contains R and Python helpers used around the main
tools. They are standalone scripts rather than installed packages.

:::{note}
Descriptions below come from each script's own header comments and
`option_list`/argument parser. Please review — they have not been
verified against actual runs.
:::

## R scripts

`preprocess_QDNASeq.R`
: Converts QDNAseq copy-number-calling output into `cnetml` input
  (`*-cn.txt`, `*-sample-ids.txt`, and optionally `*-rel-times.txt`). A
  library of functions (`transform.data`, `write_cn`, ...) rather than a
  standalone CLI — source it and call the functions you need.

`build_parsimony_tree.R`
: Builds maximum-parsimony NEWICK trees (via `phangorn`) to use as initial
  trees for `cnetml` tree search. Options: `-f/--file_cn`,
  `-b/--bootstrap`, `-i/--input_format`, `-o/--output_format`,
  `-d/--dir_nwk`, `-c/--file_bs`, `-m/--incl_normal`,
  `-a/--is_haplotype_specific`, `-n/--num_generate`, `-s/--num_select`.

`check_convergence.R`
: Checks MCMC chain convergence with
  [RWTY](https://github.com/danlwarren/RWTY). Positional arguments:
  `Rscript check_convergence.R <mcmc_dir> <output.pdf>`.

`compare_trees.R`
: A library of tree-comparison functions (Robinson-Foulds and related
  metrics via `ape`/`phangorn`/`Quartet`, plus MEDICC tree import). Its
  CLI option parsing is currently commented out in the source — it is
  sourced by other analysis scripts rather than run directly.

`compute_ci.R`
: Computes a confidence interval (normal approximation or 2.5th/97.5th
  percentile) from a two-column `ID value` file. Positional arguments:
  `Rscript compute_ci.R <file> <digits> <Y/N: normal approximation>`.

`plot-cns.R`, `plot-trees-all.R`
: Plot simulated/inferred copy-number profiles and trees (a single tree,
  all trees, or bootstrap trees with support values). `plot-util.R` is
  the shared helper library both source and is not run directly. See the
  `Rscript util/plot-trees-all.R ...` and
  `Rscript util/plot-cns.R -d $dir -b util/bin_locations_4401.Rdata`
  calls in `run-cnets.sh` / `run-cnetml.sh` for working examples.

`convert4medicc_fasta.R`, `convert4medicc_tsv.R`
: Convert simulated haplotype-specific copy numbers to MEDICC / MEDICC2
  input format. Example:
  `Rscript convert4medicc_fasta.R -i sim-data-1-allele-cn.txt -d medicc_input`.

`create_pseudo_subclone.R`
: Creates data with subclonal structure (purity-mixed samples) from
  `cnets` output; total copy number only. A library of functions, not a
  standalone CLI.

`check_site_pattern.R`
: Counts unique copy-number site patterns across samples — used to choose
  `max_site_change` for `cnetml` model 3. Example (from the original
  README):
  `Rscript util/check_site_pattern.R -c sim1-cn.txt -t sim1-patterns.txt`.

## Python scripts

`newick2elist.py`
: Converts NEWICK/NEXUS trees into the edge-list format read by the
  `svtree*` family of programs. Examples (from the script's own
  docstring):
  `python newick2elist.py -f 0 -t AllTreesNr5.txt -b 5.25` (NEXUS, may
  contain multiple trees) or `python newick2elist.py -f 1 -t tree1.nwk`
  (NEWICK, single tree).

## Reference data

`bin_locations_4401.Rdata`
: Bin-boundary table (chromosome/start/end for 4401 fixed-size bins), used
  by `preprocess_QDNASeq.R` and passed to `plot-cns.R` via `-b`.

`cytoBand_hg19.txt`, `cytoBand_hg38.txt`
: Standard UCSC cytoband tables for hg19/hg38, used by `plot-trees-all.R`
  via `--cyto_file` for chromosome-level plot annotation.
