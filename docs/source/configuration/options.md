# Command-line options reference

Every command-line option of `cnets` and `cnetml`, with its short form,
default value, accepted values and meaning. Defaults and accepted values
are taken from the `boost::program_options` definitions in
`code/cnets.cpp` and `code/cnetml.cpp`. For an explanation of *why* you
would pick one setting over another, see each tool's User guide page
([`cnets`](../user-guide/cnets.md#options),
[`cnetml`](../user-guide/cnetml.md#options)). `cnetmcmc` reads most of
its settings from a config file; see its
[Options](../user-guide/cnetmcmc.md#options).

Usage: `bin/cnets --name value` or `bin/cnets -x value`, for example
`bin/cnets -o sim/ --nregion 5 --bsr_mode 3 --seed 1`. Run any tool with
`--help` to print its full option list.

On/off switches take an integer: `0` = no, `1` = yes.

## `cnets` (simulation)

### General and required

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--help` | `-h` | — | — | Print all options |
| `--version` | `-v` | — | — | Print version |
| `--odir` | `-o` | **required** | directory path | Output directory (trailing `/` optional) |

### Run mode and model

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--mode` | | 0 | 0 / 1 / 2–4 | 0: simulate the genome in fixed-size bins (4401 bins of 500 Kbp by default); 1: simulate segments of variable size, each site being a final segment; 2–4: internal tests |
| `--method` | | 0 | 0 / 1 | 0: simulate waiting times (supports all event types); 1: simulate sequences directly at branch ends from the substitution model (duplication/deletion only; forces `mode=1`) |
| `--model` | `-d` | 2 | 0 / 1 / 2 / 3 | Model of evolution. 0: Mk; 1: one-step bounded (total); 2: one-step bounded (haplotype-specific); 3: infinite sites (Poisson process) |
| `--nsim` | `-n` | 3 | positive integer | Number of independent multi-region datasets to simulate (output `sim-data-1` … `sim-data-N`) |
| `--nregion` | `-r` | 5 | positive integer | Number of samples (regions) per dataset, i.e. the number of tips |
| `--tree_file` | | "" | file path | Simulate mutations along this tree; if empty, a random coalescent tree is generated |

### Tree generation (random coalescent tree)

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--epop` | `-e` | 2 | positive integer | Effective population size Ne of cell populations, used to convert coalescent time to years |
| `--gtime` | | 0.002739726 | > 0 | Generation time in years (default = 1 day) |
| `--beta` | `-b` | 0.0 | ≥ 0 | Population growth rate; > 0 simulates exponential growth |
| `--age` | `-a` | 100 | positive integer | Patient age in years |
| `--constrained` | | 1 | 0 / 1 | 1: scale initial branch lengths by `--age` so that the tree height is below the patient age |
| `--tdiff` | `-t` | 0 | ≥ 0 | Baseline sampling-time difference; if > 0, each sample's time = earliest sample + a random multiple of `tdiff` |
| `--stime` | | "" | space-separated integers | Sampling time of each sample in whole years, with the most recent sample at 0; give all values in one quoted string, e.g. `--stime "0 1 2 3 4"` |
| `--tree_seed` | | 0 | non-negative integer | Seed of a separate generator used only to draw trees. 0: share the `--seed` stream. Use it to get the same trees in runs that differ only in later settings such as `--bsr_mode` |

### Genome and sites

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--cn_max` | | 4 | positive integer | Maximum copy number of a site |
| `--seg_max` | | 100 | ≥ number of chromosomes | Maximum number of sites to simulate (used when `mode=1`) |
| `--fix_nseg` | | 1 | 0 / 1 | 1: number of sites fixed at `seg_max`; 0: drawn at random between 0.2 and 0.8 × `seg_max` |

### Mutation rates and sizes

All rates are per year. Set a rate to 0 to leave that event type out of
the simulation (`cnets` has no `--cn_type` option).

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--dup_rate` | | 0.001 | ≥ 0 | Site-level duplication rate (per haplotype per site per year) |
| `--del_rate` | | 0.001 | ≥ 0 | Site-level deletion rate (per haplotype per site per year) |
| `--chr_gain` | | 0.0 | ≥ 0 | Chromosome gain rate (per haplotype per chromosome per year) |
| `--chr_loss` | | 0.0 | ≥ 0 | Chromosome loss rate (per haplotype per chromosome per year) |
| `--wgd` | | 0.0 | ≥ 0 | Whole-genome doubling rate (per year) |
| `--dup_size` | | 1 | > 0 | Mean duplication size in bins; forced to 1 when `mode=1` |
| `--del_size` | | 1 | > 0 | Mean deletion size in bins; forced to 1 when `mode=1` |
| `--error_rate` | | 0 | 0–1 | Fraction of sites per sample whose copy number is randomly altered |

### Branch-specific rates

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--bsr_mode` | | 0 | 0 / 1 / 2 / 3 | Only with `method=0`; with `method=1`, `cnets` stops with an error. 0: one rate on all branches; 1: one shared multiplier per branch (all event types); 2: independent rate per branch and event type; 3: random local clock (rates inherited top-down, changing at each node with probability `bsr_p`) |
| `--bsr_dist` | | 0 | 0 / 1 | Distribution of the multipliers. 0: log-normal; 1: gamma |
| `--bsr_variance` | | 0.25 | > 0 | Target variance of multiplier m (E[m] = 1). Log-normal: σ = √log(1 + var), ζ = −σ²/2; gamma: shape = 1/var, scale = var |
| `--bsr_p` | | 0.3 | (0, 1) | `bsr_mode=3` only: probability of a rate change at each node |

### Output

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--prefix` | `-p` | "" | string | Output file prefix; `sim-data-N` if empty |
| `--print_allele` | | 1 | 0 / 1 | Write haplotype-specific copy numbers |
| `--print_mut` | | 1 | 0 / 1 | Write the list of mutations |
| `--print_nex` | | 1 | 0 / 1 | Write the tree in NEXUS format |
| `--print_relative` | | 0 | 0 / 1 | Write relative total copy numbers, as in CGHcall: −2 double deletion, −1 single deletion, 0 normal, 1 gain, 2 double gain |
| `--use_nwgd` | | 0 | 0 / 1 | Normalise by the known number of WGDs when computing relative copy numbers (with `print_relative=1`) |
| `--random_round` | | 0 | 0 / 1 | Use random rounding when computing relative copy numbers (with `print_relative=1`) |

### Other

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--seed` | | 0 | non-negative integer | Random seed; 0 = generated automatically |
| `--verbose` | | 0 | 0 / 1 / 2 | 0: default; 1: debug; 2: debug with tree details |

### Notes on `cnets` behaviour

- With `method=0` (waiting times), `model=1` is changed to `model=2`.
- The number of erroneous sites is round(number of sites × `error_rate`),
  and the number of sites is only set when `mode=1`. With `mode=0`,
  `--error_rate` therefore has no effect.
- `--stime` takes whole years only, separated by spaces. A comma or a
  decimal point stops the reading there, and the remaining samples get
  wrong times without a warning.
- `cnets` exits if `seg_max` is smaller than the number of chromosomes.
- `bin/run-cnets.sh` names the tree-height constraint `cons`; the
  `cnets` flag is `--constrained`.

## `cnetml` (maximum-likelihood inference)

### General, input and output

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--help` | `-h` | — | — | Print all options |
| `--version` | `-v` | — | — | Print version |
| `--cfile` | `-c` | **required** | file path | Input copy-number file (e.g. `*-cn.txt.gz` / `*-haplotype-cn.txt.gz` from `cnets`). A missing file can cause a segmentation fault |
| `--tfile` | `-t` | "" | file path | Sampling-time file (e.g. `*-rel-times.txt`); if empty, all samples have time 0 |
| `--meta_file` | | "" | file path | Meta information: known number of WGDs (required if the file is given), optionally per-chromosome gain/loss and per-site duplication/deletion |
| `--nsample` | `-s` | 5 | positive integer | Number of samples; must match the input file |
| `--tree_file` | | "" | file path | Given tree (needed for `mode` 1–4) |
| `--ofile` | `-o` | `maxL-tree.txt` | file path | Output maximum-likelihood tree; other outputs are `<ofile>.nex`, `.nmut.nex`, `.summary.txt`, `.edge_rates.txt`, etc. |
| `--seg_file` | | "" | file path | Write the preprocessed copy-number matrix used for tree building |

### Run mode

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--mode` | | 0 | 0–5 | 0: infer the ML tree from copy numbers; 1: test on example data; 2: likelihood of a given tree with branch lengths; 3: given topology, optimise branch lengths (and rates) for the maximum likelihood; 4: infer ancestral states of a given tree; 5: write the segment file and exit |
| `--bootstrap` | `-b` | 0 | 0 / 1 | Bootstrap resampling of sites (for mode 0 / 3) |

### Model and constraints

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--model` | `-d` | 2 | 0 / 1 / 2 / 3 | 0: Mk (deprecated); 1: one-step bounded, total (deprecated); 2: one-step bounded, haplotype-specific; 3: independent Markov chains (recommended for data with chromosome gain/loss and WGD) |
| `--constrained` | | 1 | 0 / 1 | Branch-length constraint. 0: none; 1: fixed total time (constrained by age and sampling times) |
| `--estmu` | `-u` | 0 | 0 / 1 | 0: rates fixed at the input values; 1: estimate a constant mutation rate. Reliable only with `constrained=1` and differing sampling times |
| `--bsr_mode` | | 0 | 0 / 1 / 2 / 3 | Branch-specific rate inference. 0: constant; 1: one shared multiplier per branch; 2: independent rate per branch and event type; 3: random local clock (ML-RLC). Values > 0 require `model=3`; otherwise `cnetml` stops with an error |
| `--optim` | | 1 | 0 / 1 | Optimiser. 0: Simplex; 1: L-BFGS-B |
| `--cn_type` | | **required** | 0–4 | Event types whose rates are estimated. 0: segment level only (dup/del); 1: chromosome gain/loss + WGD; 2: dup/del + WGD; 3: dup/del + chromosome gain/loss; 4: all. Required; `cnetml` stops with an error if it is missing |

### Copy-number limits (rate-matrix dimensions)

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--cn_max` | | 4 | positive integer | Maximum copy number allowed by the bounded models |
| `--max_wgd` | | 1 | ≥ 0 | Maximum number of WGDs (WGD rate-matrix dimension, model 3) |
| `--max_chr_change_haplotype` | | 1 | ≥ 0 | Maximum chromosome gains/losses per haplotype (model 3) |
| `--max_site_change_haplotype` | | 2 | ≥ 0 | Maximum site duplications/deletions per haplotype (model 3). Too low gives an incorrect likelihood; too high slows computation |
| `--m_max` | | 1 | ≥ 0 | Maximum segment copies before chromosome-level events (currently unused; set to 0 when `cn_type=2`) |

### Input data type

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--is_total` | | 1 | 0 / 1 | Input is total copy number (0 = haplotype-specific) |
| `--is_rcn` | | 0 | 0 / 1 | Input is relative copy number; if 1, `cn_type` is forced to 0 |
| `--is_bin` | | 1 | 0 / 1 | Input is given per bin. 1: merge consecutive identical bins; 0: read as is |
| `--incl_all` | | 1 | 0 / 1 | Use all input segments, including invariant sites. With `is_bin=0` and `incl_all=1`, `correct_bias` is forced to 0 |

### Likelihood computation

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--use_repeat` | | 1 | 0 / 1 | Merge repeated site patterns to speed up the likelihood |
| `--correct_bias` | | 1 | 0 / 1 | Apply acquisition-bias correction when invariant sites are excluded |

### Ancestral states (mode 4)

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--infer_marginal_state` | | 1 | 0 / 1 | Infer the marginal ancestral state of the MRCA |
| `--infer_joint_state` | | 1 | 0 / 1 | Infer the joint ancestral states of all internal nodes |

### Initial trees

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--init_tree` | | 0 | 0 / 1 | Initial tree method. 0: random coalescent tree; 1: maximum parsimony tree |
| `--dir_itrees` | | "" | directory path | Directory of user-supplied initial trees (recommended for large trees, to avoid local optima) |
| `--epop` | | 2 | positive integer | Effective population size for random initial trees |
| `--gtime` | | 1 | > 0 | Generation time in years for random initial trees (differs from the `cnets` default) |
| `--beta` | | 0 | ≥ 0 | Population growth rate for random initial trees |

### Tree search and optimisation

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--tree_search` | | 1 | 0 / 1 / 2 | 0: genetic algorithm (slow, deprecated); 1: random-restart hill climbing (≥ 5 samples); 2: exhaustive search (efficient for ≤ 7 samples) |
| `--speed_nni` | | 1 | 0 / 1 | Use reduced NNI during hill climbing |
| `--npop` | `-p` | 100 | positive integer | Genetic-algorithm population size, or maximum number of initial trees |
| `--ngen` | `-g` | 50 | positive integer | Genetic-algorithm generations, or maximum perturbations/optimisations per tree |
| `--nstop` | `-e` | 10 | positive integer | Stop after this many rounds without improvement |
| `--tolerance` | `-r` | 0.01 | > 0 | Convergence tolerance of the optimiser |
| `--miter` | `-m` | 2000 | positive integer | Maximum optimiser iterations |
| `--loglh_epsilon` | | 0.001 | > 0 | Tolerance for treating two log-likelihoods as equal |
| `--ssize` | `-z` | 0.01 | > 0 | Initial step size for GSL optimisation |
| `--scale_tobs` | | 1.0 | > 0 | Scale factor on the maximum sampling-time difference giving the lower bound of root age in constrained (L-BFGS-B) optimisation |

### Mutation rates (initial or fixed values)

Fixed values when `estmu=0`; starting values when `estmu=1` or
`bsr_mode>0`. Units are not specified in the code.

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--mu` | `-x` | 0.02 | > 0 | Overall mutation rate (overridden by 1 / number of sites under model 0) |
| `--dup_rate` | | 0.001 | ≥ 0 | Duplication rate |
| `--del_rate` | | 0.001 | ≥ 0 | Deletion rate |
| `--chr_gain_rate` | | 0 | ≥ 0 | Chromosome gain rate |
| `--chr_loss_rate` | | 0 | ≥ 0 | Chromosome loss rate |
| `--wgd_rate` | | 0 | ≥ 0 | WGD rate |

### Branch-specific rates (`bsr_mode > 0`)

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--estimate_bsr0_first` | | 0 | 0 / 1 | 0: skip calibration and estimate the trunk edge's (root → MRCA) own rate jointly with the per-branch multipliers/shifts; 1: calibrate the reference rates with `bsr_mode=0` first, then estimate per-branch multipliers relative to them |
| `--rlc_search_method` | | 1 | 0 / 1 / 2 | `bsr_mode=3` only: outer search over shift edges. 0: exhaustive over all subsets (small trees only); 1: stepwise greedy; 2: genetic algorithm (not implemented; exits with an error) |
| `--rlc_criterion` | | 0 | 0 / 1 | `bsr_mode=3` only: model-selection criterion. 0: AIC; 1: BIC |

### Other

| Option | Short | Default | Values | Description |
|---|---|---|---|---|
| `--seed` | | 0 | non-negative integer | Random seed; 0 = generated automatically |
| `--verbose` | | 0 | 0 / 1 / 2 | 0: default; 1: debug; > 1 prints more detail on the input data |

### Option combinations checked by `cnetml`

| Condition | Result |
|---|---|
| `bsr_mode>0` and `optim≠1` | Exits with an error |
| `bsr_mode>0` and `estmu=1` | Exits with an error |
| `bsr_mode>0` and `mode≠3` | Exits with an error (branch-specific rates are supported only for a given topology, mode 3) |
| `rlc_search_method=2` | Exits with an error (not implemented) |
| `optim>1` or `model>3` | Assertion failure |
| `is_rcn=1` | `cn_type` set to 0 |
| `is_bin=0` and `incl_all=1` | `correct_bias` set to 0 |
| `model=0` | `mu` set to 1 / number of sites |
| `mode=5` | Writes `seg_file`, then exits |

## Same name, different meaning in `cnets` and `cnetml`

| Option | `cnets` | `cnetml` |
|---|---|---|
| `-b` | `--beta`, population growth rate | `--bootstrap` |
| `-e` | `--epop`, effective population size | `--nstop`, stopping rule |
| `-p` | `--prefix`, output prefix | `--npop`, population size |
| `-r` | `--nregion`, number of samples | `--tolerance` |
| `-t` | `--tdiff`, sampling-time difference | `--tfile`, time file |
| `-o` | `--odir`, output directory | `--ofile`, output tree file |
| `--model 3` | infinite sites | independent Markov chains |
| `--gtime` default | 0.002739726 years (1 day) | 1 year |
| Chromosome-level rate options | `--chr_gain` / `--chr_loss` / `--wgd` | `--chr_gain_rate` / `--chr_loss_rate` / `--wgd_rate` |
| Number of samples | `--nregion` / `-r` | `--nsample` / `-s` |
