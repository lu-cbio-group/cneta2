# Architecture

The current structure: independent command-line tools sharing C++ sources
under `code/`, before the planned `libcneta` core exists. See
[Repository layout](index.md#repository-layout) for the directory-level
view.

## How the main files fit together

```mermaid
flowchart LR
    A[cnets.cpp] --> B[tree_op.cpp and evo_tree.cpp]
    A --> C[model.cpp]
    A --> D[parse_cn.cpp]
    A --> E[Simulation outputs]
    F[cnetml.cpp] --> D
    F --> B
    F --> G[nni.cpp]
    F --> H[optimization.cpp]
    H --> I[likelihood.cpp]
    I --> C
    F --> J[Tree, Nexus, summary, and edge-rate outputs]
```

## `cnetml` likelihood calculation

The current detailed likelihood path is the DECOMP path, where WGD,
chromosome-level changes, and segment-level changes are represented by
separate chains. Constant-rate and branch-specific-rate modes differ mainly
in how transition probability matrices are built.

```mermaid
flowchart TD
    A([Likelihood request]) --> B[Validate tree against timing constraints]
    B -->|Invalid| C[Return very small likelihood]
    B -->|Valid| D[Set matrix dimensions from observed changes]
    D --> E{Branch-specific rates?}
    E -->|No| F[Build global Q rate matrices for WGD, chromosome, segment]
    F --> G[Build P transition matrices keyed by branch length]
    G --> H[Run pruning over chromosomes and sites]
    E -->|Yes| I[Use edge_rates stored on tree]
    I --> J[Build per-edge P transition matrices keyed by edge id]
    J --> K[Run variable-rate pruning over chromosomes and sites]
    H --> L[Combine site-pattern log likelihoods]
    K --> L
    L --> M{Acquisition-bias correction?}
    M -->|Yes| N[Compute likelihood of invariant bin]
    N --> O[Add invariant-bin contribution]
    M -->|No| P[Keep raw log likelihood]
    O --> Q[Clamp non-finite or too-small values]
    P --> Q
    Q --> R([Return log likelihood])
```

For branch-specific-rate modes, the same three biological levels are used,
but transition matrices are built per edge from `tree.edge_rates` using
`build_transition_matrices_variable`. The traversal functions are then
`get_likelihood_wgd_variable`, `get_likelihood_per_chr_variable`, and
`get_likelihood_site_change_variable`, and the final aggregation happens in
`get_likelihood_chr_change_variable`.

## Three-level DECOMP likelihood

In the DECOMP model, observed copy-number changes are decomposed into three
biological levels. Segment/site changes are evaluated per site, chromosome
changes are evaluated once per chromosome, and WGD is evaluated once per
sample tree because it is shared across all sites.

```mermaid
flowchart TD
    A["Observed copy-number changes (CN_CHANGE) per sample and site"] --> B[Split into three DECOMP levels]
    B --> C["WGD level (num_wgd)"]
    B --> D["Chromosome level (cn_change_chr)"]
    B --> E["Segment/site level (cn_change_site)"]
    C --> F["Set WGD matrix dimension (set_pmat_decomp_dim)"]
    D --> G["Set chromosome matrix dimension (set_pmat_decomp_dim)"]
    E --> H["Set segment matrix dimension (set_pmat_decomp_dim)"]
    F --> I["Build WGD Q and P matrices (get_rate_matrix_wgd)"]
    G --> J["Build chromosome Q and P matrices (get_rate_matrix_change_haplotype), using chr_gain and chr_loss"]
    H --> K["Build segment Q and P matrices (get_rate_matrix_change_haplotype), using dup and del"]
    I --> L["Initialize WGD likelihood table (initialize_lnl_table_wgd)"]
    J --> M["Initialize chromosome likelihood table (initialize_lnl_table_chr)"]
    K --> N["Initialize site likelihood table (initialize_lnl_table_site)"]
    L --> O["Tree pruning for WGD (get_likelihood_wgd)"]
    M --> P["Tree pruning for chromosome gain/loss (get_likelihood_per_chr)"]
    N --> Q["Tree pruning for segment changes (get_likelihood_site_change)"]
    O --> R[Extract WGD log likelihood]
    P --> S[Extract chromosome log likelihood]
    Q --> T[Extract segment log likelihood]
    R --> U[Sum level-specific log likelihoods]
    S --> U
    T --> U
    U --> V[get_likelihood_chr_change]
    V --> W([Return DECOMP log likelihood])
```

## `cnetml` candidate tree optimization

For a fixed topology, `cnetml` optimizes branch lengths or node-time ratios,
and optionally rate parameters. With branch-specific-rate inference, the
optimizer updates `tree.edge_rates` before each likelihood evaluation.

```mermaid
flowchart TD
    A([Candidate tree]) --> B[Compute internal-node traversal order]
    B --> C{Optimizer}
    C -->|Simplex| D[Run GSL maximization]
    C -->|L-BFGS-B| E[Select optimization variables]
    E --> F{Constrained by sampling time?}
    F -->|Yes| G[Optimize transformed node-time ratios]
    F -->|No| H[Optimize branch lengths directly]
    G --> I{bsr_mode}
    H --> I
    I -->|0| J[Use constant global rates]
    I -->|1| K[Estimate one shared multiplier per active edge]
    I -->|2| L[Estimate event-specific multipliers per active edge]
    I -->|3| M[Estimate multipliers only for selected RLC shift edges]
    J --> N[Evaluate targetFunk]
    K --> N
    L --> N
    M --> N
    N --> O[Update tree variables from optimizer vector]
    O --> P{DECOMP model?}
    P -->|Yes, constant rate| Q[get_likelihood_change]
    P -->|Yes, variable rate| R[get_likelihood_change_variable_rate]
    P -->|No| S[get_likelihood_revised]
    Q --> T[Return negative log likelihood]
    R --> T
    S --> T
    T --> U[Optimizer iterates until convergence or limit]
    U --> V[Store optimized tree score]
    V --> W([Optimized candidate])
```

## Random local clock shift-edge search

When `--bsr_mode 3` is used with `mode 3`, CNETML performs an outer greedy
search over shift edges. The score is `logL - rlc_lambda * K`, where `K` is
the number of selected shift edges.

```mermaid
flowchart TD
    A([Start ML-RLC]) --> B[Initialize selected_edges as empty]
    B --> C[Optimize baseline model with no shift edges]
    C --> D[Compute current penalized score]
    D --> E[Build candidate pool from optimizable edges]
    E --> F{Any unselected candidate edge?}
    F -->|No| G[Finish with current selected_edges]
    F -->|Yes| H[Add one candidate edge to selected_edges]
    H --> I[Warm-start candidate shift multipliers at 1]
    I --> J[Run L-BFGS-B for branch lengths and selected shift multipliers]
    J --> K[Compute candidate score: logL - rlc_lambda * K]
    K --> L{Best candidate improves score?}
    L -->|Check more candidates| F
    L -->|No improvement after all candidates| G
    L -->|Yes| M[Accept best edge]
    M --> N[Update current tree, selected_edges, raw logL, and penalized score]
    N --> F
    G --> O[Write rlc_shift_eids, raw_logL, penalized_score]
    O --> P[Write edge-rate table with local clock IDs]
    P --> Q([End])
```

:::{note}
These five diagrams come from the research team and cover `cnets` and
`cnetml`; `cnetmcmc` is intentionally not diagrammed yet. The top-level
`cnets` and `cnetml` workflows live in [Workflows](../workflows/index.md)
since they describe user-facing data flow rather than internals.
:::
