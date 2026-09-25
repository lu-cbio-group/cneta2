# Workflows and data flow

How the tools are chained together, and what data passes between them.

:::{note}
The `cnets` and `cnetml` diagrams below are the reviewed ones — those
tools are published and current priority. The [MCMC workflow](#mcmc-workflow-cnetmcmc)
diagram for `cnetmcmc` is a draft, since that tool isn't officially
released yet.
:::

## Simulation workflow (`cnets`)

`cnets` simulates copy-number evolution along a tree and writes simulated
copy-number profiles, timing information, mutation lists, and tree files.

```mermaid
flowchart TD
    A([Start cnets]) --> B[Parse command-line options]
    B --> C[Validate simulation settings]
    C --> D[Initialize random-number generator]
    D --> E[Read optional sampling times]
    E --> F[Configure output options]
    F --> G{Simulation method}
    G -->|Waiting times| H[Use haplotype-specific CNA process]
    G -->|Sequences directly| I[Treat each site as a final segment]
    I --> J[Force segment mode]
    H --> K[Build chromosome or segment layout]
    J --> K
    K --> L[Set global mutation rates: dup, del, chr gain, chr loss, WGD]
    L --> M{Branch-specific rates?}
    M -->|bsr_mode 0| N[Use one global rate set on all branches]
    M -->|bsr_mode 1| O[Draw one multiplier per branch]
    M -->|bsr_mode 2| P[Draw event-specific rates per branch]
    M -->|bsr_mode 3| Q[Random local clock: inherit rates unless a shift occurs]
    N --> R[Generate or read tree]
    O --> R
    P --> R
    Q --> R
    R --> S[Apply patient-age and sampling-time constraints]
    S --> T[Simulate CNA events along branches]
    T --> U[Optionally add copy-number calling error]
    U --> V[Write outputs]
    V --> W[Copy-number matrices]
    V --> X[Tree files]
    V --> Y[Mutation and branch info files]
    V --> Z[Relative copy-number files if requested]
    W --> AA([End])
    X --> AA
    Y --> AA
    Z --> AA
```

## Maximum Likelihood Inference workflow (`cnetml`)

`cnetml` reads observed copy-number profiles and either searches for a
maximum-likelihood tree, scores a supplied tree, optimizes a supplied
topology, or reconstructs ancestral states.

```mermaid
flowchart TD
    A([Start cnetml]) --> B[Parse command-line options]
    B --> C[Reject incompatible options, for example bsr_mode with unsupported modes]
    C --> D[Initialize random-number generator]
    D --> E[Read sample times if provided]
    E --> F[Read and preprocess copy-number profiles]
    F --> G{Model}
    G -->|DECOMP model| H[Read copy-number changes by chromosome]
    H --> I[Build observed-change vectors]
    I --> J[Estimate dimensions needed for WGD, chromosome, and segment chains]
    G -->|Other models| K[Read copy-number states by chromosome]
    K --> L[Build observed-state vectors]
    J --> M[Build likelihood configuration]
    L --> M
    M --> N[Build optimization configuration]
    N --> O{mode}
    O -->|0: infer ML tree| P[Initialize tree-search state]
    P --> Q{tree_search}
    Q -->|0| R[Evolutionary algorithm]
    Q -->|1| S[Random-restart hill climbing with NNI]
    Q -->|2| T[Exhaustive topology search]
    R --> U[Optimize each candidate tree]
    S --> U
    T --> U
    U --> V[Keep best-scoring tree]
    O -->|2: score given tree| W[Load tree]
    W --> X[Compute likelihood once]
    O -->|3: optimize given tree| Y[Load tree]
    Y --> Z[Optimize branch lengths and rates]
    O -->|4: ancestral states| AB[Load tree]
    AB --> AC[Infer marginal and/or joint ancestral states]
    O -->|5: segment file only| AD[Write postprocessed segment file]
    V --> AE[Write tree, summary, Nexus, and edge-rate reports]
    X --> AF[Print log likelihood]
    Z --> AE
    AC --> AG[Write ancestral-state output]
    AD --> AH([End])
    AE --> AH
    AF --> AH
    AG --> AH
```

For how the DECOMP likelihood itself is computed once this workflow reaches
"Build likelihood configuration", see
[Architecture](../developer-guide/architecture.md).

## MCMC workflow (`cnetmcmc`)

:::{warning}
**Draft.** `cnetmcmc` is not yet officially released — current focus is
`cnets` and `cnetml`. This diagram was traced from `code/cnetmcmc.cpp`
(`main`, `run_mcmc`, `run_with_reference_tree`) rather than from a
reviewed design doc, and has not been checked against the running
program's actual behaviour. Treat it as a starting point, not a
reference.
:::

```mermaid
flowchart TD
    A([Start cnetmcmc]) --> B[Parse command-line options / mcmc.cfg]
    B --> C[Set up RNG with seed]
    C --> D[Read copy-number input, build observation vectors]
    D --> E[Read optional sample-timing file]
    E --> F{Reference tree given? --rtreefile}
    F -->|Yes| G[Load reference tree, compute its likelihood]
    G --> H{fix_topology}
    H -->|1| I[Start tree: random branch lengths on reference topology]
    H -->|0| J[Start tree: random coalescent tree]
    F -->|No| K{init_tree}
    K -->|0| L[Random coalescent tree]
    K -->|1| M[Provided tree, --file_itree]
    K -->|2| N[Random tree with reference topology]
    I --> O[Assign initial mutation rates, compute start-tree likelihood]
    J --> O
    L --> O
    M --> O
    N --> O
    O --> P[run_mcmc: begin chain]
    P --> Q{Next of n_draws iterations}
    Q --> R[Randomly pick a move: topology / branch length / rate]
    R --> S[Propose new state, compute Metropolis-Hastings ratio]
    S --> T[Accept or reject]
    T --> U{i > n_burnin and i % n_gap == 0?}
    U -->|Yes| V[Append sample to *.p trace and *.t tree file]
    U -->|No| Q
    V --> Q
    Q -->|n_draws reached| W([End])
```

:::{note}
Whether `--rtreefile` is given only decides where the *starting* tree
comes from. Whether the topology is actually held fixed during
sampling is a separate flag, `fix_topology`, used on both paths — this
differs from the summary in the top-level `README.md`, which describes
the reference tree itself as fixing the topology.
:::

## Which tool produces what

| file(s) | produced by | consumed by | status |
| --- | --- | --- | --- |
| `*-cn.txt.gz` / `*-haplotype-cn.txt.gz` | `cnets` (simulated); or preprocessing scripts for real data | `cnetml`, `cnetmcmc` | required |
| `*-rcn.txt.gz` / `*-haplotype-rcn.txt.gz` | `cnets` | external tools expecting relative copy number | optional |
| `*-inodes-cn.txt.gz` / `*-inodes-haplotype-cn.txt.gz` | `cnets` | downstream analysis of internal-node truth | diagnostic |
| `*-rel-times.txt` | `cnets` (simulated); or supplied directly for real data | `cnetml`, `cnetmcmc` | optional (needed for `estmu`/mutation-rate estimation) |
| `<prefix>-tree.txt`, `<prefix>-tree.nex`, `<prefix>-tree-nmut.nex` | `cnets` (ground truth) | `cnetml` (as an optional initial tree), `cnetmcmc` (`--file_itree`), accuracy comparisons | primary output |
| `<ofile>`, `<ofile>.nex`, `<ofile>.nmut.nex` | `cnetml` (reconstructed; `<ofile>` set by `-o`/`--ofile`, **not** `*-tree.*` — see [File formats](../file-formats/index.md#tree-formats)) | tree viewers (e.g. FigTree), downstream analysis, accuracy comparisons against `cnets`' output | primary output |
| `<ofile>.summary.txt`, `<ofile>.edge_rates.txt` | `cnetml` | run bookkeeping, per-branch rate analysis (`bsr_mode > 0`) | diagnostic |
| `*-info.txt`, `*-mut.txt`, `*-edge_rates.txt` | `cnets` | mutation-mapping / accuracy analysis; `*-edge_rates.txt` against `cnetml`'s reconstructed `<ofile>.edge_rates.txt` | diagnostic |
| `*-segs.txt` | `cnetml` (written when reading the input copy-number file, named by `--seg_file`) | intermediate to `cnetml`'s own likelihood computation | intermediate |
| `<ofile>.mrca.cn`, `<ofile>.joint.cn`, and (models 0-2) `<ofile>.mrca.state`/`<ofile>.joint.state`, or (model 3) `<ofile>.mrca.{seg,chr,wgd}.state`/`<ofile>.joint.state` | `cnetml` (mode 4, ancestral-state reconstruction) | downstream analysis of ancestral states | primary output (mode 4 only) |
| `mcmc.cfg` | hand-written / copied from the repository root | `cnetmcmc` (`--config_file`) | required by `cnetmcmc` |
| `*.p`, `*.t` | `cnetmcmc` (via `--trace_param_file`/`--trace_tree_file`; C++ defaults are `trace-mcmc-params.txt`/`trace-mcmc-trees.txt` — `run-cnetmcmc.sh` is what gives them `.p`/`.t` extensions) | Tracer / RWTY (`*.p`), TreeAnnotator (`*.t`) | primary output |

See [File formats](../file-formats/index.md) for column definitions of
each file.

## Future direction

TODO: in-memory exchange via `libcneta`; Nextflow as an orchestration layer.

A Nextflow pipeline from raw sequencing data, or sequence alignments, or copy number calls to phylogenetic analysis and downstream analysis (such as copy number signature attachment to the tree branches).
<!-- https://github.com/sivaranjanjohnson/building  -->
