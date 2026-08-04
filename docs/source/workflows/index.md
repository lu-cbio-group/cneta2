# Workflows and data flow

How the tools are chained together, and what data passes between them.

:::{note}
The diagrams below cover `cnets` and `cnetml`. `cnetmcmc` is less mature and
intentionally not diagrammed yet — see [MCMC workflow](#mcmc-workflow-cnetmcmc).
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

TODO + diagram. `cnetmcmc` is under active development; no diagram has been
supplied for it yet.

## Which tool produces what

TODO: table mapping outputs to the tools that consume them, and marking
each as required, optional, diagnostic, or intermediate.

## Future direction

TODO: in-memory exchange via `libcneta`; Nextflow as an orchestration layer.

A Nextflow pipeline from raw sequencing data, or sequence alignments, or copy number calls to phylogenetic analysis and downstream analysis (such as copy number signature attachment to the tree branches).
<!-- https://github.com/sivaranjanjohnson/building  -->
