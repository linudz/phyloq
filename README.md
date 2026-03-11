# Introducing PhyloQ, a test to identify significant phenotype shifts in phylogeny.

## About the project
Evolutionary shifts—large phenotypic deviations between species—can reveal convergence and divergence in trait
evolution, yet are rarely modeled explicitly in standard comparative approaches. Stochastic models such as
Brownian Motion (BM), Ornstein–Uhlenbeck (OU), or acceleration–deceleration (AC/DC) describe neutral evolution
but do not pinpoint where or when shifts occur, which is key for accurate evolutionary inference.

We present PhyloQ, a model-aware test that detects shifts by comparing observed interspecific differences to trait-
specific null distributions under the best-fitting evolutionary model. Applied to 77 continuous primate traits—
including life-history (e.g., maximum lifespan), metabolic, and morphological traits (e.g., body mass, brain size)—and
using the most dense phylogeny available to date (Kuderna et al., 2023), we identified the most likely process (BM,
OU, Early Burst) for each trait, simulated null expectations, and detected significant shifts. Compared to existing
methods (Felsenstein, 1985; Garland et al., 1992; Stayton, 2015; Castiglione et al., 2019), PhyloQ flexibly accounts
for trait-specific dynamics, enabling detection of both convergent and divergent patterns.

Our results highlight the influence of evolutionary model choice on inference and show how model-based
simulations improve macroevolutionary resolution.

## About this repository
This repository distributes the code and data associated with **PhyloQ**, accompanying a manuscript currently in preparation.

Specifically, the repository contains:

1. **Source code of PhyloQ**, (./ folder) including the scripts used to run the test and generate null distributions under different evolutionary models.

2. **A convergence test on real data**, using phenotypic traits obtained from the **Primate Genome-Phenome Archive (PGA)**  (./supplementary/convergence.test folder)
   https://pgarchive.github.io

3. **A validation framework**, (./bootstrap folder) including excess tests and **over-representation analyses** used to evaluate the biological signal detected by PhyloQ.

Together, these materials provide the full computational workflow necessary to reproduce the analyses presented in the manuscript.

For any further question, please write to fabio.barteri@upf.edu
