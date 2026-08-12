# PhyloQ

PhyloQ is a comparative workflow for testing whether molecular evolutionary
rates and amino-acid substitutions are associated with predefined
phylogenetic phenotype contrasts. The pipeline runs two independent analysis
branches:

1. **CAAStools discovery** searches protein alignments for amino-acid states
   associated with foreground/background species groups.
2. **RERconverge** tests whether relative evolutionary rates across genes are
   associated with the same binary phenotype groups.

The two branches start in parallel. RERconverge does not wait for CAAStools,
and CAAStools does not consume RERconverge output.

## Repository layout

```text
phyloq/
├── README.md
├── bootstrap/
│   ├── launch_pipeline.nf          Nextflow workflow
│   ├── nextflow.config             default pipeline parameters
│   ├── caastools/                  bundled CAAStools checkout
│   ├── alignments/                 complete alignment collection
│   ├── cfg.files/                  complete phenotype configurations
│   ├── fast.run/
│   │   ├── alignments/             five test alignments
│   │   └── configs/                three test configurations
│   └── rerconverge/
│       ├── build_gene_trees.R      constructs gene trees from alignments
│       ├── run_rerconverge.R       RERconverge command-line wrapper
│       └── inputs/
│           ├── science.abn7829_data_s4.nex.tree
│           ├── taxon_name_map.tsv
│           └── trees/
│               └── rerconverge_gene_trees.tsv
├── results/                       timestamped run outputs (Git-ignored)
└── supplementary/
```

Nextflow's active cache is `bootstrap/work/`, because the launcher always runs
the workflow from `bootstrap/`. The matching `bootstrap/.nextflow/`,
`bootstrap/.nextflow.log*`, and `bootstrap/.nextflow/history` paths contain
execution metadata. Keep `bootstrap/work/` and `bootstrap/.nextflow/` if an
interrupted run may need to be resumed. A `work/` directory in the repository
root is left over from an old launch method and is not used by the current
launcher.

## Requirements

The workflow requires:

- Java compatible with the installed Nextflow version;
- Nextflow with DSL2 support;
- Python 3;
- Biopython, NumPy, and SciPy for CAAStools;
- R;
- the R packages `ape`, `phangorn`, and `RERconverge`;
- a C/C++ and Fortran toolchain when RERconverge must be compiled from source.

Check the main programs:

```bash
java -version
nextflow -version
python3 --version
R --version
```

Check the Python libraries:

```bash
python3 -c "import Bio, numpy, scipy; print('Python dependencies: OK')"
```

## Conda environment on a cluster

The repository provides `environment.yml`, which installs Python, R, Java, Nextflow, compilers, and the available RERconverge dependencies in one isolated environment.

Many clusters provide Conda through their module system. Module names differ between institutions, so first search the software catalogue:

```bash
module spider conda
module spider anaconda
module spider miniconda
```

Read the result and load the module available on that cluster, for example:

```bash
module load Miniconda3
```

`Miniconda3` is only an example; use the exact module name and version shown by the cluster. If none of the searches returns a module, consult the local cluster documentation or support team.

From the repository root, create the environment:

```bash
conda env create -f environment.yml
```

If the `phyloq` environment already exists, update it instead:

```bash
conda env update -n phyloq -f environment.yml --prune
```

Activate it:

```bash
conda activate phyloq
```

If Conda reports `Run 'conda init' before 'conda activate'`, initialize Conda for the current shell without changing the account configuration:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phyloq
```

Confirm that the executables come from the environment:

```bash
which python
which Rscript
python --version
Rscript --version
nextflow -version
```

CAAStools currently imports the legacy `pkg_resources` module, so `environment.yml` pins `setuptools<82`. RERconverge itself is checked once by the pipeline and installed automatically from GitHub if it is missing.

## Manual installation of RERconverge

RERconverge is installed from its official GitHub repository rather than
from CRAN. The commands below are intended to be run manually in a terminal.

Official repository:

<https://github.com/nclark-lab/RERconverge>

### 1. Install system compilers

RERconverge contains compiled C++ code and links against a Fortran runtime.
The compiler must therefore be installed before the R package.

On macOS, first install the Xcode command-line tools:

```bash
xcode-select --install
```

Then install a Fortran compiler. The most direct option for an R installation
downloaded from CRAN is the matching GNU Fortran installer provided at:

<https://mac.r-project.org/tools/>

It normally installs the libraries under `/opt/gfortran`, which is where the
official macOS R build expects them.

Alternatively, on an Apple Silicon Mac with Homebrew:

```bash
brew install gcc
which gfortran
```

If `gfortran` is found under `/opt/homebrew` but R still searches
`/opt/gfortran`, create or edit `~/.R/Makevars`:

```bash
mkdir -p "$HOME/.R"
nano "$HOME/.R/Makevars"
```

Add this line on Apple Silicon:

```make
FLIBS = -L/opt/homebrew/opt/gcc/lib/gcc/current -lgfortran -lquadmath
```

Do not replace an existing `Makevars` blindly; add the setting while
preserving any compiler options already present. On Intel Macs, obtain the
Homebrew prefix with `brew --prefix gcc` and adapt the path accordingly.

On Debian or Ubuntu, install the build tools with:

```bash
sudo apt update
sudo apt install build-essential gfortran libcurl4-openssl-dev libssl-dev libxml2-dev
```

### 2. Install the R installation helpers

Start R or run:

```bash
Rscript -e 'install.packages(c("remotes", "BiocManager"), repos="https://cloud.r-project.org")'
```

### 3. Install the Bioconductor dependency

`impute` is required by RERconverge but may not be installed automatically by
the GitHub installer:

```bash
Rscript -e 'BiocManager::install("impute", ask=FALSE, update=FALSE)'
```

### 4. Install RERconverge from GitHub

```bash
Rscript -e 'remotes::install_github("nclark-lab/RERconverge", ref="master", dependencies=TRUE, upgrade="never")'
```

Warnings about package documentation, long paths in bundled vignettes, or
dependencies built with a nearby R patch release do not necessarily indicate
failure. A successful installation ends with:

```text
* DONE (RERconverge)
```

### 5. Verify the installation

```bash
Rscript -e 'library(RERconverge); cat("RERconverge", as.character(packageVersion("RERconverge")), "\n")'
```

For a more complete check:

```bash
Rscript -e 'library(RERconverge); stopifnot(exists("readTrees"), exists("getAllResiduals"), requireNamespace("impute", quietly=TRUE)); cat("RERconverge installation: OK\n")'
```

The pipeline was developed and tested with RERconverge `0.3.0`.

### Common installation errors

#### `dependency 'impute' is not available`

Install `impute` through Bioconductor and repeat the GitHub installation:

```bash
Rscript -e 'BiocManager::install("impute", ask=FALSE, update=FALSE)'
```

#### `library 'gfortran' not found`

The Fortran compiler is absent or R is searching the wrong library path.
Install the CRAN-compatible macOS compiler or configure `FLIBS` in
`~/.R/Makevars` as described above, then repeat the installation.

#### A package library is not writable

Install into a personal R library instead of using administrator privileges:

```bash
mkdir -p "$HOME/R/library"
R_LIBS_USER="$HOME/R/library" Rscript -e 'remotes::install_github("nclark-lab/RERconverge", dependencies=TRUE, upgrade="never")'
```

Use the same `R_LIBS_USER` setting when launching Nextflow, or export it in
the shell before running the pipeline.

## Input data

### Alignments

Each file contains one amino-acid multiple-sequence alignment. The default
format is relaxed PHYLIP (`phylip-relaxed`). Species identifiers must be
compatible with the phenotype configuration files and reference tree.

The default validation run reads five alignments from:

```text
bootstrap/fast.run/alignments/
```

The complete collection is under `bootstrap/alignments/`.

### Phenotype configurations

Each `.cfg` file is a headerless, two-column table:

```text
Species_name    0
Species_name    0
Species_name    1
Species_name    1
```

`0` identifies background species and `1` foreground species. The same files
are consumed independently by CAAStools and RERconverge. Conflicting duplicate
assignments for a species are invalid.

The default run reads three files from:

```text
bootstrap/fast.run/configs/
```

### Gene trees for RERconverge

RERconverge requires several gene trees simultaneously. It cannot estimate a
relative evolutionary rate from a single gene. Two genes are the technical
minimum, but small sets such as the three-gene default are only pipeline
tests; a scientific analysis should use tens, hundreds, or more genes.

The pipeline reads the two-column manifest:

```text
bootstrap/rerconverge/inputs/trees/rerconverge_gene_trees.tsv
```

Each row contains a gene identifier and its Newick tree. The trees are built
on the Kuderna reference topology, with gene-specific branch lengths.

To regenerate trees for the five fast-run alignments:

```bash
cd bootstrap
Rscript rerconverge/build_gene_trees.R
```

To build trees from another alignment directory:

```bash
Rscript rerconverge/build_gene_trees.R \
  "/absolute/path/to/alignments" \
  "rerconverge/inputs/science.abn7829_data_s4.nex.tree" \
  "/absolute/path/to/output/trees" \
  '\.phy$' \
  "rerconverge/inputs/taxon_name_map.tsv"
```

The output directory will contain individual `.nw` files,
`gene_tree_build_summary.tsv`, and `rerconverge_gene_trees.tsv`.

## Pipeline configuration

Defaults are defined in `bootstrap/nextflow.config`.

| Parameter | Purpose | Fast-run default |
|---|---|---|
| `caastools` | CAAStools launcher | bundled `caastools/ct` |
| `alignments` | alignment glob | `fast.run/alignments/*` |
| `config_dir` | phenotype configuration directory | `fast.run/configs` |
| `config_pattern` | pattern inside `config_dir` | `*.cfg` |
| `run_id` | timestamped run identifier supplied by the launcher | generated automatically |
| `results_root` | root directory for all published results | repository `results` directory |
| `fmt` | CAAStools alignment format | `phylip-relaxed` |
| `rerconverge_script` | RERconverge wrapper | `rerconverge/run_rerconverge.R` |
| `rerconverge_trees` | gene-tree manifest | fast-run manifest |
| `rerconverge_master` | topology used to validate gene trees | Kuderna tree |
| `rerconverge_max_trees` | first N manifest trees; `0` means all | `3` |
| `rerconverge_min_trees` | minimum trees required by `readTrees` | `1` |
| `rerconverge_min_valid` | minimum valid genes per branch | `1` |
| `rerconverge_min_species` | minimum usable species per gene | `5` |
| `rerconverge_min_foreground` | minimum foreground branches | `2` |

The deliberately permissive RERconverge thresholds are for the three-gene
technical test. They are not recommended thresholds for inference.

## Running the pipeline

### 1. Open the project directory

Starting from the repository root:

```bash
cd bootstrap
```

Running from `bootstrap/` is recommended because Nextflow then discovers
`nextflow.config` automatically and keeps its execution metadata in one place.

### 2. Optional preflight checks

```bash
nextflow -version
python3 -c "import Bio, numpy, scipy"
Rscript -e 'library(RERconverge); library(ape); library(phangorn)'
test -x caastools/ct
test -s rerconverge/inputs/trees/rerconverge_gene_trees.tsv
test -s rerconverge/inputs/science.abn7829_data_s4.nex.tree
```

No output from the three `test` commands means that the checks succeeded.

### 3. Run both branches

```bash
bash run_pipeline.sh
```

The launcher creates a new `results/run-YYMMDDHHMM/` directory and records its ID in the local, Git-ignored file `bootstrap/.last_run_id`. All published outputs from both branches go into that run directory.

At the start of each workflow run, `INSTALL_RERCONVERGE` checks whether the active R library contains RERconverge. If it is missing, this single setup process installs it from the official GitHub repository before any RERconverge analysis starts. The process requires network access to GitHub only when installation is necessary. CAAStools is independent and can run in parallel while this setup is taking place.

The automatic check does not replace a reproducible software environment: `remotes` and the compiled R dependencies must already be available, as provided by `environment.yml`. On clusters whose compute nodes cannot access GitHub, install RERconverge manually on the login node by following the instructions above before launching Nextflow.

The default validation run starts, in parallel:

- 15 CAAStools tasks: 5 alignments × 3 phenotype configurations;
- 3 RERconverge tasks: one for each configuration, each analysing the first
  3 gene trees in the manifest.

The phrase “three RERconverge tasks” does not mean one task per gene. Each
RERconverge task must receive several genes simultaneously and tests one
phenotype configuration.

### Run only RERconverge

```bash
bash run_pipeline.sh -entry RERCONVERGE_FAST
```

This is useful for checking the R installation and gene-tree inputs without
starting CAAStools.

### Resume an interrupted run

```bash
bash run_pipeline.sh -resume
```

The launcher reads `bootstrap/.last_run_id`, so resume publishes into the same timestamped directory rather than creating a new run.

Use the same inputs and parameter values. Do not remove `bootstrap/work/` or
`bootstrap/.nextflow/` before resuming.

### Override parameters on the command line

For example, use every tree in the manifest:

```bash
bash run_pipeline.sh --rerconverge_max_trees 0
```

Use external configurations:

```bash
bash run_pipeline.sh \
  --config_dir "/absolute/path/to/configurations" \
  --config_pattern "**/*.cfg"
```

Use external alignments and configurations:

```bash
bash run_pipeline.sh \
  --alignments "/absolute/path/to/alignments/*.phy" \
  --config_dir "/absolute/path/to/configurations" \
  --rerconverge_trees "/absolute/path/to/rerconverge_gene_trees.tsv"
```

Quote absolute paths and glob expressions, especially when a path contains
spaces.

### Use an additional Nextflow configuration

```bash
bash run_pipeline.sh -c /absolute/path/to/analysis.config
```

`-c` supplies a Nextflow settings file. `--config_dir` instead identifies the
directory containing biological phenotype `.cfg` inputs.

### Suggested settings for a larger RERconverge analysis

The appropriate thresholds depend on gene and taxon coverage. A reasonable
starting command for a manifest containing many genes is:

```bash
bash run_pipeline.sh \
  --rerconverge_max_trees 0 \
  --rerconverge_min_trees 20 \
  --rerconverge_min_valid 20 \
  --rerconverge_min_species 10
```

Inspect missingness and usable branch coverage before treating these values
as final analytical choices.

## Parallel execution

The workflow has two independent branches. Nextflow schedules CAAStools and
RERconverge as soon as their respective inputs are available. Within the
CAAStools branch, alignment/configuration combinations are independent. Within
the RERconverge branch, phenotype configurations are independent, but genes
inside one RERconverge analysis must remain together.

Nextflow uses the local executor by default and schedules tasks up to the
available resource limit. Executor, CPU, memory, and queue settings can be
added in an external Nextflow configuration when running on a cluster.

## Outputs

### Run directory

Each new launch creates `results/run-YYMMDDHHMM/`. A resumed launch reuses the last run ID stored locally in `bootstrap/.last_run_id`; both the state file and the complete `results/` tree are ignored by Git. Input directories are never used as output destinations.

```text
results/
└── run-YYMMDDHHMM/
    ├── caas.results/
    │   └── <gene>_results/
    │       └── <gene>.<configuration>.caas
    └── rer.results/
        └── rerconverge.<configuration>/
            ├── associations.tsv
            ├── phenotype_paths.tsv
            ├── rer_matrix.tsv
            ├── rerconverge_objects.rds
            ├── run_summary.tsv
            └── missing_foreground_species.txt   optional
```

### CAAStools

Each `.caas` file represents one alignment/configuration combination.
CAAStools compares amino-acid states with the species grouping encoded by the
phenotype `.cfg` file and reports candidate positions matching that grouping.
The filename preserves both input names. An empty file means that CAAStools
completed without reporting a hit; it does not by itself indicate failure.

### RERconverge

- `associations.tsv` contains one row per gene. `Rho` gives the direction and
  strength of association with foreground branches (positive is relatively
  faster; negative is relatively slower), `N` is the number of branch
  observations, `P` is unadjusted, and `p.adj` is package-adjusted.
- `rer_matrix.tsv` contains relative evolutionary rates.
- `phenotype_paths.tsv` records the phenotype mapped onto branches.
- `run_summary.tsv` records effective inputs and thresholds.
- `missing_foreground_species.txt` appears when foreground species are absent
  from the operational RERconverge tree.
- `rerconverge_objects.rds` preserves the R objects for downstream inspection.

With only three genes, `NA` associations are expected when taxon or branch
coverage is insufficient. The fast run validates software integration; it is
not the final biological analysis.

The complete `results/` directory is excluded from Git. Preserve or copy any
run needed for later interpretation before cleaning local outputs.

## Troubleshooting pipeline runs

### RERconverge works in R but not in Nextflow

Confirm that Nextflow sees the same R executable and library:

```bash
which Rscript
Rscript -e '.libPaths(); library(RERconverge)'
```

If a personal library is used, export it before launching Nextflow:

```bash
export R_LIBS_USER="$HOME/R/library"
bash run_pipeline.sh
```

### `Only N phenotype species occur in the tree set`

Too few configuration species are present in the operational gene-tree set.
Check taxon names and gene coverage. Lower `rerconverge_min_species` only for
a technical test, not simply to force an inferential run to complete.

### Foreground species are reported missing

Compare spelling and taxonomic synonyms across the `.cfg` file, alignment,
`taxon_name_map.tsv`, and reference tree. Missing foreground taxa reduce the
effective phenotype sample.

### A run fails after changing inputs

Inspect `bootstrap/.nextflow.log` and the failing task's `.command.err` under
`bootstrap/work/`.
After correcting the problem, use `-resume`. Nextflow will reuse tasks whose
inputs and commands have not changed.

## Reproducibility checklist

Record for every analysis:

- repository commit;
- Nextflow, Java, Python, R, and RERconverge versions;
- exact pipeline parameters and external configuration files;
- checksums or immutable versions of alignments, configurations, and trees;
- numbers of alignments, phenotype configurations, and gene trees;
- effective species and foreground counts from `run_summary.tsv`;
- whether the run was resumed from cache.

## Status

This repository accompanies a manuscript in preparation. Interfaces and file
organization may change while the analysis is being finalized.

For questions, contact fabio.barteri@upf.edu.
