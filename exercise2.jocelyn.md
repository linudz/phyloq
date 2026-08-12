# Exercise 2: Install RERconverge and Test the New Bootstrap Pipeline

## Goal of the Exercise

The goal of this exercise is to install `RERconverge` and run the new test
version of the PhyloQ pipeline.

This is a technical exercise. You do not need to interpret the biological
results yet.

The exercise is successful if:

1. `RERconverge` is installed correctly.
2. The Nextflow pipeline starts.
3. The CAAStools and RERconverge tasks finish without errors.
4. The expected output folders are created.

## 1. Update the Repository

Open a terminal and move into your local `phyloq` repository:

```bash
cd /path/to/phyloq
```

Replace `/path/to/phyloq` with the real path on your computer.

Update the repository:

```bash
git pull
```

This is important because the exercise uses the new RERconverge scripts,
configuration file, and Nextflow workflow.

## 2. Make Conda Available

The repository contains an `environment.yml` file with the software needed for this exercise, including Python, R, Java, Nextflow, compilers, and the R dependencies.

On a computing cluster, Conda is often provided as a software module. Your cluster is different from Fabio's cluster, so the exact module name may also be different. Search for the available module:

```bash
module spider conda
module spider anaconda
module spider miniconda
```

The command opens a scrollable help page. Press `q` to leave it. Load the Conda or Miniconda module shown by your cluster. For example:

```bash
module load Miniconda3
```

This is only an example. Use the exact module name and version shown on your server. If you cannot find one, check your cluster documentation or ask its support team.

Confirm that Conda is available:

```bash
conda --version
```

## 3. Create the PhyloQ Environment

Run this command from the root of the `phyloq` repository:

```bash
conda env create -f environment.yml
```

Conda may need several minutes to solve and install the environment. If an environment called `phyloq` already exists, update it instead:

```bash
conda env update -n phyloq -f environment.yml --prune
```

## 4. Activate and Check the Environment

Activate it:

```bash
conda activate phyloq
```

If Conda says `Run 'conda init' before 'conda activate'`, run:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phyloq
```

Then check the main programs:

```bash
which python
which Rscript
python --version
Rscript --version
nextflow -version
python -c "import Bio, numpy, scipy, pkg_resources; print('Python dependencies: OK')"
```

The `python` and `Rscript` paths should point inside an environment named `phyloq`. A deprecation warning from `pkg_resources` is acceptable; the final `Python dependencies: OK` message should still appear.

## 5. Read the README and Check the Inputs

Read the main `README.md`, especially the sections about the Conda environment, RERconverge, running the pipeline, and troubleshooting.

RERconverge does not need to be installed by hand before the first test. Nextflow checks it once and installs it automatically from GitHub if it is missing. If that automatic step fails because compute nodes cannot access GitHub, follow the manual installation section in the README from the login node.

Move into the bootstrap directory:

```bash
cd bootstrap
```

Check the input files:

```bash
test -s rerconverge/inputs/trees/rerconverge_gene_trees.tsv
test -s rerconverge/inputs/science.abn7829_data_s4.nex.tree
```

No output means that both files exist and are not empty.

## 6. Run the New Bootstrap Pipeline

Start the pipeline with the default test settings:

```bash
bash run_pipeline.sh
```

The launcher prints the new run name, for example `run-2608121236`. It creates that directory under `results/` in the repository root and stores the name in `.last_run_id` for the checks below and for `-resume`.

The two parts of the pipeline run in parallel:

- CAAStools analyses five test alignments with three phenotype
  configurations;
- RERconverge performs one analysis for each of the three phenotype
  configurations, using three test genes together.

Nextflow first checks RERconverge once. If the package is missing, the setup process tries to install it automatically from GitHub before starting the RERconverge analyses. CAAStools can start while this check is running. If the cluster compute node cannot access GitHub, use the manual installation steps in the README from the login node and then resume the pipeline.

The RERconverge analysis must use several genes at the same time. It is not
launched separately for each gene.

Wait until Nextflow finishes. A successful run should end without a red error
message and should show the processes as completed.

## 7. Check the Output

All outputs are outside the input folders. They are stored under:

```text
results/run-YYMMDDHHMM/
├── caas.results/
└── rer.results/
```

From `bootstrap/`, `$(cat .last_run_id)` prints the current run directory name.
To print the complete path to the current results:

```bash
echo ../results/"$(cat .last_run_id)"
```

The `results/` directory is local and is not uploaded to GitHub. A new run
gets a new timestamped directory; `-resume` continues in the same directory.

### CAAStools output

Count the CAAStools result files:

```bash
find ../results/"$(cat .last_run_id)"/caas.results -type f -name '*.caas' | wc -l
```

The expected number is:

```text
15
```

This corresponds to:

```text
5 alignments x 3 phenotype configurations = 15 CAAStools tasks
```

Each `.caas` file is one alignment tested with one phenotype configuration.
The `.cfg` file divides the species into groups, and CAAStools reports
amino-acid positions whose pattern matches that grouping. The filename names
both the gene/alignment and the configuration. An empty file means that the
task finished but found no candidate position; it is not automatically an
error.

List all CAAStools outputs:

```bash
find ../results/"$(cat .last_run_id)"/caas.results -type f -name '*.caas' | sort
```

### RERconverge output

List the RERconverge result directories:

```bash
find ../results/"$(cat .last_run_id)"/rer.results -maxdepth 1 -type d -name 'rerconverge.*' | sort
```

You should see three directories:

```text
rerconverge.ctrl.group.001
rerconverge.div.group.001
rerconverge.div.group.002
```

Check that each analysis produced an association table:

```bash
find ../results/"$(cat .last_run_id)"/rer.results -type f -name 'associations.tsv' | wc -l
```

The expected number is:

```text
3
```

Each `rerconverge.<configuration>/` directory is one phenotype configuration
analysed across the three test genes. In its main file, `associations.tsv`,
each row represents one gene:

- `Rho`: direction and strength of the association with foreground branches.
  Positive means relatively faster evolution in the foreground; negative means
  relatively slower evolution.
- `N`: number of branch observations used.
- `P`: unadjusted association p-value.
- `p.adj`: adjusted p-value produced by RERconverge.

The directory also contains `rer_matrix.tsv` (relative evolutionary rates),
`phenotype_paths.tsv` (phenotype mapped onto branches), `run_summary.tsv`
(inputs and thresholds), and `rerconverge_objects.rds` (saved R objects).
`missing_foreground_species.txt` may appear if foreground species are missing
from the tree.

To inspect one association table:

```bash
column -t -s $'\t' ../results/"$(cat .last_run_id)"/rer.results/rerconverge.ctrl.group.001/associations.tsv | less -S
```

Press `q` to leave `less`.

It is normal for some values inside these test tables to be `NA`. This is a
very small three-gene test. In this exercise we are checking that the software
and workflow run correctly, not whether the results are biologically
meaningful.

## 8. When Is the Exercise Complete?

The exercise is complete if all of the following are true:

- the `INSTALL_RERCONVERGE` process and all RERconverge analyses finish successfully;
- Nextflow finishes without errors;
- there are 15 `.caas` files;
- there are 3 RERconverge `associations.tsv` files.

If all four checks pass, the new bootstrap pipeline works correctly on your
computer.

## If Something Goes Wrong

Do not delete `bootstrap/work/` immediately. It is Nextflow's active cache and
contains useful information about failed tasks. A `work/` directory in the
repository root is obsolete and is not used by the current launcher.

First read the troubleshooting sections in `README.md`. If the problem is not
resolved, save and report:

1. the exact command you ran;
2. the complete error message;
3. the output of `R --version`;
4. the output of `nextflow -version`;
5. the last part of `.nextflow.log`.

After fixing the problem, you can continue the same run with:

```bash
bash run_pipeline.sh -resume
```
