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

## 2. Read the README

Before installing anything, open the main file:

```text
README.md
```

Follow the section called:

```text
Manual installation of RERconverge
```

The README contains the complete instructions for:

- installing the required compiler;
- installing the R packages `remotes` and `BiocManager`;
- installing the Bioconductor package `impute`;
- installing RERconverge from GitHub;
- solving common installation errors.

Please follow those instructions in order. Do not skip the compiler or
`impute` steps.

## 3. Check the RERconverge Installation

After the installation, run this command:

```bash
Rscript -e 'library(RERconverge); stopifnot(exists("readTrees"), exists("getAllResiduals"), requireNamespace("impute", quietly=TRUE)); cat("RERconverge installation: OK\n")'
```

The expected final message is:

```text
RERconverge installation: OK
```

R may also print some warnings about packages built with a slightly different
R version. Warnings are not always errors. The important points are that the
command finishes and the final `OK` message appears.

If the command stops with `Error`, return to the troubleshooting section of
the README before continuing.

## 4. Move into the Bootstrap Directory

From the repository root, run:

```bash
cd bootstrap
```

Run all the remaining commands from inside this directory.

## 5. Check the Main Programs and Input Files

Check that Nextflow, Python, and R are available:

```bash
nextflow -version
python3 --version
R --version
```

Check the main RERconverge input files:

```bash
test -s rerconverge/inputs/trees/rerconverge_gene_trees.tsv
test -s rerconverge/inputs/science.abn7829_data_s4.nex.tree
```

The two `test` commands normally print nothing. No output means that the files
exist and are not empty.

## 6. Run the New Bootstrap Pipeline

Start the pipeline with the default test settings:

```bash
nextflow run launch_pipeline.nf
```

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

### CAAStools output

Count the CAAStools result files:

```bash
find fast.run/results -type f -name '*.caas' | wc -l
```

The expected number is:

```text
15
```

This corresponds to:

```text
5 alignments x 3 phenotype configurations = 15 CAAStools tasks
```

### RERconverge output

List the RERconverge result directories:

```bash
find fast.run/rerconverge_results -maxdepth 1 -type d -name 'rerconverge.*' | sort
```

You should see three directories:

```text
rerconverge.ctrl.group.001
rerconverge.div.group.001
rerconverge.div.group.002
```

Check that each analysis produced an association table:

```bash
find fast.run/rerconverge_results -type f -name 'associations.tsv' | wc -l
```

The expected number is:

```text
3
```

It is normal for some values inside these test tables to be `NA`. This is a
very small three-gene test. In this exercise we are checking that the software
and workflow run correctly, not whether the results are biologically
meaningful.

## 8. When Is the Exercise Complete?

The exercise is complete if all of the following are true:

- the RERconverge installation check prints `RERconverge installation: OK`;
- Nextflow finishes without errors;
- there are 15 `.caas` files;
- there are 3 RERconverge `associations.tsv` files.

If all four checks pass, the new bootstrap pipeline works correctly on your
computer.

## If Something Goes Wrong

Do not delete the `work` directory immediately. It contains useful information
about failed tasks.

First read the troubleshooting sections in `README.md`. If the problem is not
resolved, save and report:

1. the exact command you ran;
2. the complete error message;
3. the output of `R --version`;
4. the output of `nextflow -version`;
5. the last part of `.nextflow.log`.

After fixing the problem, you can continue the same run with:

```bash
nextflow run launch_pipeline.nf -resume
```
