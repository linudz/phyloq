# Versioned workflow inputs

## Observed states

`states/*.observed.state.rds` are the frozen outputs of workflow 01 in the
local mammalian brain/body case study. Each state contains the matched
1,427-species phylogeny, trait values, taxonomy, observed BM/OU fits, selected
model, observed opportunity depths and tail settings required by workflow 02.

Traits:

- `log_body_mass`
- `log_brain_mass`
- `relative_brain_mass`

## PSS implementation

`core/pss.core.R` is the exact PSS implementation used for the observed run.
It is copied into this analysis to prevent cluster execution from depending on
paths outside the `phyloq` repository.

Run `shasum -a 256 -c MANIFEST.sha256` from the workflow directory to verify
these files and the repository-level Conda environment.
