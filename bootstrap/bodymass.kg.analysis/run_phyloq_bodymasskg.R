#!/usr/bin/env Rscript

# Run the standalone R PhyloQ implementation on the BodyMass_kg trait.
# This mirrors the previous Body_mass bootstrap run: BM model, 1000 simulations,
# and tab-separated outputs written into this directory.

source("../R.package/phyloq.main.R")

# geiger can read mc.cores from global options. In some desktop/sandbox
# sessions that value can resolve to zero, so force a valid single-core fit.
options(mc.cores = 1)

result <- phyloq(
  data = "../250612.trait.bestmodel.fig1/inputs/nhp.phenomic.dataset.tsv",
  tree = "../250612.trait.bestmodel.fig1/inputs/phylogeny.nw",
  trait = "BodyMass_kg",
  species_col = "SpeciesBROAD",
  model = "BM",
  nsim = 1000,
  outdir = ".",
  write_outputs = TRUE,
  verbose = TRUE
)

saveRDS(result, file = "phyloq_bodymasskg_result.rds")

write.table(
  result$nu,
  file = "pairwise_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  result$ranked_nu,
  file = "ranked_pairwise_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote PhyloQ BodyMass_kg results to: ", normalizePath("."))
