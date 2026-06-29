#!/usr/bin/env Rscript

# Build a high/low phenotype table for within-family BodyMass_kg divergence
# contrasts. A contrast is retained when both species are in the same family and
# Nu > 0.95. Within each retained pair, species are ordered by observed
# BodyMass_kg value.

pairwise_file <- "ranked_pairwise_results.with_families.tsv"
phenotype_file <- file.path("BodyMass_kg", "BodyMass_kg.BM.single_values.tsv")
output_file <- "bodymasskg_internal_family_divergence_high_low.tsv"

pairwise <- read.delim(
  pairwise_file,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

phenotypes <- read.delim(
  phenotype_file,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_pairwise <- c("Species1", "Species2", "Family1", "Family2", "Nu")
missing_pairwise <- setdiff(required_pairwise, names(pairwise))
if (length(missing_pairwise) > 0) {
  stop("Missing pairwise column(s): ", paste(missing_pairwise, collapse = ", "))
}

required_phenotypes <- c("SpeciesBROAD", "BodyMass_kg")
missing_phenotypes <- setdiff(required_phenotypes, names(phenotypes))
if (length(missing_phenotypes) > 0) {
  stop("Missing phenotype column(s): ", paste(missing_phenotypes, collapse = ", "))
}

selected <- pairwise[pairwise$Family1 == pairwise$Family2 & pairwise$Nu > 0.95, ]
phenotype_lookup <- setNames(phenotypes$BodyMass_kg, phenotypes$SpeciesBROAD)

rows <- lapply(seq_len(nrow(selected)), function(i) {
  row <- selected[i, ]
  value1 <- unname(phenotype_lookup[row$Species1])
  value2 <- unname(phenotype_lookup[row$Species2])

  if (is.na(value1) || is.na(value2)) {
    stop("Missing phenotype value for row ", i)
  }
  if (value1 == value2) {
    stop("Tied phenotype values for row ", i, ": ", row$Species1, " and ", row$Species2)
  }

  if (value1 > value2) {
    species_high <- row$Species1
    species_low <- row$Species2
    value_high <- value1
    value_low <- value2
  } else {
    species_high <- row$Species2
    species_low <- row$Species1
    value_high <- value2
    value_low <- value1
  }

  data.frame(
    nome_famiglia = row$Family1,
    specie_con_valore_alto = species_high,
    specie_con_valore_basso = species_low,
    valore_specie_con_valore_alto = value_high,
    valore_specie_con_valore_basso = value_low,
    nu = row$Nu,
    stringsAsFactors = FALSE
  )
})

out <- do.call(rbind, rows)
out <- out[order(out$nome_famiglia, -out$nu, out$specie_con_valore_alto, out$specie_con_valore_basso), ]

write.table(
  out,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote: ", output_file)
message("Rows: ", nrow(out))
