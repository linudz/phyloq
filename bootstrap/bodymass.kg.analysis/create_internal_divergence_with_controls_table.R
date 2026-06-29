#!/usr/bin/env Rscript

# Build an expanded high/low phenotype table for families that contain at least
# one within-family divergent contrast. All within-family species pairs from
# those families are retained and labelled as DIVERGENT when Nu > 0.95, otherwise
# CONTROL.

pairwise_file <- "ranked_pairwise_results.with_families.tsv"
phenotype_file <- file.path("BodyMass_kg", "BodyMass_kg.BM.single_values.tsv")
output_file <- "bodymasskg_internal_family_divergence_with_controls.tsv"

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

within_family <- pairwise[pairwise$Family1 == pairwise$Family2, ]
divergent_families <- sort(unique(within_family$Family1[within_family$Nu > 0.95]))

selected <- within_family[within_family$Family1 %in% divergent_families, ]
phenotype_lookup <- setNames(phenotypes$BodyMass_kg, phenotypes$SpeciesBROAD)

rows <- lapply(seq_len(nrow(selected)), function(i) {
  row <- selected[i, ]
  value1 <- unname(phenotype_lookup[row$Species1])
  value2 <- unname(phenotype_lookup[row$Species2])

  if (is.na(value1) || is.na(value2)) {
    stop("Missing phenotype value for row ", i)
  }

  if (value1 >= value2) {
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
    contrast_type = ifelse(row$Nu > 0.95, "DIVERGENT", "CONTROL"),
    specie_con_valore_alto = species_high,
    specie_con_valore_basso = species_low,
    valore_specie_con_valore_alto = value_high,
    valore_specie_con_valore_basso = value_low,
    nu = row$Nu,
    stringsAsFactors = FALSE
  )
})

out <- do.call(rbind, rows)
out <- out[order(
  out$nome_famiglia,
  out$contrast_type,
  -out$nu,
  out$specie_con_valore_alto,
  out$specie_con_valore_basso
), ]

write.table(
  out,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote: ", output_file)
message("Families: ", paste(divergent_families, collapse = ", "))
message("Rows: ", nrow(out))
message("DIVERGENT: ", sum(out$contrast_type == "DIVERGENT"))
message("CONTROL: ", sum(out$contrast_type == "CONTROL"))
