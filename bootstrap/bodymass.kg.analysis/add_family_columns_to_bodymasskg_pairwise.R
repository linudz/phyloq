#!/usr/bin/env Rscript

# Add family annotations to BodyMass_kg PhyloQ pairwise output.

base_dir <- "."
taxonomy_file <- "../250612.trait.bestmodel.fig1/inputs/nhp.phenomic.dataset.tsv"

extract_family <- function(group_name) {
  parts <- strsplit(group_name, "_", fixed = TRUE)[[1]]
  if (length(parts) < 2) {
    return(NA_character_)
  }
  parts[[2]]
}

taxonomy <- read.delim(
  taxonomy_file,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

species_to_family <- setNames(
  vapply(taxonomy$GroupName, extract_family, character(1)),
  taxonomy$SpeciesBROAD
)

add_families <- function(input_file, output_file) {
  tab <- read.delim(
    input_file,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  family1 <- unname(species_to_family[tab$Species1])
  family2 <- unname(species_to_family[tab$Species2])

  missing_species <- sort(unique(c(
    tab$Species1[is.na(family1)],
    tab$Species2[is.na(family2)]
  )))
  if (length(missing_species) > 0) {
    warning("Missing family for: ", paste(missing_species, collapse = ", "))
  }

  species2_idx <- match("Species2", names(tab))
  out <- cbind(
    tab[seq_len(species2_idx)],
    Family1 = family1,
    Family2 = family2,
    tab[(species2_idx + 1):ncol(tab)]
  )

  write.table(
    out,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

add_families(
  file.path(base_dir, "pairwise_results.tsv"),
  file.path(base_dir, "pairwise_results.with_families.tsv")
)
add_families(
  file.path(base_dir, "ranked_pairwise_results.tsv"),
  file.path(base_dir, "ranked_pairwise_results.with_families.tsv")
)

message("Wrote family-annotated BodyMass_kg pairwise outputs.")
