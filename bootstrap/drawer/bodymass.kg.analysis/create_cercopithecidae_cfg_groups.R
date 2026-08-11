#!/usr/bin/env Rscript

# Create binary configuration files from Cercopithecidae BodyMass_kg contrasts.
# Each file contains four non-overlapping species pairs: the high-phenotype
# species from each pair is assigned 1 and the low-phenotype species is assigned
# 0. No species can appear more than once within a file.

input_file <- "bodymasskg_internal_family_divergence_with_controls.tsv"
out_root <- "cfg"
set.seed(20260629)

read_contrasts <- function(input_file) {
  tab <- read.delim(
    input_file,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required <- c(
    "nome_famiglia",
    "contrast_type",
    "specie_con_valore_alto",
    "specie_con_valore_basso"
  )
  missing_cols <- setdiff(required, names(tab))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  tab[tab$nome_famiglia == "Cercopithecidae", ]
}

sample_non_overlapping_pairs <- function(tab, n_pairs = 4, max_attempts = 10000) {
  if (nrow(tab) < n_pairs) {
    stop("Not enough pairs to sample ", n_pairs, " pairs.")
  }

  for (attempt in seq_len(max_attempts)) {
    candidate_idx <- sample(seq_len(nrow(tab)), n_pairs)
    candidate <- tab[candidate_idx, ]
    species <- c(candidate$specie_con_valore_alto, candidate$specie_con_valore_basso)

    if (!anyDuplicated(species)) {
      return(candidate)
    }
  }

  stop("Could not sample ", n_pairs, " non-overlapping pairs after ", max_attempts, " attempts.")
}

write_cfg <- function(pairs, output_file) {
  out <- rbind(
    data.frame(
      species = pairs$specie_con_valore_alto,
      value = 1L,
      stringsAsFactors = FALSE
    ),
    data.frame(
      species = pairs$specie_con_valore_basso,
      value = 0L,
      stringsAsFactors = FALSE
    )
  )

  out <- out[order(out$value, out$species), ]

  write.table(
    out,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

write_group_files <- function(tab, contrast_type, output_dir, prefix, n_files = 100) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  contrast_tab <- tab[tab$contrast_type == contrast_type, ]
  if (nrow(contrast_tab) == 0) {
    stop("No rows available for contrast type: ", contrast_type)
  }

  for (i in seq_len(n_files)) {
    pairs <- sample_non_overlapping_pairs(contrast_tab, n_pairs = 4)
    output_file <- file.path(output_dir, sprintf("%s.group.%03d.cfg", prefix, i))
    write_cfg(pairs, output_file)
  }

  invisible(n_files)
}

contrasts <- read_contrasts(input_file)

n_div <- write_group_files(
  tab = contrasts,
  contrast_type = "DIVERGENT",
  output_dir = file.path(out_root, "div"),
  prefix = "div",
  n_files = 100
)

n_ctrl <- write_group_files(
  tab = contrasts,
  contrast_type = "CONTROL",
  output_dir = file.path(out_root, "ctrl"),
  prefix = "ctrl",
  n_files = 100
)

message("Wrote divergent cfg files: ", n_div)
message("Wrote control cfg files: ", n_ctrl)
