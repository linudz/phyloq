#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})

usage <- function() {
  paste(
    "Usage:",
    "Rscript build_gene_trees.R [alignment_dir] [reference_tree] [output_dir] [pattern] [taxon_map]",
    "",
    "Defaults:",
    "  alignment_dir  ../fast.run/alignments",
    "  reference_tree inputs/science.abn7829_data_s4.nex.tree",
    "  output_dir     inputs/trees",
    "  pattern        \\.phy$",
    "  taxon_map      inputs/taxon_name_map.tsv",
    sep = "\n"
  )
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(script_arg)) stop("Cannot determine script location.", call. = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) && args[[1]] %in% c("-h", "--help")) {
  cat(usage(), "\n")
  quit(status = 0)
}

alignment_dir <- if (length(args) >= 1L) args[[1]] else file.path(script_dir, "..", "fast.run", "alignments")
reference_file <- if (length(args) >= 2L) args[[2]] else file.path(script_dir, "inputs", "science.abn7829_data_s4.nex.tree")
output_dir <- if (length(args) >= 3L) args[[3]] else file.path(script_dir, "inputs", "trees")
pattern <- if (length(args) >= 4L) args[[4]] else "\\.phy$"
taxon_map_file <- if (length(args) >= 5L) args[[5]] else file.path(script_dir, "inputs", "taxon_name_map.tsv")

alignment_dir <- normalizePath(alignment_dir, mustWork = TRUE)
reference_file <- normalizePath(reference_file, mustWork = TRUE)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
taxon_map_file <- normalizePath(taxon_map_file, mustWork = TRUE)

taxon_map <- read.delim(taxon_map_file, stringsAsFactors = FALSE, check.names = FALSE)
required_map_columns <- c("AlignmentName", "ReferenceTreeName")
if (!all(required_map_columns %in% names(taxon_map))) {
  stop("Taxon map must contain: ", paste(required_map_columns, collapse = ", "), call. = FALSE)
}
if (anyDuplicated(taxon_map$AlignmentName)) stop("Duplicated AlignmentName entries in taxon map.", call. = FALSE)
name_lookup <- setNames(taxon_map$ReferenceTreeName, taxon_map$AlignmentName)

read_relaxed_phylip_aa <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (!length(lines)) stop("Empty alignment: ", path, call. = FALSE)

  header <- strsplit(trimws(lines[[1]]), "[[:space:]]+")[[1]]
  if (length(header) < 2L || anyNA(suppressWarnings(as.integer(header[1:2])))) {
    stop("Invalid PHYLIP header in ", path, call. = FALSE)
  }
  expected_taxa <- as.integer(header[[1]])
  expected_sites <- as.integer(header[[2]])

  records <- lines[-1L]
  records <- records[nzchar(trimws(records))]
  parsed <- regexec("^[[:space:]]*([^[:space:]]+)[[:space:]]+(.+?)[[:space:]]*$", records)
  fields <- regmatches(records, parsed)
  if (any(lengths(fields) != 3L)) {
    stop("Cannot parse one or more PHYLIP records in ", path, call. = FALSE)
  }

  taxa <- vapply(fields, `[[`, character(1), 2L)
  sequences <- vapply(fields, `[[`, character(1), 3L)
  sequences <- gsub("[[:space:]]", "", sequences)
  sequences <- toupper(sequences)

  if (length(taxa) != expected_taxa) {
    stop("Expected ", expected_taxa, " taxa but found ", length(taxa), " in ", path, call. = FALSE)
  }
  if (anyDuplicated(taxa)) stop("Duplicated taxon labels in ", path, call. = FALSE)
  lengths_seq <- nchar(sequences)
  if (any(lengths_seq != expected_sites)) {
    bad <- taxa[lengths_seq != expected_sites]
    stop("Sequence-length mismatch in ", path, ": ", paste(bad, collapse = ", "), call. = FALSE)
  }

  allowed <- strsplit("ARNDCQEGHILKMFPSTWYVBZX*-?", "", fixed = TRUE)[[1]]
  observed <- unique(unlist(strsplit(sequences, "", fixed = TRUE), use.names = FALSE))
  invalid <- setdiff(observed, allowed)
  if (length(invalid)) {
    stop("Unsupported amino-acid symbols in ", path, ": ", paste(invalid, collapse = ", "), call. = FALSE)
  }

  matrix_data <- do.call(rbind, strsplit(sequences, "", fixed = TRUE))
  rownames(matrix_data) <- taxa
  phangorn::phyDat(matrix_data, type = "AA")
}

read_reference_tree <- function(path) {
  tree <- tryCatch(ape::read.tree(path), error = function(e) NULL)
  if (is.null(tree)) tree <- ape::read.nexus(path)
  if (inherits(tree, "multiPhylo")) tree <- tree[[1L]]
  if (!inherits(tree, "phylo")) stop("Reference tree could not be read.", call. = FALSE)
  if (anyDuplicated(tree$tip.label)) stop("Duplicated labels in reference tree.", call. = FALSE)
  tree
}

fit_gene_tree <- function(alignment_file, reference_tree, output_dir, name_lookup) {
  alignment <- read_relaxed_phylip_aa(alignment_file)
  taxa <- names(alignment)
  mapped <- unname(name_lookup[taxa])
  taxa[!is.na(mapped)] <- mapped[!is.na(mapped)]
  if (anyDuplicated(taxa)) {
    stop("Taxon-name reconciliation created duplicate labels in ", basename(alignment_file), call. = FALSE)
  }
  names(alignment) <- taxa
  missing_from_tree <- setdiff(taxa, reference_tree$tip.label)
  if (length(missing_from_tree)) {
    stop(
      "Alignment taxa absent from reference tree in ", basename(alignment_file), ": ",
      paste(missing_from_tree, collapse = ", "), call. = FALSE
    )
  }

  gene_tree <- ape::keep.tip(reference_tree, taxa)
  alignment <- alignment[gene_tree$tip.label]

  initial <- phangorn::pml(gene_tree, alignment, model = "LG")
  fitted <- phangorn::optim.pml(
    initial,
    model = "LG",
    optEdge = TRUE,
    optNni = FALSE,
    optBf = FALSE,
    optQ = FALSE,
    optInv = FALSE,
    optGamma = FALSE,
    control = phangorn::pml.control(trace = 0)
  )

  result <- fitted$tree
  if (any(!is.finite(result$edge.length)) || any(result$edge.length < 0)) {
    stop("Invalid fitted branch lengths for ", basename(alignment_file), call. = FALSE)
  }
  result$edge.length[result$edge.length == 0] <- 1e-08

  output_file <- file.path(output_dir, paste0(basename(alignment_file), ".nw"))
  ape::write.tree(result, file = output_file, digits = 10)

  data.frame(
    Alignment = basename(alignment_file),
    Tree = basename(output_file),
    Taxa = length(result$tip.label),
    Sites = sum(attr(alignment, "weight")),
    Model = "LG",
    LogLikelihood = fitted$logLik,
    stringsAsFactors = FALSE
  )
}

alignment_files <- sort(list.files(alignment_dir, pattern = pattern, full.names = TRUE))
if (!length(alignment_files)) {
  stop("No alignments matched pattern '", pattern, "' in ", alignment_dir, call. = FALSE)
}

reference_tree <- read_reference_tree(reference_file)
cat("Reference tree:", reference_file, "\n")
cat("Alignment directory:", alignment_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Taxon map:", taxon_map_file, "\n")
cat("Alignments:", length(alignment_files), "\n\n")

summary_rows <- vector("list", length(alignment_files))
for (i in seq_along(alignment_files)) {
  cat(sprintf("[%d/%d] %s\n", i, length(alignment_files), basename(alignment_files[[i]])))
  summary_rows[[i]] <- fit_gene_tree(alignment_files[[i]], reference_tree, output_dir, name_lookup)
}

summary_table <- do.call(rbind, summary_rows)
summary_file <- file.path(output_dir, "gene_tree_build_summary.tsv")
write.table(summary_table, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)

tree_lines <- vapply(
  file.path(output_dir, summary_table$Tree),
  function(path) paste(readLines(path, warn = FALSE), collapse = ""),
  character(1)
)
rer_input <- data.frame(Gene = summary_table$Alignment, Tree = tree_lines, stringsAsFactors = FALSE)
rer_input_file <- file.path(output_dir, "rerconverge_gene_trees.tsv")
write.table(rer_input, rer_input_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("\nCreated", nrow(summary_table), "gene trees.\n")
cat("Summary:", summary_file, "\n")
cat("RERconverge input:", rer_input_file, "\n")
