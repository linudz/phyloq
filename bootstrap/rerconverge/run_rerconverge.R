#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(RERconverge)
})

usage <- function(status = 0L) {
  cat(paste0(
    "Usage:\n",
    "  Rscript run_rerconverge.R --trees FILE --master-tree FILE \\\n",
    "    --phenotype FILE --outdir DIR [options]\n\n",
    "Required:\n",
    "  --trees FILE          Two-column RERconverge gene-tree file\n",
    "  --master-tree FILE    Master tree in Newick format\n",
    "  --phenotype FILE      Two-column species/0-or-1 phenotype file\n",
    "  --outdir DIR          Output directory\n\n",
    "Options:\n",
    "  --clade VALUE         terminal, ancestral, or all [terminal]\n",
    "  --transition VALUE    unidirectional or bidirectional [unidirectional]\n",
    "  --weighted BOOL       Weight foreground clades [false]\n",
    "  --transform VALUE     RER transformation [sqrt]\n",
    "  --impute BOOL         Impute missing branch values [true]\n",
    "  --min-trees INT       Trees required by readTrees [1]\n",
    "  --max-trees INT       Maximum trees to read; 0 means all [0]\n",
    "  --min-species INT     Species required per gene [10]\n",
    "  --min-valid INT       Valid genes required per branch [1]\n",
    "  --min-foreground INT  Foreground branches required [2]\n",
    "  --bootstrap BOOL      Bootstrap weighted correlations [false]\n",
    "  --bootn INT           Bootstrap replicates [1000]\n"
  ))
  quit(status = status)
}

parse_args <- function(x) {
  if (length(x) == 0L || any(x %in% c("-h", "--help"))) usage(0L)
  if (length(x) %% 2L != 0L) stop("Every option must have a value.")
  keys <- sub("^--", "", x[seq(1L, length(x), 2L)])
  if (any(!startsWith(x[seq(1L, length(x), 2L)], "--"))) {
    stop("Options must start with --.")
  }
  vals <- x[seq(2L, length(x), 2L)]
  setNames(as.list(vals), gsub("-", "_", keys, fixed = TRUE))
}

as_bool <- function(x, name) {
  value <- tolower(x)
  if (!value %in% c("true", "false")) stop("--", name, " must be true or false.")
  identical(value, "true")
}

as_int <- function(x, name, minimum = 1L) {
  value <- suppressWarnings(as.integer(x))
  if (is.na(value) || value < minimum) stop("--", name, " must be >= ", minimum, ".")
  value
}

opt <- modifyList(list(
  clade = "terminal",
  transition = "unidirectional",
  weighted = "false",
  transform = "sqrt",
  impute = "true",
  min_trees = "1",
  max_trees = "0",
  min_species = "10",
  min_valid = "1",
  min_foreground = "2",
  bootstrap = "false",
  bootn = "1000"
), parse_args(commandArgs(trailingOnly = TRUE)))

required <- c("trees", "master_tree", "phenotype", "outdir")
missing <- required[!required %in% names(opt)]
if (length(missing)) stop("Missing required option(s): ", paste(missing, collapse = ", "))

allowed <- c(required, "clade", "transition", "weighted", "transform", "impute",
             "min_trees", "max_trees", "min_species", "min_valid", "min_foreground",
             "bootstrap", "bootn")
unknown <- setdiff(names(opt), allowed)
if (length(unknown)) stop("Unknown option(s): ", paste(unknown, collapse = ", "))
if (!opt$clade %in% c("terminal", "ancestral", "all")) stop("Invalid --clade.")
if (!opt$transition %in% c("unidirectional", "bidirectional")) stop("Invalid --transition.")

weighted <- as_bool(opt$weighted, "weighted")
impute <- as_bool(opt$impute, "impute")
bootstrap <- as_bool(opt$bootstrap, "bootstrap")
min_trees <- as_int(opt$min_trees, "min-trees")
max_trees <- as_int(opt$max_trees, "max-trees", minimum = 0L)
min_species <- as_int(opt$min_species, "min-species")
min_valid <- as_int(opt$min_valid, "min-valid")
min_foreground <- as_int(opt$min_foreground, "min-foreground")
bootn <- as_int(opt$bootn, "bootn")

for (file in c(opt$trees, opt$master_tree, opt$phenotype)) {
  if (!file.exists(file)) stop("Input file does not exist: ", file)
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

phenotype <- read.table(opt$phenotype, header = FALSE, sep = "", stringsAsFactors = FALSE,
                        col.names = c("species", "state"), comment.char = "#")
if (!nrow(phenotype)) stop("Phenotype file is empty.")
if (any(!phenotype$state %in% c(0, 1))) stop("Phenotype states must be 0 or 1.")
conflicts <- aggregate(state ~ species, phenotype, function(z) length(unique(z)))
if (any(conflicts$state > 1L)) {
  stop("Conflicting duplicated phenotype states for: ",
       paste(conflicts$species[conflicts$state > 1L], collapse = ", "))
}
phenotype <- phenotype[!duplicated(phenotype$species), , drop = FALSE]
foreground <- phenotype$species[phenotype$state == 1]
if (length(foreground) < 2L) stop("At least two foreground species are required.")

master_tree <- read.tree(opt$master_tree)
if (is.null(master_tree)) stop("Could not read the master tree as Newick.")

tree_tokens <- scan(opt$trees, what = "character", quiet = TRUE)
if (length(tree_tokens) %% 2L != 0L) stop("The gene-tree file must contain name/tree pairs.")
tree_text <- tree_tokens[seq(2L, length(tree_tokens), 2L)]
if (max_trees > 0L) tree_text <- head(tree_text, max_trees)
if (length(tree_text) < 2L) {
  stop("RERconverge requires at least two genes to estimate relative evolutionary rates.")
}
selected_gene_trees <- lapply(tree_text, function(z) read.tree(text = z))
selected_species <- unique(unlist(lapply(selected_gene_trees, function(z) z$tip.label)))
missing_master_species <- setdiff(selected_species, master_tree$tip.label)
if (length(missing_master_species)) stop("Gene-tree species absent from master tree: ", paste(missing_master_species, collapse = ", "))
master_tree <- keep.tip(master_tree, selected_species)
max_tree_species <- max(vapply(selected_gene_trees, function(z) length(z$tip.label), integer(1)))
topology_ok <- vapply(selected_gene_trees, function(gene_tree) {
  reference <- keep.tip(master_tree, gene_tree$tip.label)
  isTRUE(all.equal.phylo(unroot(gene_tree), unroot(reference), use.edge.length = FALSE))
}, logical(1))
if (any(!topology_ok)) {
  stop("Gene trees discordant with the supplied master topology: ",
       paste(which(!topology_ok), collapse = ", "))
}

message("Reading gene trees...")
trees_obj <- readTrees(
  file = opt$trees,
  max.read = if (max_trees == 0L) NA else max_trees,
  minTreesAll = min_trees
)

analysis_species <- intersect(phenotype$species, trees_obj$masterTree$tip.label)
if (length(analysis_species) < min_species) stop("Only ", length(analysis_species), " phenotype species occur in the tree set; --min-species is ", min_species, ".")

message("Calculating relative evolutionary rates...")
rer_matrix <- getAllResiduals(
  treesObj = trees_obj,
  transform = opt$transform,
  impute = impute,
  min.sp = min_species,
  min.valid = min_valid,
  useSpecies = analysis_species
)

available <- intersect(foreground, trees_obj$masterTree$tip.label)
missing_foreground <- setdiff(foreground, available)
if (length(available) < 2L) stop("Fewer than two foreground species occur in the analysed tree set.")

message("Mapping binary phenotype to branches...")
phenotype_tree <- foreground2Tree(
  foreground = available,
  treesObj = trees_obj,
  plotTree = FALSE,
  clade = opt$clade,
  weighted = weighted,
  transition = opt$transition,
  useSpecies = analysis_species
)
phenotype_paths <- tree2Paths(phenotype_tree, trees_obj)

message("Testing phenotype/RER association...")
association <- correlateWithBinaryPhenotype(
  RERmat = rer_matrix,
  charP = phenotype_paths,
  min.sp = min_species,
  min.pos = min_foreground,
  weighted = if (weighted) TRUE else FALSE,
  bootstrap = bootstrap,
  bootn = bootn,
  sort = TRUE
)

write.table(as.data.frame(rer_matrix), file.path(opt$outdir, "rer_matrix.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
path_names <- names(phenotype_paths)
if (is.null(path_names)) path_names <- colnames(rer_matrix)
write.table(data.frame(branch = path_names, value = as.numeric(phenotype_paths)),
            file.path(opt$outdir, "phenotype_paths.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
write.table(as.data.frame(association), file.path(opt$outdir, "associations.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(data.frame(
  metric = c("genes", "phenotype_species", "foreground_input", "foreground_used",
             "foreground_missing", "clade", "transition", "weighted", "transform",
             "impute", "min_trees", "max_trees", "min_species", "min_valid", "min_foreground",
             "bootstrap", "bootn"),
  value = c(nrow(rer_matrix), nrow(phenotype), length(foreground), length(available),
            length(missing_foreground), opt$clade, opt$transition, weighted, opt$transform,
            impute, min_trees, max_trees, min_species, min_valid, min_foreground, bootstrap, bootn)
), file.path(opt$outdir, "run_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (length(missing_foreground)) {
  writeLines(missing_foreground, file.path(opt$outdir, "missing_foreground_species.txt"))
}
saveRDS(list(trees = trees_obj, rer = rer_matrix, phenotype_tree = phenotype_tree,
             phenotype_paths = phenotype_paths, association = association),
        file.path(opt$outdir, "rerconverge_objects.rds"))

message("RERconverge completed: ", normalizePath(opt$outdir))
