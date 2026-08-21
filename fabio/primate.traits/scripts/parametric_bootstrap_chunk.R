#!/usr/bin/env Rscript

# Chunked trait-wise parametric bootstrap for the analytical phylogenetic shift
# score. Simulation numbers are global, so independent Slurm chunks reproduce
# exactly the seeds and results of a single uninterrupted run.
# Each simulated trait is analysed twice: with the observed model/parameters
# held fixed (conditional bootstrap) and after BM/OU refitting and selection
# with the same AIC rule used by pss.core.R (full bootstrap).

parse_arguments <- function(arguments) {
  valued <- c(
    "--trait", "--trait-file", "--tree", "--taxonomy", "--pss-core",
    "--observed-classified", "--simulation-start", "--replicates",
    "--tail-proportion", "--seed",
    "--max-attempts", "--simulations-output", "--observed-output",
    "--fit-output"
  )
  result <- list()
  index <- 1L
  while (index <= length(arguments)) {
    key <- arguments[[index]]
    if (!key %in% valued || index == length(arguments)) {
      stop("Invalid or incomplete argument: ", key, call. = FALSE)
    }
    result[[substring(key, 3L)]] <- arguments[[index + 1L]]
    index <- index + 2L
  }
  missing <- substring(valued, 3L)[!substring(valued, 3L) %in% names(result)]
  if (length(missing)) {
    stop("Missing arguments: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  result
}

stable_string_hash <- function(value) {
  hash <- 0
  for (byte in utf8ToInt(enc2utf8(value))) {
    hash <- (hash * 131 + byte) %% 2147483000
  }
  as.integer(hash)
}

fit_models <- function(tree, values) {
  fits <- list(
    BM = suppressWarnings(
      geiger::fitContinuous(tree, values, model = "BM", ncores = 1)
    ),
    OU = suppressWarnings(
      geiger::fitContinuous(tree, values, model = "OU", ncores = 1)
    )
  )
  validate_fit(fits$BM, "BM")
  validate_fit(fits$OU, "OU")
  fits
}

covariances_from_fits <- function(tree, fits) {
  list(
    BM = bm_covariance(tree, fits$BM$opt$sigsq),
    OU = ou_covariance(tree, fits$OU$opt$alpha, fits$OU$opt$sigsq)
  )
}

simulate_gaussian_trait <- function(mean_value, covariance) {
  covariance <- (covariance + t(covariance)) / 2
  decomposition <- eigen(covariance, symmetric = TRUE)
  tolerance <- 1e-9 * max(1, max(abs(decomposition$values)))
  if (min(decomposition$values) < -tolerance) {
    stop("Generating covariance is not positive semidefinite.", call. = FALSE)
  }
  simulated <- as.numeric(
    mean_value + decomposition$vectors %*%
      (sqrt(pmax(decomposition$values, 0)) * stats::rnorm(nrow(covariance)))
  )
  stats::setNames(simulated, rownames(covariance))
}

load_taxonomic_levels <- function(taxonomy_file, species_1, species_2) {
  taxonomy <- utils::read.delim(
    taxonomy_file, check.names = FALSE, stringsAsFactors = FALSE
  )
  required <- c("species", "family", "primate_group_code")
  if (!all(required %in% names(taxonomy))) {
    stop("Taxonomy is missing required columns.", call. = FALSE)
  }
  if (anyDuplicated(taxonomy$species)) {
    stop("Taxonomy contains duplicated species.", call. = FALSE)
  }
  family <- stats::setNames(taxonomy$family, taxonomy$species)
  superfamily <- stats::setNames(
    taxonomy$primate_group_code, taxonomy$species
  )
  if (anyNA(family[c(species_1, species_2)]) ||
      anyNA(superfamily[c(species_1, species_2)])) {
    stop("At least one matched species is absent from taxonomy.", call. = FALSE)
  }
  ifelse(
    family[species_1] == family[species_2], "family",
    ifelse(
      superfamily[species_1] == superfamily[species_2],
      "superfamily", "order"
    )
  )
}

tail_statistics <- function(scores, selected_model, taxonomic_levels,
                            tail_proportion) {
  required <- c(
    "NormalizedTraitDifference", "NormalizedPatristicDistance",
    "Sbm", "Sou", "FinalScore"
  )
  if (!all(required %in% names(scores))) {
    stop("Score table is missing required columns.", call. = FALSE)
  }
  if (nrow(scores) != length(taxonomic_levels)) {
    stop("Taxonomic levels do not align with score rows.", call. = FALSE)
  }
  n_pairs <- nrow(scores)
  n_tail <- max(1L, ceiling(n_pairs * tail_proportion))
  ranking <- order(scores$FinalScore, decreasing = TRUE, method = "radix")
  depth_quantile <- (
    rank(scores$NormalizedPatristicDistance, ties.method = "average") - 0.5
  ) / n_pairs
  chosen_s <- if (selected_model == "BM") scores$Sbm else scores$Sou
  opportunity <- vapply(
    c("family", "superfamily", "order"),
    function(level) mean(taxonomic_levels == level), numeric(1)
  )
  make_row <- function(indices, tail_name) {
    observed_share <- vapply(
      names(opportunity),
      function(level) mean(taxonomic_levels[indices] == level), numeric(1)
    )
    excess <- observed_share - opportunity
    excess[opportunity == 0] <- NA_real_
    data.frame(
      tail = tail_name,
      n_pairs = n_pairs,
      n_tail = n_tail,
      family_excess = unname(excess[["family"]]),
      superfamily_excess = unname(excess[["superfamily"]]),
      order_excess = unname(excess[["order"]]),
      median_depth_quantile = stats::median(depth_quantile[indices]),
      median_S = stats::median(chosen_s[indices]),
      median_Delta = stats::median(
        scores$NormalizedTraitDifference[indices]
      ),
      median_T = stats::median(scores$NormalizedPatristicDistance[indices]),
      median_score = stats::median(scores$FinalScore[indices]),
      stringsAsFactors = FALSE
    )
  }
  rbind(
    make_row(ranking[seq_len(n_tail)], "top"),
    make_row(tail(ranking, n_tail), "bottom")
  )
}

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
trait <- arguments$trait
simulation_start <- as.integer(arguments$`simulation-start`)
replicates <- as.integer(arguments$replicates)
tail_proportion <- as.numeric(arguments$`tail-proportion`)
base_seed <- as.integer(arguments$seed)
max_attempts <- as.integer(arguments$`max-attempts`)
if (!is.finite(simulation_start) || simulation_start < 1L ||
    !is.finite(replicates) || replicates < 1L ||
    !is.finite(tail_proportion) || tail_proportion <= 0 ||
    tail_proportion >= 0.5 || !is.finite(base_seed) ||
    !is.finite(max_attempts) || max_attempts < 1L) {
  stop("Invalid numeric bootstrap setting.", call. = FALSE)
}

source(arguments$`pss-core`)
require_package("ape")
require_package("geiger")

trait_data <- utils::read.delim(
  arguments$`trait-file`, check.names = FALSE, stringsAsFactors = FALSE
)
if (!all(c("species", "trait_value") %in% names(trait_data))) {
  stop("Filtered trait input must contain species and trait_value.", call. = FALSE)
}
trait_data[[trait]] <- trait_data$trait_value
trait_data$trait_value <- NULL
matched <- match_trait_to_tree(
  read_trait_data(trait_data, "species"),
  read_phylogeny(arguments$tree), trait, "species"
)
tree <- matched$tree
observed_values <- stats::setNames(matched$data[[trait]], matched$data$species)
observed_values <- observed_values[tree$tip.label]

message("Fitting observed BM and OU models for ", trait, ".")
observed_fits <- fit_models(tree, observed_values)
observed_covariances <- covariances_from_fits(tree, observed_fits)
generating_model <- select_model(observed_fits)
generating_mean <- observed_fits[[generating_model]]$opt$z0
if (is.null(generating_mean) || !is.finite(generating_mean)) {
  stop("Selected model has an invalid fitted mean.", call. = FALSE)
}

observed_scores <- utils::read.delim(
  arguments$`observed-classified`, check.names = FALSE,
  stringsAsFactors = FALSE
)
required_observed <- c(
  "Species1", "Species2", "taxonomic_level", "SelectedModel"
)
if (!all(required_observed %in% names(observed_scores))) {
  stop("Observed classified result is missing required columns.", call. = FALSE)
}
recorded_models <- unique(observed_scores$SelectedModel)
if (length(recorded_models) != 1L || recorded_models != generating_model) {
  stop(
    "New observed fit does not reproduce the recorded selected model for ",
    trait, ".", call. = FALSE
  )
}
observed_stats <- tail_statistics(
  observed_scores, generating_model, observed_scores$taxonomic_level,
  tail_proportion
)
observed_stats <- cbind(
  data.frame(
    trait = trait, generating_model = generating_model,
    n_species = length(observed_values), stringsAsFactors = FALSE
  ),
  observed_stats
)

pair_template <- utils::combn(seq_along(tree$tip.label), 2L)
species_1 <- tree$tip.label[pair_template[1L, ]]
species_2 <- tree$tip.label[pair_template[2L, ]]
template_taxonomic_levels <- load_taxonomic_levels(
  arguments$taxonomy, species_1, species_2
)

fit_rows <- model_fit_table(observed_fits, trait, generating_model)
names(fit_rows) <- tolower(names(fit_rows))
fit_rows$n_species <- length(observed_values)
fit_rows$generating_model <- generating_model

RNGkind("L'Ecuyer-CMRG")
trait_seed <- stable_string_hash(trait)
simulation_ids <- seq.int(simulation_start, length.out = replicates)
simulation_rows <- vector("list", replicates * 2L)
row_index <- 1L
for (local_index in seq_along(simulation_ids)) {
  simulation <- simulation_ids[[local_index]]
  success <- FALSE
  last_error <- NULL
  for (attempt in seq_len(max_attempts)) {
    replicate_seed <- (
      as.double(base_seed) + trait_seed + simulation * 1009 + attempt * 1000003
    ) %% 2147483000
    set.seed(as.integer(replicate_seed))
    candidate <- tryCatch({
      simulated_values <- simulate_gaussian_trait(
        generating_mean, observed_covariances[[generating_model]]
      )
      conditional_scores <- calculate_pairwise_scores(
        simulated_values, tree, observed_covariances$BM,
        observed_covariances$OU, generating_model
      )
      conditional_levels <- load_taxonomic_levels(
        arguments$taxonomy,
        conditional_scores$Species1, conditional_scores$Species2
      )
      conditional_stats <- tail_statistics(
        conditional_scores, generating_model, conditional_levels,
        tail_proportion
      )

      simulated_fits <- fit_models(tree, simulated_values)
      selected_model <- select_model(simulated_fits)
      simulated_covariances <- covariances_from_fits(tree, simulated_fits)
      full_scores <- calculate_pairwise_scores(
        simulated_values, tree, simulated_covariances$BM,
        simulated_covariances$OU, selected_model
      )
      full_levels <- load_taxonomic_levels(
        arguments$taxonomy, full_scores$Species1, full_scores$Species2
      )
      full_stats <- tail_statistics(
        full_scores, selected_model, full_levels, tail_proportion
      )
      list(
        conditional = conditional_stats, full = full_stats,
        selected_model = selected_model,
        aic_bm = fit_aic(simulated_fits$BM),
        aic_ou = fit_aic(simulated_fits$OU)
      )
    }, error = function(error) {
      last_error <<- conditionMessage(error)
      NULL
    })
    if (!is.null(candidate)) {
      success <- TRUE
      break
    }
  }
  if (!success) {
    stop(
      "Simulation ", simulation, " failed after ", max_attempts,
      " attempts: ", last_error, call. = FALSE
    )
  }
  common <- data.frame(
    trait = trait, simulation = simulation, generating_model = generating_model,
    n_species = length(observed_values), generation_attempt = attempt,
    stringsAsFactors = FALSE
  )
  conditional <- cbind(
    common,
    data.frame(
      bootstrap_type = "conditional", selected_model = generating_model,
      AIC_BM = NA_real_, AIC_OU = NA_real_, stringsAsFactors = FALSE
    ),
    candidate$conditional
  )
  full <- cbind(
    common,
    data.frame(
      bootstrap_type = "full", selected_model = candidate$selected_model,
      AIC_BM = candidate$aic_bm, AIC_OU = candidate$aic_ou,
      stringsAsFactors = FALSE
    ),
    candidate$full
  )
  simulation_rows[[row_index]] <- conditional
  simulation_rows[[row_index + 1L]] <- full
  row_index <- row_index + 2L
  if (local_index %% 10L == 0L || local_index == replicates) {
    message(
      "Completed ", local_index, "/", replicates, " simulations for ", trait,
      " (global simulation ", simulation, ")."
    )
  }
}

simulation_table <- do.call(rbind, simulation_rows)
expected_rows <- replicates * 2L * 2L
if (nrow(simulation_table) != expected_rows) {
  stop("Unexpected number of simulation rows.", call. = FALSE)
}
write_table <- function(data, path) {
  utils::write.table(
    data, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}
write_table(simulation_table, arguments$`simulations-output`)
write_table(observed_stats, arguments$`observed-output`)
write_table(fit_rows, arguments$`fit-output`)

message(
  "Bootstrap chunk complete for ", trait, ": simulations ",
  min(simulation_ids), "-", max(simulation_ids), "."
)
