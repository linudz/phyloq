#!/usr/bin/env Rscript

# Chunked conditional and full parametric bootstrap for continuous
# phylogenetic-depth profiles of PSS tails. The simulation seed formula and
# PSS calculations intentionally match the validated bootstrap workflow.
# Instead of discarding the selected tails, this script retains their fixed,
# opportunity-relative patristic-depth percentiles in a compact RDS file.

parse_arguments <- function(arguments) {
  valued <- c(
    "--trait", "--trait-file", "--tree", "--taxonomy", "--pss-core",
    "--observed-fits",
    "--observed-classified", "--simulation-start", "--replicates",
    "--tail-proportion", "--seed", "--max-attempts", "--profile-output"
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

load_fixed_observed_fits <- function(path, trait, n_species) {
  table <- utils::read.delim(
    path, check.names = FALSE, stringsAsFactors = FALSE
  )
  required <- c(
    "trait", "model", "sigma2", "alpha", "aic", "selected",
    "n_species", "generating_model"
  )
  if (!all(required %in% names(table))) {
    stop("Frozen observed-fit table is missing required columns.", call. = FALSE)
  }
  selected <- table[table$trait == trait, required, drop = FALSE]
  selected <- selected[match(c("BM", "OU"), selected$model), , drop = FALSE]
  if (nrow(selected) != 2L || anyNA(selected$model) ||
      !identical(selected$model, c("BM", "OU")) ||
      sum(selected$selected) != 1L ||
      length(unique(selected$generating_model)) != 1L ||
      selected$generating_model[[1L]] != selected$model[selected$selected] ||
      any(selected$n_species != n_species) ||
      any(!is.finite(selected$sigma2)) || any(selected$sigma2 <= 0) ||
      !is.finite(selected$alpha[selected$model == "OU"]) ||
      selected$alpha[selected$model == "OU"] < 0) {
    stop("Invalid frozen observed fits for ", trait, ".", call. = FALSE)
  }
  selected
}

estimate_process_mean <- function(values, covariance) {
  covariance <- (covariance + t(covariance)) / 2
  one <- rep(1, length(values))
  weights <- tryCatch(
    solve(covariance, one),
    error = function(error) qr.solve(covariance, one, tol = 1e-12)
  )
  denominator <- sum(weights)
  estimate <- sum(weights * as.numeric(values)) / denominator
  if (!is.finite(estimate) || !is.finite(denominator) || denominator == 0) {
    estimate <- mean(values)
  }
  as.numeric(estimate)
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

canonical_pair_key <- function(species_1, species_2) {
  first <- ifelse(species_1 <= species_2, species_1, species_2)
  second <- ifelse(species_1 <= species_2, species_2, species_1)
  paste(first, second, sep = "\r")
}

tail_profile <- function(scores, selected_model, taxonomic_levels,
                         tail_proportion) {
  required <- c(
    "Species1", "Species2", "NormalizedTraitDifference",
    "NormalizedPatristicDistance", "Sbm", "Sou", "FinalScore"
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

  make_tail <- function(indices, tail_name) {
    observed_share <- vapply(
      names(opportunity),
      function(level) mean(taxonomic_levels[indices] == level), numeric(1)
    )
    excess <- observed_share - opportunity
    excess[opportunity == 0] <- NA_real_
    depths <- as.numeric(depth_quantile[indices])
    if (length(depths) != n_tail || any(!is.finite(depths)) ||
        any(depths < 0 | depths > 1)) {
      stop("Invalid selected depth percentiles.", call. = FALSE)
    }
    list(
      statistics = data.frame(
        tail = tail_name,
        n_pairs = n_pairs,
        n_tail = n_tail,
        family_excess = unname(excess[["family"]]),
        superfamily_excess = unname(excess[["superfamily"]]),
        order_excess = unname(excess[["order"]]),
        median_depth_quantile = stats::median(depths),
        median_S = stats::median(chosen_s[indices]),
        median_Delta = stats::median(
          scores$NormalizedTraitDifference[indices]
        ),
        median_T = stats::median(
          scores$NormalizedPatristicDistance[indices]
        ),
        median_score = stats::median(scores$FinalScore[indices]),
        stringsAsFactors = FALSE
      ),
      depths = depths,
      pair_keys = canonical_pair_key(
        scores$Species1[indices], scores$Species2[indices]
      )
    )
  }

  top <- make_tail(ranking[seq_len(n_tail)], "top")
  bottom <- make_tail(tail(ranking, n_tail), "bottom")
  list(
    statistics = rbind(top$statistics, bottom$statistics),
    depths = list(top$depths, bottom$depths),
    pair_keys = list(top$pair_keys, bottom$pair_keys)
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

message("Loading frozen observed BM and OU fits for ", trait, ".")
fit_rows <- load_fixed_observed_fits(
  arguments$`observed-fits`, trait, length(observed_values)
)
generating_model <- fit_rows$generating_model[[1L]]
observed_covariances <- list(
  BM = bm_covariance(
    tree, fit_rows$sigma2[fit_rows$model == "BM"]
  ),
  OU = ou_covariance(
    tree,
    fit_rows$alpha[fit_rows$model == "OU"],
    fit_rows$sigma2[fit_rows$model == "OU"]
  )
)
generating_mean <- estimate_process_mean(
  observed_values, observed_covariances[[generating_model]]
)

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
    "Frozen generating model does not match the recorded selected model for ",
    trait, ".", call. = FALSE
  )
}
observed_profile <- tail_profile(
  observed_scores, generating_model, observed_scores$taxonomic_level,
  tail_proportion
)
observed_table <- cbind(
  data.frame(
    trait = trait, generating_model = generating_model,
    n_species = length(observed_values), stringsAsFactors = FALSE
  ),
  observed_profile$statistics
)

RNGkind("L'Ecuyer-CMRG")
trait_seed <- stable_string_hash(trait)
simulation_ids <- seq.int(simulation_start, length.out = replicates)
simulation_rows <- vector("list", replicates * 2L)
simulation_depths <- vector("list", replicates * 4L)
row_index <- 1L
depth_index <- 1L

for (local_index in seq_along(simulation_ids)) {
  simulation <- simulation_ids[[local_index]]
  success <- FALSE
  last_error <- NULL
  for (attempt in seq_len(max_attempts)) {
    replicate_seed <- (
      as.double(base_seed) + trait_seed + simulation * 1009 +
        attempt * 1000003
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
      conditional_profile <- tail_profile(
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
      full_profile <- tail_profile(
        full_scores, selected_model, full_levels, tail_proportion
      )
      list(
        conditional = conditional_profile,
        full = full_profile,
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
    candidate$conditional$statistics
  )
  full <- cbind(
    common,
    data.frame(
      bootstrap_type = "full", selected_model = candidate$selected_model,
      AIC_BM = candidate$aic_bm, AIC_OU = candidate$aic_ou,
      stringsAsFactors = FALSE
    ),
    candidate$full$statistics
  )
  simulation_rows[[row_index]] <- conditional
  simulation_rows[[row_index + 1L]] <- full
  row_index <- row_index + 2L

  simulation_depths[depth_index:(depth_index + 3L)] <- c(
    candidate$conditional$depths, candidate$full$depths
  )
  depth_index <- depth_index + 4L

  if (local_index %% 10L == 0L || local_index == replicates) {
    message(
      "Completed ", local_index, "/", replicates, " simulations for ", trait,
      " (global simulation ", simulation, ")."
    )
  }
}

simulation_table <- do.call(rbind, simulation_rows)
expected_rows <- replicates * 2L * 2L
if (nrow(simulation_table) != expected_rows ||
    length(simulation_depths) != expected_rows) {
  stop("Unexpected number of simulation profile rows.", call. = FALSE)
}
depth_lengths <- vapply(simulation_depths, length, integer(1))
depth_medians <- vapply(simulation_depths, stats::median, numeric(1))
if (any(depth_lengths != simulation_table$n_tail) ||
    any(abs(depth_medians - simulation_table$median_depth_quantile) > 1e-12)) {
  stop("Selected-depth vectors do not reproduce tail statistics.", call. = FALSE)
}

profile_object <- list(
  schema_version = 1L,
  depth_definition = paste0(
    "(average rank of NormalizedPatristicDistance - 0.5) / n_pairs"
  ),
  trait = trait,
  simulation_start = min(simulation_ids),
  simulation_end = max(simulation_ids),
  tail_proportion = tail_proportion,
  simulations = simulation_table,
  simulation_depths = simulation_depths,
  observed = observed_table,
  observed_depths = observed_profile$depths,
  observed_pair_keys = observed_profile$pair_keys,
  fits = fit_rows
)
saveRDS(profile_object, arguments$`profile-output`, compress = "gzip")

message(
  "Depth-profile chunk complete for ", trait, ": simulations ",
  min(simulation_ids), "-", max(simulation_ids), "."
)
