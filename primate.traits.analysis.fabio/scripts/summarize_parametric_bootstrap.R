#!/usr/bin/env Rscript

# Combine trait-wise bootstrap results and compare observed statistics with
# their conditional and full parametric null distributions.

parse_arguments <- function(arguments) {
  result <- list()
  for (index in seq(1L, length(arguments), by = 2L)) {
    if (index == length(arguments) || !startsWith(arguments[[index]], "--")) {
      stop("Invalid command-line arguments.", call. = FALSE)
    }
    result[[substring(arguments[[index]], 3L)]] <- arguments[[index + 1L]]
  }
  result
}

read_and_bind <- function(paths) {
  if (!length(paths)) stop("No input tables found.", call. = FALSE)
  do.call(rbind, lapply(paths, function(path) {
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  }))
}

write_tsv <- function(data, path) {
  utils::write.table(
    data, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}

collapse_repeated <- function(data, key_columns, label, tolerance = 1e-7) {
  key <- do.call(paste, c(data[key_columns], sep = "\r"))
  groups <- split(seq_len(nrow(data)), key)
  collapsed <- lapply(groups, function(indices) {
    group <- data[indices, , drop = FALSE]
    for (column in setdiff(names(group), key_columns)) {
      values <- group[[column]]
      if (is.numeric(values)) {
        if (length(unique(is.na(values))) > 1L) {
          stop(label, " differs between chunks in column ", column, call. = FALSE)
        }
        finite <- values[is.finite(values)]
        if (length(finite) > 1L) {
          scale <- max(1, abs(finite[[1L]]))
          if (max(abs(finite - finite[[1L]])) > tolerance * scale) {
            stop(label, " differs between chunks in column ", column, call. = FALSE)
          }
        }
      } else if (length(unique(values)) > 1L) {
        stop(label, " differs between chunks in column ", column, call. = FALSE)
      }
    }
    group[1L, , drop = FALSE]
  })
  result <- do.call(rbind, collapsed)
  rownames(result) <- NULL
  result
}

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
required <- c(
  "search-dir", "simulations-output", "observed-output", "fits-output",
  "summary-output", "stability-output", "expected-replicates",
  "expected-traits"
)
if (!all(required %in% names(arguments))) {
  stop("Missing summary arguments.", call. = FALSE)
}

simulation_files <- list.files(
  arguments$`search-dir`, pattern = "[.]bootstrap_simulations[.]tsv$",
  recursive = TRUE, full.names = TRUE
)
observed_files <- list.files(
  arguments$`search-dir`, pattern = "[.]bootstrap_observed[.]tsv$",
  recursive = TRUE, full.names = TRUE
)
fit_files <- list.files(
  arguments$`search-dir`, pattern = "[.]bootstrap_fit[.]tsv$",
  recursive = TRUE, full.names = TRUE
)
simulations <- read_and_bind(simulation_files)
observed <- read_and_bind(observed_files)
fits <- read_and_bind(fit_files)

expected_replicates <- as.integer(arguments$`expected-replicates`)
expected_traits <- as.integer(arguments$`expected-traits`)
if (!is.finite(expected_replicates) || expected_replicates < 1L ||
    !is.finite(expected_traits) || expected_traits < 1L) {
  stop("Expected counts must be positive integers.", call. = FALSE)
}

traits <- sort(unique(simulations$trait))
if (length(traits) != expected_traits) {
  stop(
    "Expected ", expected_traits, " traits but found ", length(traits), ".",
    call. = FALSE
  )
}

full_key <- paste(
  simulations$trait, simulations$bootstrap_type, simulations$simulation,
  simulations$tail, sep = "\r"
)
if (anyDuplicated(full_key)) {
  stop("Duplicate trait/type/simulation/tail rows found.", call. = FALSE)
}
expected_ids <- seq_len(expected_replicates)
simulation_groups <- split(
  simulations,
  interaction(
    simulations$trait, simulations$bootstrap_type, simulations$tail,
    drop = TRUE
  )
)
bad_groups <- vapply(simulation_groups, function(group) {
  !identical(sort(unique(as.integer(group$simulation))), expected_ids)
}, logical(1))
if (any(bad_groups)) {
  stop(
    "Incomplete simulation IDs in ", sum(bad_groups),
    " trait/bootstrap/tail groups.", call. = FALSE
  )
}

# Observed statistics and fitted models are repeated by every chunk. Very small
# optimizer-level numerical differences in fitted parameters are tolerated,
# while substantive disagreement stops aggregation.
observed <- collapse_repeated(observed, c("trait", "tail"), "Observed output")
fits <- collapse_repeated(fits, c("trait", "model"), "Fit output")
observed_counts <- table(observed$trait)
fit_counts <- table(fits$trait)
if (!identical(sort(names(observed_counts)), traits) ||
    any(observed_counts != 2L)) {
  stop("Observed output is not exactly two tails per trait.", call. = FALSE)
}
if (!identical(sort(names(fit_counts)), traits) || any(fit_counts != 2L)) {
  stop("Fit output is not exactly two models per trait.", call. = FALSE)
}

simulation_key <- paste(
  simulations$trait, simulations$bootstrap_type, simulations$simulation,
  sep = "\r"
)
model_rows <- simulations[!duplicated(simulation_key), c(
  "trait", "bootstrap_type", "simulation", "generating_model",
  "selected_model", "generation_attempt"
)]
stability_groups <- split(
  model_rows,
  interaction(model_rows$trait, model_rows$bootstrap_type, drop = TRUE)
)
stability <- do.call(rbind, lapply(stability_groups, function(group) {
  data.frame(
    trait = group$trait[[1L]],
    bootstrap_type = group$bootstrap_type[[1L]],
    generating_model = group$generating_model[[1L]],
    n_simulations = nrow(group),
    proportion_same_model = mean(
      group$selected_model == group$generating_model
    ),
    proportion_selected_BM = mean(group$selected_model == "BM"),
    proportion_selected_OU = mean(group$selected_model == "OU"),
    mean_generation_attempts = mean(group$generation_attempt),
    max_generation_attempts = max(group$generation_attempt),
    stringsAsFactors = FALSE
  )
}))
rownames(stability) <- NULL

metrics <- c(
  "family_excess", "superfamily_excess", "order_excess",
  "median_depth_quantile", "median_S", "median_Delta", "median_T",
  "median_score"
)
summary_rows <- list()
row_index <- 1L
group_keys <- unique(simulations[c("trait", "bootstrap_type", "tail")])
for (index in seq_len(nrow(group_keys))) {
  key <- group_keys[index, ]
  selected <- simulations$trait == key$trait &
    simulations$bootstrap_type == key$bootstrap_type &
    simulations$tail == key$tail
  group <- simulations[selected, , drop = FALSE]
  observed_row <- observed[
    observed$trait == key$trait & observed$tail == key$tail, , drop = FALSE
  ]
  if (nrow(observed_row) != 1L) {
    stop("Observed statistic is not unique for ", key$trait, "/", key$tail)
  }
  stability_row <- stability[
    stability$trait == key$trait &
      stability$bootstrap_type == key$bootstrap_type, , drop = FALSE
  ]
  for (metric in metrics) {
    simulated_values <- group[[metric]]
    simulated_values <- simulated_values[is.finite(simulated_values)]
    observed_value <- observed_row[[metric]]
    if (!is.finite(observed_value) || !length(simulated_values)) {
      null_mean <- null_median <- null_sd <- null_low <- null_high <- NA_real_
      empirical_percentile <- standardized_deviation <- empirical_p <- NA_real_
    } else {
      null_mean <- mean(simulated_values)
      null_median <- stats::median(simulated_values)
      null_sd <- stats::sd(simulated_values)
      limits <- stats::quantile(
        simulated_values, c(0.025, 0.975), type = 8, names = FALSE
      )
      null_low <- limits[[1L]]
      null_high <- limits[[2L]]
      empirical_percentile <- (
        sum(simulated_values < observed_value) +
          0.5 * sum(simulated_values == observed_value) + 0.5
      ) / (length(simulated_values) + 1)
      standardized_deviation <- if (is.finite(null_sd) && null_sd > 0) {
        (observed_value - null_mean) / null_sd
      } else {
        NA_real_
      }
      empirical_p <- min(
        1, 2 * min(empirical_percentile, 1 - empirical_percentile)
      )
    }
    summary_rows[[row_index]] <- data.frame(
      trait = key$trait,
      bootstrap_type = key$bootstrap_type,
      tail = key$tail,
      statistic = metric,
      generating_model = observed_row$generating_model,
      n_simulations = nrow(group),
      n_valid = length(simulated_values),
      observed_value = observed_value,
      null_mean = null_mean,
      null_median = null_median,
      null_sd = null_sd,
      null_interval_low = null_low,
      null_interval_high = null_high,
      empirical_percentile = empirical_percentile,
      empirical_two_sided_p = empirical_p,
      standardized_deviation = standardized_deviation,
      model_stability = stability_row$proportion_same_model,
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L
  }
}
summary_table <- do.call(rbind, summary_rows)

simulations <- simulations[order(
  simulations$trait, simulations$simulation,
  simulations$bootstrap_type, simulations$tail
), ]
observed <- observed[order(observed$trait, observed$tail), ]
fits <- fits[order(fits$trait, fits$model), ]
summary_table <- summary_table[order(
  summary_table$trait, summary_table$bootstrap_type,
  summary_table$tail, summary_table$statistic
), ]
stability <- stability[order(stability$trait, stability$bootstrap_type), ]

write_tsv(simulations, arguments$`simulations-output`)
write_tsv(observed, arguments$`observed-output`)
write_tsv(fits, arguments$`fits-output`)
write_tsv(summary_table, arguments$`summary-output`)
write_tsv(stability, arguments$`stability-output`)

message(
  "Combined ", length(unique(simulations$trait)), " traits and ",
  nrow(simulations), " simulation-tail rows."
)
