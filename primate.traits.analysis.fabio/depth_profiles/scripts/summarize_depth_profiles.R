#!/usr/bin/env Rscript

# Validate compact depth-profile chunks and aggregate observed and bootstrap
# empirical CDFs. Inference uses pointwise intervals, an unstudentized
# simultaneous maximum-deviation envelope, calibrated global curve statistics,
# and BH correction across traits within bootstrap type and tail.

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

write_tsv <- function(data, path) {
  utils::write.table(
    data, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}

write_tsv_gz <- function(data, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  utils::write.table(
    data, connection, sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA"
  )
}

curve_from_depths <- function(depths, grid) {
  if (!length(depths) || any(!is.finite(depths)) ||
      any(depths < 0 | depths > 1)) {
    stop("Invalid depth vector.", call. = FALSE)
  }
  curve <- vapply(grid, function(point) mean(depths <= point), numeric(1))
  if (any(diff(curve) < -1e-14) || any(curve < 0 | curve > 1) ||
      abs(tail(curve, 1L) - 1) > 1e-12) {
    stop("Invalid cumulative depth curve.", call. = FALSE)
  }
  curve
}

integration_weights <- function(grid) {
  differences <- diff(grid)
  if (length(grid) < 2L || any(!is.finite(grid)) ||
      any(differences <= 0)) {
    stop("Depth grid must be finite and strictly increasing.", call. = FALSE)
  }
  c(
    differences[[1L]] / 2,
    (head(differences, -1L) + tail(differences, -1L)) / 2,
    tail(differences, 1L) / 2
  )
}

empirical_upper_p <- function(null_values, observed_value) {
  (1 + sum(null_values >= observed_value)) / (length(null_values) + 1)
}

empirical_lower_p <- function(null_values, observed_value) {
  (1 + sum(null_values <= observed_value)) / (length(null_values) + 1)
}

empirical_two_sided_p <- function(null_values, observed_value) {
  empirical_upper_p(abs(null_values), abs(observed_value))
}

empirical_percentile <- function(null_values, observed_value) {
  (
    sum(null_values < observed_value) +
      0.5 * sum(null_values == observed_value) + 0.5
  ) / (length(null_values) + 1)
}

standardized_deviation <- function(null_values, observed_value) {
  deviation <- stats::sd(null_values)
  if (!is.finite(deviation) || deviation <= 0) return(NA_real_)
  (observed_value - mean(null_values)) / deviation
}

extract_excursions <- function(observed, lower, upper, grid) {
  direction <- ifelse(observed < lower, -1L, ifelse(observed > upper, 1L, 0L))
  encoded <- rle(direction)
  run_end <- cumsum(encoded$lengths)
  run_start <- c(1L, head(run_end, -1L) + 1L)
  selected <- which(encoded$values != 0L)
  empty <- list(
    crossed = FALSE,
    n_excursions = 0L,
    total_excursion_length = 0,
    first_excursion_start = NA_real_,
    first_excursion_end = NA_real_,
    first_excursion_sign = NA_character_,
    longest_excursion_start = NA_real_,
    longest_excursion_end = NA_real_,
    longest_excursion_length = 0,
    longest_excursion_sign = NA_character_,
    excursion_signs = NA_character_
  )
  if (!length(selected)) return(empty)

  starts <- grid[run_start[selected]]
  ends <- grid[run_end[selected]]
  lengths <- ends - starts
  signs <- ifelse(encoded$values[selected] < 0L, "deep", "shallow")
  longest <- which.max(lengths)
  list(
    crossed = TRUE,
    n_excursions = length(selected),
    total_excursion_length = sum(lengths),
    first_excursion_start = starts[[1L]],
    first_excursion_end = ends[[1L]],
    first_excursion_sign = signs[[1L]],
    longest_excursion_start = starts[[longest]],
    longest_excursion_end = ends[[longest]],
    longest_excursion_length = lengths[[longest]],
    longest_excursion_sign = signs[[longest]],
    excursion_signs = paste(unique(signs), collapse = ";")
  )
}

validate_profile_object <- function(object, path, expected_tail_proportion) {
  required <- c(
    "schema_version", "trait", "simulation_start", "simulation_end",
    "tail_proportion", "simulations", "simulation_depths", "observed",
    "observed_depths", "observed_pair_keys", "fits"
  )
  if (!is.list(object) || !all(required %in% names(object))) {
    stop("Invalid profile object: ", path, call. = FALSE)
  }
  if (!identical(as.integer(object$schema_version), 1L)) {
    stop("Unsupported profile schema in ", path, call. = FALSE)
  }
  if (!isTRUE(all.equal(
    as.numeric(object$tail_proportion), expected_tail_proportion,
    tolerance = 1e-14
  ))) {
    stop("Tail proportion differs in ", path, call. = FALSE)
  }
  simulations <- object$simulations
  if (nrow(simulations) != length(object$simulation_depths)) {
    stop("Depth vectors do not align with simulation rows in ", path)
  }
  if (nrow(object$observed) != length(object$observed_depths) ||
      nrow(object$observed) != length(object$observed_pair_keys)) {
    stop("Observed depth vectors do not align in ", path)
  }
  simulation_lengths <- vapply(object$simulation_depths, length, integer(1))
  simulation_medians <- vapply(
    object$simulation_depths, stats::median, numeric(1)
  )
  if (any(simulation_lengths != simulations$n_tail) ||
      any(abs(
        simulation_medians - simulations$median_depth_quantile
      ) > 1e-12)) {
    stop("Simulation-depth QC failed in ", path, call. = FALSE)
  }
  observed_lengths <- vapply(object$observed_depths, length, integer(1))
  observed_medians <- vapply(object$observed_depths, stats::median, numeric(1))
  if (any(observed_lengths != object$observed$n_tail) ||
      any(abs(
        observed_medians - object$observed$median_depth_quantile
      ) > 1e-12)) {
    stop("Observed-depth QC failed in ", path, call. = FALSE)
  }
  invisible(TRUE)
}

validate_repeated_table <- function(tables, key_columns, label,
                                    relative_tolerance = 1e-3,
                                    absolute_tolerance = 1e-8) {
  combined <- do.call(rbind, tables)
  key <- do.call(paste, c(combined[key_columns], sep = "\r"))
  groups <- split(seq_len(nrow(combined)), key)
  collapsed <- lapply(groups, function(indices) {
    group <- combined[indices, , drop = FALSE]
    for (column in setdiff(names(group), key_columns)) {
      values <- group[[column]]
      if (is.numeric(values)) {
        if (length(unique(is.na(values))) > 1L) {
          stop(label, " differs between chunks in ", column, call. = FALSE)
        }
        finite <- values[is.finite(values)]
        if (length(finite) > 1L) {
          scale <- max(abs(finite))
          maximum_difference <- max(abs(finite - finite[[1L]]))
          allowed <- absolute_tolerance + relative_tolerance * scale
          if (maximum_difference > allowed) {
            stop(label, " differs between chunks in ", column, call. = FALSE)
          }
        }
      } else if (length(unique(values)) > 1L) {
        stop(label, " differs between chunks in ", column, call. = FALSE)
      }
    }
    group[1L, , drop = FALSE]
  })
  result <- do.call(rbind, collapsed)
  rownames(result) <- NULL
  result
}

adjust_grouped <- function(probabilities, groups) {
  adjusted <- rep(NA_real_, length(probabilities))
  for (indices in split(seq_along(probabilities), groups)) {
    valid <- indices[is.finite(probabilities[indices])]
    if (length(valid)) {
      adjusted[valid] <- stats::p.adjust(probabilities[valid], method = "BH")
    }
  }
  adjusted
}

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
required_arguments <- c(
  "search-dir", "expected-replicates", "expected-traits",
  "tail-proportion", "grid-points", "replicates-output",
  "observed-output", "fits-output", "curves-output", "summary-output"
)
if (!all(required_arguments %in% names(arguments))) {
  stop("Missing depth-profile summary arguments.", call. = FALSE)
}

expected_replicates <- as.integer(arguments$`expected-replicates`)
expected_traits <- as.integer(arguments$`expected-traits`)
tail_proportion <- as.numeric(arguments$`tail-proportion`)
grid_points <- as.integer(arguments$`grid-points`)
if (!is.finite(expected_replicates) || expected_replicates < 1L ||
    !is.finite(expected_traits) || expected_traits < 1L ||
    !is.finite(tail_proportion) || tail_proportion <= 0 ||
    tail_proportion >= 0.5 || !is.finite(grid_points) || grid_points < 3L) {
  stop("Invalid numeric aggregation setting.", call. = FALSE)
}
depth_grid <- seq(0, 1, length.out = grid_points)
area_weights <- integration_weights(depth_grid)

profile_files <- list.files(
  arguments$`search-dir`, pattern = "[.]depth_profiles[.]rds$",
  recursive = TRUE, full.names = TRUE
)
if (!length(profile_files)) {
  stop("No depth-profile RDS chunks found.", call. = FALSE)
}
trait_names <- sub(
  "[.]chunk_[0-9]+[.]depth_profiles[.]rds$", "", basename(profile_files)
)
if (any(trait_names == basename(profile_files))) {
  stop("At least one chunk filename does not follow the expected pattern.")
}
files_by_trait <- split(profile_files, trait_names)
traits <- sort(names(files_by_trait))
if (length(traits) != expected_traits) {
  stop(
    "Expected ", expected_traits, " traits but found ", length(traits), ".",
    call. = FALSE
  )
}

replicate_rows <- list()
observed_rows <- list()
fit_rows <- list()
curve_rows <- list()
summary_rows <- list()
replicate_index <- observed_index <- fit_index <- 1L
curve_index <- summary_index <- 1L

for (trait in traits) {
  paths <- sort(files_by_trait[[trait]])
  objects <- lapply(paths, readRDS)
  for (index in seq_along(objects)) {
    validate_profile_object(objects[[index]], paths[[index]], tail_proportion)
    if (!identical(objects[[index]]$trait, trait)) {
      stop("Trait name does not match filename: ", paths[[index]])
    }
  }

  simulations <- do.call(rbind, lapply(objects, `[[`, "simulations"))
  simulation_depths <- unlist(
    lapply(objects, `[[`, "simulation_depths"), recursive = FALSE
  )
  rownames(simulations) <- NULL
  key <- paste(
    simulations$trait, simulations$bootstrap_type,
    simulations$simulation, simulations$tail, sep = "\r"
  )
  if (anyDuplicated(key)) {
    stop("Duplicated simulation profile for ", trait, call. = FALSE)
  }
  expected_ids <- seq_len(expected_replicates)
  groups <- split(
    simulations,
    interaction(simulations$bootstrap_type, simulations$tail, drop = TRUE)
  )
  complete <- vapply(groups, function(group) {
    identical(sort(as.integer(group$simulation)), expected_ids)
  }, logical(1))
  if (length(groups) != 4L || any(!complete)) {
    stop("Incomplete simulation IDs for ", trait, call. = FALSE)
  }

  observed <- validate_repeated_table(
    lapply(objects, `[[`, "observed"), c("trait", "tail"),
    paste0("Observed output for ", trait)
  )
  observed <- observed[match(c("top", "bottom"), observed$tail), , drop = FALSE]
  if (anyNA(observed$tail)) {
    stop("Observed tails are incomplete for ", trait, call. = FALSE)
  }
  observed_depths <- objects[[1L]]$observed_depths
  observed_pair_keys <- objects[[1L]]$observed_pair_keys
  for (object in objects[-1L]) {
    if (!isTRUE(all.equal(
      object$observed_depths, observed_depths, tolerance = 1e-14
    )) || !identical(object$observed_pair_keys, observed_pair_keys)) {
      stop("Observed selected pairs differ between chunks for ", trait)
    }
  }
  fits <- validate_repeated_table(
    lapply(objects, `[[`, "fits"), c("trait", "model"),
    paste0("Observed fits for ", trait)
  )

  replicate_rows[[replicate_index]] <- simulations
  replicate_index <- replicate_index + 1L
  observed_rows[[observed_index]] <- observed
  observed_index <- observed_index + 1L
  fit_rows[[fit_index]] <- fits
  fit_index <- fit_index + 1L

  for (bootstrap_type in c("conditional", "full")) {
    model_rows <- simulations[
      simulations$bootstrap_type == bootstrap_type &
        simulations$tail == "top", , drop = FALSE
    ]
    model_stability <- mean(
      model_rows$selected_model == model_rows$generating_model
    )
    for (tail_name in c("top", "bottom")) {
      selected <- which(
        simulations$bootstrap_type == bootstrap_type &
          simulations$tail == tail_name
      )
      selected <- selected[order(simulations$simulation[selected])]
      null_curves <- t(vapply(
        simulation_depths[selected], curve_from_depths,
        numeric(length(depth_grid)), grid = depth_grid
      ))
      observed_index_in_object <- match(tail_name, observed$tail)
      observed_curve <- curve_from_depths(
        observed_depths[[observed_index_in_object]], depth_grid
      )

      null_center <- apply(null_curves, 2L, stats::median)
      point_low <- apply(
        null_curves, 2L, stats::quantile,
        probs = 0.025, type = 8, names = FALSE
      )
      point_high <- apply(
        null_curves, 2L, stats::quantile,
        probs = 0.975, type = 8, names = FALSE
      )
      null_deviation <- sweep(null_curves, 2L, null_center, "-")
      null_max_abs <- apply(abs(null_deviation), 1L, max)
      critical <- unname(stats::quantile(
        null_max_abs, probs = 0.95, type = 8, names = FALSE
      ))
      simultaneous_low <- pmax(0, null_center - critical)
      simultaneous_high <- pmin(1, null_center + critical)

      observed_deviation <- observed_curve - null_center
      observed_signed_area <- sum(observed_deviation * area_weights)
      observed_absolute_area <- sum(abs(observed_deviation) * area_weights)
      observed_max_abs <- max(abs(observed_deviation))
      observed_max_deep <- max(-observed_deviation)
      observed_max_shallow <- max(observed_deviation)

      null_signed_area <- as.vector(null_deviation %*% area_weights)
      null_absolute_area <- as.vector(abs(null_deviation) %*% area_weights)
      null_max_deep <- apply(-null_deviation, 1L, max)
      null_max_shallow <- apply(null_deviation, 1L, max)
      excursion <- extract_excursions(
        observed_curve, simultaneous_low, simultaneous_high, depth_grid
      )
      resolution <- if (observed$n_tail[[observed_index_in_object]] < 5L) {
        "very_low"
      } else if (observed$n_tail[[observed_index_in_object]] < 10L) {
        "cautious"
      } else {
        "shape_suitable"
      }

      curve_rows[[curve_index]] <- data.frame(
        trait = trait,
        bootstrap_type = bootstrap_type,
        tail = tail_name,
        tail_proportion = tail_proportion,
        generating_model = observed$generating_model[[observed_index_in_object]],
        model_stability = model_stability,
        n_species = observed$n_species[[observed_index_in_object]],
        n_pairs = observed$n_pairs[[observed_index_in_object]],
        n_tail = observed$n_tail[[observed_index_in_object]],
        resolution = resolution,
        depth_grid = depth_grid,
        observed_curve = observed_curve,
        null_median = null_center,
        pointwise_low = point_low,
        pointwise_high = point_high,
        simultaneous_low = simultaneous_low,
        simultaneous_high = simultaneous_high,
        deviation = observed_deviation,
        outside_simultaneous = observed_curve < simultaneous_low |
          observed_curve > simultaneous_high,
        excursion_direction = ifelse(
          observed_curve < simultaneous_low, "deep",
          ifelse(observed_curve > simultaneous_high, "shallow", "inside")
        ),
        stringsAsFactors = FALSE
      )
      curve_index <- curve_index + 1L

      summary_rows[[summary_index]] <- data.frame(
        trait = trait,
        bootstrap_type = bootstrap_type,
        tail = tail_name,
        tail_proportion = tail_proportion,
        generating_model = observed$generating_model[[observed_index_in_object]],
        model_stability = model_stability,
        n_simulations = nrow(null_curves),
        n_species = observed$n_species[[observed_index_in_object]],
        n_pairs = observed$n_pairs[[observed_index_in_object]],
        n_tail = observed$n_tail[[observed_index_in_object]],
        resolution = resolution,
        signed_area = observed_signed_area,
        absolute_area = observed_absolute_area,
        max_absolute_deviation = observed_max_abs,
        max_deep_deviation = observed_max_deep,
        max_shallow_deviation = observed_max_shallow,
        signed_area_percentile = empirical_percentile(
          null_signed_area, observed_signed_area
        ),
        signed_area_standardized_deviation = standardized_deviation(
          null_signed_area, observed_signed_area
        ),
        p_signed_area_two_sided = empirical_two_sided_p(
          null_signed_area, observed_signed_area
        ),
        p_signed_area_deep = empirical_lower_p(
          null_signed_area, observed_signed_area
        ),
        p_signed_area_shallow = empirical_upper_p(
          null_signed_area, observed_signed_area
        ),
        p_absolute_area = empirical_upper_p(
          null_absolute_area, observed_absolute_area
        ),
        p_max_absolute_deviation = empirical_upper_p(
          null_max_abs, observed_max_abs
        ),
        p_max_deep_deviation = empirical_upper_p(
          null_max_deep, observed_max_deep
        ),
        p_max_shallow_deviation = empirical_upper_p(
          null_max_shallow, observed_max_shallow
        ),
        simultaneous_critical_value = critical,
        crossed_simultaneous_band = excursion$crossed,
        n_excursions = excursion$n_excursions,
        total_excursion_length = excursion$total_excursion_length,
        first_excursion_start = excursion$first_excursion_start,
        first_excursion_end = excursion$first_excursion_end,
        first_excursion_sign = excursion$first_excursion_sign,
        longest_excursion_start = excursion$longest_excursion_start,
        longest_excursion_end = excursion$longest_excursion_end,
        longest_excursion_length = excursion$longest_excursion_length,
        longest_excursion_sign = excursion$longest_excursion_sign,
        excursion_signs = excursion$excursion_signs,
        dominant_direction = if (
          observed_signed_area < 0
        ) "deep" else if (observed_signed_area > 0) "shallow" else "neutral",
        stringsAsFactors = FALSE
      )
      summary_index <- summary_index + 1L
    }
  }
  message("Aggregated depth profiles for ", trait, ".")
}

replicates <- do.call(rbind, replicate_rows)
observed <- do.call(rbind, observed_rows)
fits <- do.call(rbind, fit_rows)
curves <- do.call(rbind, curve_rows)
summary_table <- do.call(rbind, summary_rows)
rownames(replicates) <- rownames(observed) <- rownames(fits) <- NULL
rownames(curves) <- rownames(summary_table) <- NULL

fdr_groups <- interaction(
  summary_table$bootstrap_type, summary_table$tail, drop = TRUE
)
p_columns <- c(
  "p_signed_area_two_sided", "p_signed_area_deep",
  "p_signed_area_shallow", "p_absolute_area",
  "p_max_absolute_deviation", "p_max_deep_deviation",
  "p_max_shallow_deviation"
)
for (column in p_columns) {
  q_column <- sub("^p_", "q_", column)
  summary_table[[q_column]] <- adjust_grouped(
    summary_table[[column]], fdr_groups
  )
}

replicates <- replicates[order(
  replicates$trait, replicates$simulation,
  replicates$bootstrap_type, replicates$tail
), ]
observed <- observed[order(observed$trait, observed$tail), ]
fits <- fits[order(fits$trait, fits$model), ]
curves <- curves[order(
  curves$trait, curves$bootstrap_type, curves$tail, curves$depth_grid
), ]
summary_table <- summary_table[order(
  summary_table$trait, summary_table$bootstrap_type, summary_table$tail
), ]

write_tsv_gz(replicates, arguments$`replicates-output`)
write_tsv(observed, arguments$`observed-output`)
write_tsv(fits, arguments$`fits-output`)
write_tsv(curves, arguments$`curves-output`)
write_tsv(summary_table, arguments$`summary-output`)

message(
  "Combined ", length(traits), " traits, ", nrow(replicates),
  " simulation-tail rows, and ", nrow(curves), " plotting rows."
)
