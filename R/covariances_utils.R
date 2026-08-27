# Pooled variance estimator
estimate_covariance_pooled <- function(x, means, lengths) {
  n <- nrow(x)
  signal_dim <- ncol(x)
  n_segments <- nrow(means)
  residual_df <- n - n_segments

  if (residual_df < signal_dim) {
    log_n_stop(
      paste(
        "Not enough residual degrees of freedom to estimate",
        "a positive-definite pooled covariance matrix."
      )
    )
  }

  segment_id <- rep(
    seq_len(n_segments),
    times = lengths
  )

  fitted_values <- means[
    segment_id,
    ,
    drop = FALSE
  ]

  residuals <- x - fitted_values

  crossprod(residuals) / residual_df
}

# MoM based on first differences
estimate_covariance_diff_mom <- function(
    x,
    number_of_blocks = NULL,
    aggregation = c("entrywise", "geometric"),
    geometric_tolerance = 1e-8,
    geometric_max_iterations = 1000L,
    eigenvalue_floor = 1e-8
) {
  aggregation <- match.arg(aggregation)

  n <- nrow(x)
  signal_dim <- ncol(x)
  number_of_differences <- n - 1L

  if (number_of_differences < 1L) {
    log_n_stop(
      "At least two observations are required to compute first differences."
    )
  }

  if (is.null(number_of_blocks)) {
    number_of_blocks <- ceiling(
      sqrt(number_of_differences)
    )
  }

  if (length(number_of_blocks) != 1L ||
      !is.finite(number_of_blocks) ||
      number_of_blocks < 1 ||
      number_of_blocks != as.integer(number_of_blocks)) {
    log_n_stop(
      "covariance_blocks must be a positive integer."
    )
  }

  number_of_blocks <- as.integer(number_of_blocks)

  if (number_of_blocks > number_of_differences) {
    log_n_stop(
      paste(
        "covariance_blocks cannot exceed the number",
        "of available first differences."
      )
    )
  }

  # First differences.
  first_differences <-
    x[-1L, , drop = FALSE] -
    x[-n, , drop = FALSE]

  # Divide the first differences into consecutive blocks.
  block_indices <- split_into_consecutive_blocks(
    number_of_observations = number_of_differences,
    number_of_blocks = number_of_blocks
  )

  # Method-of-moments covariance estimator within each block.
  block_covariances <- lapply(
    block_indices,
    function(indices) {
      current_differences <- first_differences[
        indices,
        ,
        drop = FALSE
      ]

      crossprod(current_differences) /
        (2 * length(indices))
    }
  )

  # Robust aggregation of the block covariance estimators.
  covariance <- switch(
    aggregation,

    entrywise = {
      entrywise_matrix_median(
        block_covariances
      )
    },

    geometric = {
      geometric_matrix_median(
        matrices = block_covariances,
        tolerance = geometric_tolerance,
        max_iterations = geometric_max_iterations
      )
    }
  )

  # Always return a strictly positive-definite covariance matrix.
  covariance <- make_positive_definite(
    covariance = covariance,
    relative_floor = eigenvalue_floor
  )

  if (!all(
    dim(covariance) ==
    c(signal_dim, signal_dim)
  )) {
    log_n_stop(
      "Invalid covariance dimensions in estimate_covariance_diff_mom()."
    )
  }

  covariance
}

# utils: divide observations in blocks
split_into_consecutive_blocks <- function(
    number_of_observations,
    number_of_blocks
) {
  base_block_size <- number_of_observations %/% number_of_blocks
  remainder <- number_of_observations %% number_of_blocks

  block_sizes <- rep(
    base_block_size,
    number_of_blocks
  )

  if (remainder > 0L) {
    block_sizes[seq_len(remainder)] <-
      block_sizes[seq_len(remainder)] + 1L
  }

  block_endpoints <- cumsum(block_sizes)
  block_starts <- c(
    1L,
    head(block_endpoints, -1L) + 1L
  )

  Map(
    seq.int,
    from = block_starts,
    to = block_endpoints
  )
}

# Entrywise matrix median
entrywise_matrix_median <- function(matrices) {
  if (length(matrices) < 1L) {
    log_n_stop("At least one matrix is required.")
  }

  matrix_dimensions <- dim(matrices[[1L]])

  valid_dimensions <- vapply(
    matrices,
    function(current_matrix) {
      identical(dim(current_matrix), matrix_dimensions)
    },
    logical(1)
  )

  if (!all(valid_dimensions)) {
    log_n_stop("All matrices must have the same dimensions.")
  }

  matrix_array <- simplify2array(matrices)

  covariance <- apply(
    matrix_array,
    MARGIN = c(1L, 2L),
    FUN = stats::median
  )

  # Remove possible numerical asymmetry.
  (covariance + t(covariance)) / 2
}

# Geometric median in matrix space
geometric_matrix_median <- function(
    matrices,
    tolerance = 1e-8,
    max_iterations = 1000L
) {
  if (length(matrices) < 1L) {
    log_n_stop("At least one matrix is required.")
  }

  matrix_dimensions <- dim(matrices[[1L]])

  valid_dimensions <- vapply(
    matrices,
    function(current_matrix) {
      identical(dim(current_matrix), matrix_dimensions)
    },
    logical(1)
  )

  if (!all(valid_dimensions)) {
    log_n_stop("All matrices must have the same dimensions.")
  }

  # Start from the arithmetic mean.
  current_estimate <- Reduce(`+`, matrices) / length(matrices)

  converged <- FALSE

  for (iteration in seq_len(max_iterations)) {

    distances <- vapply(
      matrices,
      function(current_matrix) {
        sqrt(
          sum(
            (current_matrix - current_estimate)^2
          )
        )
      },
      numeric(1)
    )

    coincident <- distances <= tolerance
    noncoincident <- !coincident

    # ------------------------------------------------------------
    # Standard Weiszfeld step if the iterate is not at a data point.
    # ------------------------------------------------------------
    if (!any(coincident)) {

      weights <- 1 / distances

      new_estimate <- Reduce(
        `+`,
        Map(
          function(current_matrix, current_weight) {
            current_weight * current_matrix
          },
          matrices,
          weights
        )
      ) / sum(weights)

    } else {

      # ----------------------------------------------------------
      # Modified Weiszfeld step.
      #
      # If the current iterate coincides with one or more input
      # matrices, check the geometric-median optimality condition.
      # ----------------------------------------------------------

      coincident_multiplicity <- sum(coincident)

      if (!any(noncoincident)) {
        # All matrices coincide.
        converged <- TRUE
        break
      }

      noncoincident_matrices <- matrices[noncoincident]
      noncoincident_distances <- distances[noncoincident]

      residual_direction <- Reduce(
        `+`,
        Map(
          function(current_matrix, current_distance) {
            (current_matrix - current_estimate) /
              current_distance
          },
          noncoincident_matrices,
          noncoincident_distances
        )
      )

      residual_norm <- sqrt(
        sum(residual_direction^2)
      )

      # Subgradient optimality condition.
      if (residual_norm <= coincident_multiplicity) {
        converged <- TRUE
        break
      }

      # Ordinary Weiszfeld map excluding coincident matrices.
      weights <- 1 / noncoincident_distances

      weiszfeld_estimate <- Reduce(
        `+`,
        Map(
          function(current_matrix, current_weight) {
            current_weight * current_matrix
          },
          noncoincident_matrices,
          weights
        )
      ) / sum(weights)

      # Modified Weiszfeld correction.
      beta <- coincident_multiplicity / residual_norm

      new_estimate <-
        beta * current_estimate +
        (1 - beta) * weiszfeld_estimate
    }

    relative_change <- sqrt(
      sum(
        (new_estimate - current_estimate)^2
      )
    ) / max(
      1,
      sqrt(sum(current_estimate^2))
    )

    current_estimate <- new_estimate

    if (relative_change <= tolerance) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    log_warn(
      paste(
        "The geometric-median algorithm did not converge after",
        max_iterations,
        "iterations."
      )
    )
  }

  # Remove possible numerical asymmetry.
  (current_estimate + t(current_estimate)) / 2
}

# Correction for positive definitiveness
make_positive_definite <- function(
    covariance,
    relative_floor = 1e-8
) {
  # Remove possible numerical asymmetry.
  covariance <- (covariance + t(covariance)) / 2

  covariance_eigen <- eigen(
    covariance,
    symmetric = TRUE
  )

  eigenvalues <- covariance_eigen$values

  # Set the eigenvalue floor relative to the spectral scale
  # of the estimated covariance matrix.
  spectral_scale <- max(abs(eigenvalues))

  minimum_eigenvalue <- if (spectral_scale > 0) {
    relative_floor * spectral_scale
  } else {
    relative_floor
  }

  # Return unchanged if the matrix is already sufficiently
  # positive definite.
  if (all(eigenvalues >= minimum_eigenvalue)) {
    return(covariance)
  }

  log_warn(
    paste(
      "The covariance estimate was not sufficiently positive definite.",
      "Its eigenvalues were truncated."
    )
  )

  corrected_values <- pmax(
    eigenvalues,
    minimum_eigenvalue
  )

  corrected_covariance <-
    covariance_eigen$vectors %*%
    diag(
      corrected_values,
      nrow = length(corrected_values)
    ) %*%
    t(covariance_eigen$vectors)

  # Remove possible numerical asymmetry.
  (corrected_covariance + t(corrected_covariance)) / 2
}
