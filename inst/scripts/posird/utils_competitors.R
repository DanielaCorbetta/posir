# Utils for data simulation
# Function to get changepoints at at least ndelta distance
generate_changepoints <- function(
    n,
    n_chg,
    delta
) {
  n_segments <- n_chg + 1
  min_length <- ceiling(delta * n)

  remaining_length <-
    n - n_segments * min_length

  if (remaining_length < 0) {
    stop(
      paste(
        "The requested number of changepoints is incompatible",
        "with n and delta."
      )
    )
  }

  # Give each segment the minimum admissible length,
  # then randomly distribute the remaining observations.
  segment_lengths <- min_length +
    as.vector(
      rmultinom(
        n = 1,
        size = remaining_length,
        prob = rep(1, n_segments)
      )
    )

  # Changepoints are the final observations of each segment.
  changepoints <- cumsum(
    segment_lengths
  )[-n_segments]

  list(
    changepoints = changepoints,
    lengths = segment_lengths
  )
}

# Function to generate normal with piecewise constant mean, with jumps at changepoints
generate_piecewise_normal_2d <- function(
    n,
    changepoints,
    Sigma,
    possible_states = c(-2, 0, 2)
) {
  if (!all(dim(Sigma) == c(2, 2))) {
    stop("Sigma must be a 2 x 2 matrix.")
  }

  boundaries <- c(
    0,
    changepoints,
    n
  )

  n_segments <- length(boundaries) - 1

  # Generate piecewise-constant means.
  # Both coordinates change at every changepoint.
  means <- matrix(
    NA_real_,
    nrow = n_segments,
    ncol = 2
  )

  means[1, ] <- sample(
    possible_states,
    size = 2,
    replace = TRUE
  )

  if (n_segments > 1) {
    for (segment in 2:n_segments) {
      means[segment, 1] <- sample(
        setdiff(
          possible_states,
          means[segment - 1, 1]
        ),
        size = 1
      )

      means[segment, 2] <- sample(
        setdiff(
          possible_states,
          means[segment - 1, 2]
        ),
        size = 1
      )
    }
  }

  # Generate observations.
  x <- matrix(
    NA_real_,
    nrow = n,
    ncol = 2
  )

  for (segment in seq_len(n_segments)) {
    current_rows <- (
      boundaries[segment] + 1
    ):boundaries[segment + 1]

    x[current_rows, ] <- mvtnorm::rmvnorm(
      n = length(current_rows),
      mean = means[segment, ],
      sigma = Sigma
    )
  }

  colnames(x) <- c("Y1", "Y2")
  colnames(means) <- c("mu1", "mu2")

  list(
    x = x,
    means = means,
    changepoints = changepoints,
    lengths = diff(boundaries)
  )
}

# Data thinning function for 2 dim normal with known Sigma
thin_gaussian_2d <- function(
    x,
    Sigma,
    epsilon = 0.5
) {
  n <- nrow(x)

  thinning_covariance <-
    epsilon * (1 - epsilon) * Sigma

  z <- mvtnorm::rmvnorm(
    n = n,
    mean = c(0, 0),
    sigma = thinning_covariance
  )

  x1 <- epsilon * x + z
  x2 <- (1 - epsilon) * x - z

  colnames(x1) <- colnames(x)
  colnames(x2) <- colnames(x)

  list(
    x1 = x1,
    x2 = x2,
    epsilon = epsilon
  )
}

# get confidence regions based on chi2 quantile with Bonferroni correction
confidence_regions_thinning <- function(
    x2,
    breakpoints = numeric(),
    Sigma,
    epsilon = 0.5,
    alpha = 0.05,
    delta = 0.05
) {
  n <- nrow(x2)
  signal_dim <- ncol(x2)
  n_segments <- length(breakpoints) + 1

  segment_boundaries <- c(
    0,
    breakpoints,
    n
  )

  lengths <- diff(segment_boundaries)

  # ----------------------------------------------------------
  # Centers of the confidence regions for the ORIGINAL means
  # ----------------------------------------------------------

  means <- matrix(
    NA_real_,
    nrow = n_segments,
    ncol = signal_dim
  )

  colnames(means) <- colnames(x2)

  for (i in seq_len(n_segments)) {
    current_rows <- seq.int(
      from = segment_boundaries[i] + 1L,
      to = segment_boundaries[i + 1L]
    )

    means[i, ] <-
      colMeans(
        x2[current_rows, , drop = FALSE]
      ) /
      (1 - epsilon)
  }

  # ----------------------------------------------------------
  # Effective covariance for inference on the ORIGINAL means
  #
  # X2 / (1-epsilon) has covariance Sigma / (1-epsilon)
  # ----------------------------------------------------------

  covariance <- Sigma / (1 - epsilon)

  if (!is.null(colnames(x2))) {
    dimnames(covariance) <- list(
      colnames(x2),
      colnames(x2)
    )
  }

  # ----------------------------------------------------------
  # Bonferroni simultaneous quantile
  # ----------------------------------------------------------

  quantile <- sqrt(
    qchisq(
      p = 1 - alpha / n_segments,
      df = signal_dim
    )
  )

  # Use the same admissibility criterion as PoSIRd so that
  # comparisons are made on the same type of segments.
  minimum_length <- max(
    2L,
    ceiling(n * delta)
  )

  valid <- lengths >= minimum_length

  posir1DCR(
    x = x2,
    breakpoints = breakpoints,
    means = means,
    lengths = lengths,
    covariance = covariance,
    quantile = quantile,
    alpha = alpha,
    delta = delta,
    valid = valid
  )
}

get_true_targets <- function(
    true_signal,
    breakpoints
) {
  n <- nrow(true_signal)

  boundaries <- c(
    0,
    breakpoints,
    n
  )

  n_segments <- length(boundaries) - 1L

  targets <- matrix(
    NA_real_,
    nrow = n_segments,
    ncol = ncol(true_signal)
  )

  for (i in seq_len(n_segments)) {
    current_rows <- seq.int(
      from = boundaries[i] + 1L,
      to = boundaries[i + 1L]
    )

    targets[i, ] <- colMeans(
      true_signal[
        current_rows,
        ,
        drop = FALSE
      ]
    )
  }

  colnames(targets) <- colnames(true_signal)

  targets
}

confidence_region_areas <- function(region) {
  covariance_factor <- sqrt(
    det(region$covariance)
  )

  pi *
    region$quantile^2 /
    region$lengths *
    covariance_factor
}

make_true_signal <- function(n, changepoints, segment_means) {
  boundaries <- c(0L, changepoints, n)
  lengths <- diff(boundaries)
  signal <- segment_means[rep(seq_len(nrow(segment_means)), lengths), , drop = FALSE]
  colnames(signal) <- paste0("mu", seq_len(ncol(signal)))
  signal
}

replace_region_covariance <- function(regions, Sigma) {
  # confidence_regions_1D() supplies the PoSIRd quantile, means, lengths and
  # validity indicators. In the oracle experiment, only its estimated pooled
  # covariance is replaced by the known covariance.
  stopifnot(
    is.matrix(Sigma),
    nrow(Sigma) == ncol(Sigma),
    ncol(regions$data) == nrow(Sigma)
  )
  chol(Sigma) # fail early if Sigma is not positive definite
  regions$covariance <- Sigma
  regions
}

segment_coverage <- function(regions, targets) {
  vapply(
    seq_len(nrow(targets)),
    function(k) {
      is_in_confidence_region(
        point = targets[k, ],
        region = regions,
        segment = k
      )
    },
    logical(1)
  )
}

safe_summary <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(c(mean = NA_real_, median = NA_real_))
  }
  c(mean = mean(x), median = median(x))
}

summarise_method <- function(method, regions, targets, estimated_changepoints,
                             true_changepoints) {
  covered <- segment_coverage(regions, targets)
  areas <- confidence_region_areas(regions)
  valid <- regions$valid

  # With the same number of estimated and true changepoints, these two
  # diagnostics have simple interpretations. The Hausdorff distance also
  # remains well-defined if their numbers differ in future experiments.
  nearest_true <- vapply(
    estimated_changepoints,
    function(x) min(abs(x - true_changepoints)),
    numeric(1)
  )
  directed_est_to_true <- max(nearest_true)
  directed_true_to_est <- max(vapply(
    true_changepoints,
    function(x) min(abs(x - estimated_changepoints)),
    numeric(1)
  ))

  area_stats <- safe_summary(areas[valid])

  data.frame(
    method = method,
    simultaneous_coverage = if (any(valid)) all(covered[valid]) else NA,
    mean_area = unname(area_stats["mean"]),
    median_area = unname(area_stats["median"]),
    n_valid = sum(valid),
    proportion_valid = mean(valid),
    mean_cp_error = mean(nearest_true),
    hausdorff_cp_error = max(directed_est_to_true, directed_true_to_est),
    stringsAsFactors = FALSE
  )
}

confidence_regions_overlap <- function(region, i, j) {

  # Undefined confidence regions
  if (!region$valid[i] || !region$valid[j]) {
    return(NA)
  }

  # Squared Mahalanobis distance between centers
  distance_squared <- as.numeric(
    mahalanobis(
      x = matrix(region$means[i, ], nrow = 1),
      center = region$means[j, ],
      cov = region$covariance
    )
  )

  # Sum of the two Mahalanobis radii
  radius_sum <- region$quantile * (
    1 / sqrt(region$lengths[i]) +
      1 / sqrt(region$lengths[j])
  )

  sqrt(distance_squared) <= radius_sum
}

prune_breakpoints_from_posir <- function(posir_obj) {

  bp <- c(0, posir_obj$breakpoints)

  keep <- rep(TRUE, length(bp))

  for (i in seq_len(length(bp) - 1L)) {

    overlaps <- confidence_regions_overlap(
      posir_obj, i, i + 1L
    )

    if (isTRUE(overlaps)) {
      keep[i + 1L] <- FALSE
    }
  }

  bp[keep]
}
