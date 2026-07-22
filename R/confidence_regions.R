#' Check data arguments for one-dimensional POSIR confidence regions
#'
#' @param x A numeric matrix whose rows are observations and whose columns
#'   are signal dimensions.
#' @param breakpoints A numeric vector containing breakpoint positions,
#'   given in increasing order. Each breakpoint is the last row of a segment.
#'
#' @return A boolean.
#'
#' @keywords internal
CR1Ddatavalidation <- function(x, breakpoints) {
  if (!is.matrix(x) ||
      !is.numeric(x) ||
      nrow(x) < 1 ||
      ncol(x) < 2 ||
      any(!is.finite(x))) {
    return(FALSE)
  }

  if (!is.null(breakpoints) &&
      (!is.numeric(breakpoints) || !is.vector(breakpoints))) {
    return(FALSE)
  }

  n <- nrow(x)
  aux_bkp <- c(0, breakpoints, n)

  all(
    aux_bkp == round(aux_bkp),
    diff(aux_bkp) > 0
  )
}


#' Check parameters for one-dimensional POSIR confidence regions
#'
#' @param delta Lower bound on the relative segment length.
#' @param alpha Simultaneous error level.
#'
#' @return A boolean.
#'
#' @keywords internal
CR1Dparamvalidation <- function(delta, alpha) {
  is.numeric(delta) &&
    is.vector(delta) &&
    length(delta) == 1 &&
    is.finite(delta) &&
    delta >= 0.005 &&
    delta <= 1 &&
    is.numeric(alpha) &&
    is.vector(alpha) &&
    length(alpha) == 1 &&
    is.finite(alpha) &&
    alpha >= 0.001 &&
    alpha <= 0.5
}

#' Validate an object of class posir1DCR
#'
#' @param y An object of class posir1DCR.
#'
#' @return A boolean.
#' @export
validate_posir1DCR <- function(y) {
  if (!is.list(y) ||
      !CR1Ddatavalidation(y$data, y$breakpoints)) {
    return(FALSE)
  }

  n <- nrow(y$data)
  signal_dim <- ncol(y$data)
  n_segments <- length(y$breakpoints) + 1

  if (!is.matrix(y$means) ||
      !is.numeric(y$means) ||
      !all(dim(y$means) == c(n_segments, signal_dim)) ||
      any(!is.finite(y$means))) {
    return(FALSE)
  }

  expected_lengths <- diff(c(0, y$breakpoints, n))

  if (!is.numeric(y$lengths) ||
      length(y$lengths) != n_segments ||
      !isTRUE(all.equal(
        as.numeric(y$lengths),
        as.numeric(expected_lengths)
      ))) {
    return(FALSE)
  }

  if (!is.matrix(y$covariance) ||
      !is.numeric(y$covariance) ||
      !all(dim(y$covariance) == c(signal_dim, signal_dim)) ||
      any(!is.finite(y$covariance)) ||
      !isSymmetric(y$covariance)) {
    return(FALSE)
  }

  covariance_chol <- try(
    chol(y$covariance),
    silent = TRUE
  )

  if (inherits(covariance_chol, "try-error")) {
    return(FALSE)
  }

  if (!is.numeric(y$quantile) ||
      length(y$quantile) != 1 ||
      !is.finite(y$quantile) ||
      y$quantile <= 0) {
    return(FALSE)
  }

  if (!CR1Dparamvalidation(y$delta, y$alpha)) {
    return(FALSE)
  }

  if (!is.logical(y$valid) ||
      length(y$valid) != n_segments ||
      anyNA(y$valid)) {
    return(FALSE)
  }

  TRUE
}

#' Create an object of class posir1DCR
#'
#' @param x A numeric matrix of observations.
#' @param breakpoints Breakpoint positions.
#' @param means Matrix of segment means.
#' @param lengths Vector of segment lengths.
#' @param covariance Global pooled covariance matrix.
#' @param quantile POSIR quantile.
#' @param alpha Simultaneous error level.
#' @param delta Lower bound on relative segment length.
#' @param valid Logical vector indicating which segments are sufficiently long.
#'
#' @return An object of class posir1DCR.
#' @export
posir1DCR <- function(x, breakpoints, means, lengths,
                      covariance, quantile, alpha, delta, valid) {
  result <- list(
    data = x,
    breakpoints = breakpoints,
    means = means,
    lengths = lengths,
    covariance = covariance,
    quantile = quantile,
    alpha = alpha,
    delta = delta,
    valid = valid
  )

  class(result) <- "posir1DCR"

  if (!validate_posir1DCR(result)) {
    log_n_stop("Invalid arguments in posir1DCR().")
  }

  result
}

#' Compute simultaneous confidence regions for multivariate segment means
#'
#' Compute POSIR confidence regions for the mean vectors associated with
#' a segmentation of a multivariate signal indexed by one dimension.
#'
#' @inheritParams CR1Ddatavalidation
#' @inheritParams CR1Dparamvalidation
#' @param Qtable A matrix containing the estimated POSIR quantiles for the
#'   appropriate signal dimension. Row names contain alpha values and column
#'   names contain delta values.
#'
#' @return An object of class posir1DCR containing the segment means,
#'   segment lengths, global pooled covariance matrix and POSIR quantile.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' x <- cbind(
#'   rnorm(100) + rep(c(0, 2), each = 50),
#'   rnorm(100) + rep(c(0, -1), each = 50)
#' )
#'
#' Qtable <- read_quantiles(
#'   "Simulations/1D_Discretization_20000_Gaussian_grid_1to.005by.005_signal_dim_2/Table_quantiles.txt"
#' )
#'
#' regions <- confidence_regions_1D(
#'   x = x,
#'   breakpoints = 50,
#'   delta = 0.1,
#'   alpha = 0.05,
#'   Qtable = Qtable
#' )
#' }
confidence_regions_1D <- function(x, breakpoints = numeric(),
                                  delta, alpha, Qtable = NULL) {
  if (!CR1Ddatavalidation(x, breakpoints) ||
      !CR1Dparamvalidation(delta, alpha)) {
    log_n_stop("Invalid argument in confidence_regions_1D().")
  }

  if (is.null(Qtable)) {
    log_n_stop(
      paste(
        "Qtable must be supplied in confidence_regions_1D()",
        "for a multidimensional signal."
      )
    )
  }

  n <- nrow(x)
  signal_dim <- ncol(x)
  n_segments <- length(breakpoints) + 1

  if (n_segments >= n) {
    log_n_stop("Too many breakpoints in confidence_regions_1D().")
  }

  residual_df <- n - n_segments

  if (residual_df < signal_dim) {
    log_n_stop(
      paste(
        "Not enough residual degrees of freedom to estimate",
        "a positive-definite covariance matrix."
      )
    )
  }

  segment_boundaries <- c(0, breakpoints, n)
  lengths <- diff(segment_boundaries)

  means <- matrix(
    NA_real_,
    nrow = n_segments,
    ncol = signal_dim
  )

  colnames(means) <- colnames(x)

  for (i in seq_len(n_segments)) {
    current_rows <- (segment_boundaries[i] + 1):segment_boundaries[i + 1]

    means[i, ] <- colMeans(
      x[current_rows, , drop = FALSE]
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

  covariance <- crossprod(residuals) / residual_df

  if (!is.null(colnames(x))) {
    dimnames(covariance) <- list(
      colnames(x),
      colnames(x)
    )
  }

  covariance_chol <- try(
    chol(covariance),
    silent = TRUE
  )

  if (inherits(covariance_chol, "try-error")) {
    log_n_stop(
      paste(
        "The pooled covariance matrix estimated in",
        "confidence_regions_1D() is not positive definite."
      )
    )
  }

  quantile <- as.numeric(
    posir_quantiles(
      alpha = alpha,
      delta = delta,
      d = 1,
      Qtable = Qtable
    )
  )

  if (length(quantile) != 1 ||
      !is.finite(quantile) ||
      quantile <= 0) {
    log_n_stop(
      "Invalid quantile obtained in confidence_regions_1D()."
    )
  }

  minimum_length <- max(2, ceiling(n * delta))
  valid <- lengths >= minimum_length

  if (any(!valid)) {
    log_warn(
      paste(
        "In confidence_regions_1D(): some segments were not long enough,",
        "so the corresponding confidence regions could not be computed."
      )
    )
  }

  posir1DCR(
    x = x,
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

#' Test membership in a POSIR confidence region
#'
#' Test whether a candidate mean vector belongs to the simultaneous
#' confidence region associated with a selected segment.
#'
#' @param point A numeric vector representing the candidate mean.
#' @param region An object of class posir1DCR.
#' @param segment Index of the segment whose confidence region is considered.
#'
#' @return `TRUE` if `point` belongs to the confidence region, `FALSE` if it
#'   does not, and `NA` if the segment is too short for its confidence region
#'   to be defined.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' is_in_confidence_region(
#'   point = c(0, 0),
#'   region = regions,
#'   segment = 1
#' )
#' }
is_in_confidence_region <- function(point, region, segment) {
  if (!inherits(region, "posir1DCR") ||
      !validate_posir1DCR(region)) {
    log_n_stop(
      "Invalid region object in is_in_confidence_region()."
    )
  }

  signal_dim <- ncol(region$data)
  n_segments <- nrow(region$means)

  if (!is.numeric(point) ||
      !is.vector(point) ||
      length(point) != signal_dim ||
      any(!is.finite(point))) {
    log_n_stop(
      paste(
        "In is_in_confidence_region(), point must be a finite numeric",
        "vector with one entry per signal dimension."
      )
    )
  }

  if (!is.numeric(segment) ||
      length(segment) != 1 ||
      !is.finite(segment) ||
      segment != round(segment) ||
      segment < 1 ||
      segment > n_segments) {
    log_n_stop(
      "Invalid segment index in is_in_confidence_region()."
    )
  }

  segment <- as.integer(segment)

  if (!region$valid[segment]) {
    return(NA)
  }

  distance_squared <- as.numeric(
    mahalanobis(
      x = matrix(point, nrow = 1),
      center = region$means[segment, ],
      cov = region$covariance
    )
  )

  statistic <- region$lengths[segment] * distance_squared

  statistic <= region$quantile^2
}

#' Plot two-dimensional POSIR confidence regions
#'
#' Plot the observations, estimated segment means, and boundaries of the
#' simultaneous confidence regions for a two-dimensional signal.
#'
#' @param x An object of class posir1DCR with signal dimension equal to 2.
#' @param mu1_grid Increasing numeric vector defining the grid for the first
#'   signal component.
#' @param mu2_grid Increasing numeric vector defining the grid for the second
#'   signal component.
#' @param segments Integer vector indicating the segments to plot. If `NULL`,
#'   all segments are plotted.
#' @param col Colours assigned to the segments. The same colours are used for
#'   observations, confidence-region boundaries, and estimated segment means.
#' @param lwd Line width used for the confidence-region boundaries.
#' @param show_data Whether to display the observed data points.
#' @param color_data_by_segment Whether to colour observations according to
#'   their segment.
#' @param data_col Colour used for observations when
#'   `color_data_by_segment = FALSE`.
#' @param data_alpha Transparency of segment-coloured observations.
#' @param data_pch Plotting symbol used for the observations.
#' @param data_cex Size of the observations.
#' @param show_centers Whether to display the estimated segment means.
#' @param show_labels Whether to label the estimated segment means.
#' @param xlab Label for the horizontal axis.
#' @param ylab Label for the vertical axis.
#' @param asp Aspect ratio of the plot. The default value of 1 uses the same
#'   scale on both axes.
#' @param ... Further arguments passed to [graphics::plot()].
#'
#' @return The object `x`, invisibly.
#'
#' @export
#' @importFrom graphics contour plot points text
#' @importFrom grDevices adjustcolor hcl.colors
#'
#' @examples
#' \dontrun{
#' mu1_grid <- seq(-1, 3, length.out = 300)
#' mu2_grid <- seq(-2, 1, length.out = 300)
#'
#' plot(
#'   regions,
#'   mu1_grid = mu1_grid,
#'   mu2_grid = mu2_grid,
#'   main = "Simultaneous confidence regions"
#' )
#' }
plot.posir1DCR <- function(
    x,
    mu1_grid,
    mu2_grid,
    segments = NULL,
    col = NULL,
    lwd = 2,
    show_data = TRUE,
    color_data_by_segment = TRUE,
    data_col = "#33333340",
    data_alpha = 0.35,
    data_pch = 19,
    data_cex = 0.6,
    show_centers = TRUE,
    show_labels = TRUE,
    xlab = expression(Y[1]),
    ylab = expression(Y[2]),
    asp = 1,
    ...
) {
  if (!inherits(x, "posir1DCR") ||
      !validate_posir1DCR(x)) {
    log_n_stop("Invalid object in plot.posir1DCR().")
  }

  if (ncol(x$data) != 2) {
    log_n_stop(
      paste(
        "plot.posir1DCR() is currently available only",
        "for signal dimension 2."
      )
    )
  }

  valid_grid <- function(grid) {
    is.numeric(grid) &&
      is.vector(grid) &&
      length(grid) >= 2 &&
      all(is.finite(grid)) &&
      all(diff(grid) > 0)
  }

  if (!valid_grid(mu1_grid) ||
      !valid_grid(mu2_grid)) {
    log_n_stop(
      paste(
        "mu1_grid and mu2_grid must be finite, strictly increasing",
        "numeric vectors containing at least two values."
      )
    )
  }

  n_segments <- nrow(x$means)

  if (is.null(segments)) {
    segments <- seq_len(n_segments)
  }

  if (!is.numeric(segments) ||
      length(segments) == 0 ||
      any(!is.finite(segments)) ||
      any(segments != round(segments)) ||
      any(segments < 1) ||
      any(segments > n_segments) ||
      anyDuplicated(segments)) {
    log_n_stop(
      "Invalid segment indices in plot.posir1DCR()."
    )
  }

  segments <- as.integer(segments)
  plotted_segments <- segments[x$valid[segments]]

  if (length(plotted_segments) == 0) {
    log_n_stop(
      "None of the selected segments has a valid confidence region."
    )
  }

  if (length(plotted_segments) < length(segments)) {
    log_warn(
      paste(
        "In plot.posir1DCR(): confidence regions for segments shorter",
        "than the minimum admissible length were not plotted."
      )
    )
  }

  # Assign one colour to every segment.
  if (is.null(col)) {
    segment_col <- hcl.colors(
      n_segments,
      palette = "Dark 3"
    )
  } else {
    segment_col <- rep(
      col,
      length.out = n_segments
    )
  }

  plotted_col <- segment_col[plotted_segments]

  # Segment membership of each observation.
  segment_id <- rep(
    seq_len(n_segments),
    times = x$lengths
  )

  # Cartesian grid on which the quadratic statistics are evaluated.
  grid_points <- as.matrix(
    expand.grid(
      mu1 = mu1_grid,
      mu2 = mu2_grid
    )
  )

  plot(
    NA,
    xlim = range(mu1_grid),
    ylim = range(mu2_grid),
    xlab = xlab,
    ylab = ylab,
    asp = asp,
    ...
  )

  # Plot observations.
  if (show_data) {
    if (color_data_by_segment) {
      point_col <- adjustcolor(
        segment_col[segment_id],
        alpha.f = data_alpha
      )
    } else {
      point_col <- data_col
    }

    points(
      x$data[, 1],
      x$data[, 2],
      col = point_col,
      pch = data_pch,
      cex = data_cex
    )
  }

  boundary_level <- x$quantile^2

  # Plot confidence-region boundaries.
  for (j in seq_along(plotted_segments)) {
    segment <- plotted_segments[j]

    distance_squared <- mahalanobis(
      x = grid_points,
      center = x$means[segment, ],
      cov = x$covariance
    )

    statistic <- x$lengths[segment] *
      distance_squared

    statistic_matrix <- matrix(
      statistic,
      nrow = length(mu1_grid),
      ncol = length(mu2_grid)
    )

    if (boundary_level >= min(statistic_matrix) &&
        boundary_level <= max(statistic_matrix)) {
      contour(
        x = mu1_grid,
        y = mu2_grid,
        z = statistic_matrix,
        levels = boundary_level,
        add = TRUE,
        drawlabels = FALSE,
        col = plotted_col[j],
        lwd = lwd
      )
    } else {
      log_warn(
        paste(
          "The confidence-region boundary for segment",
          segment,
          "falls outside the supplied plotting grid."
        )
      )
    }
  }

  # Plot estimated segment means.
  if (show_centers) {
    points(
      x$means[plotted_segments, 1],
      x$means[plotted_segments, 2],
      col = plotted_col,
      pch = 19,
      cex = 1.2
    )
  }

  # Label estimated segment means.
  if (show_labels) {
    text(
      x$means[plotted_segments, 1],
      x$means[plotted_segments, 2],
      labels = plotted_segments,
      col = plotted_col,
      pos = 3
    )
  }

  invisible(x)
}
