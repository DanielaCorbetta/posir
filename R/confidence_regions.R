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
#'
#' @param Qtable A matrix containing the estimated POSIR quantiles for the
#'   appropriate signal dimension. Row names contain alpha values and column
#'   names contain delta values.
#'
#' @param covariance_method Method used to estimate the covariance matrix.
#'   One of `"pooled"`, `"diff_mom_entrywise"` or `"diff_mom_geometric"`.
#'   The pooled estimator uses residuals from the supplied segmentation.
#'   The other two methods use median-of-means estimators based on first
#'   differences.
#'
#' @param covariance_blocks Number of blocks used by the difference-based
#'   median-of-means estimators. If `NULL`, the number of blocks is set to
#'   `ceiling(sqrt(n - 1))`. Ignored when `covariance_method = "pooled"`.
#'
#' @param geometric_tolerance Convergence tolerance for the geometric-median
#'   algorithm.
#'
#' @param geometric_max_iterations Maximum number of iterations for the
#'   geometric-median algorithm.
#'
#' @param eigenvalue_floor Relative floor used to make a covariance estimate
#'   positive definite when necessary.
#'
#' @return An object of class `posir1DCR` containing the segment means,
#'   segment lengths, estimated covariance matrix and POSIR quantile.
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
#'   paste0(
#'     "Simulations/",
#'     "1D_Discretization_20000_Gaussian_",
#'     "grid_1to.005by.005_signal_dim_2/",
#'     "Table_quantiles.txt"
#'   )
#' )
#'
#' regions <- confidence_regions_1D(
#'   x = x,
#'   breakpoints = 50,
#'   delta = 0.1,
#'   alpha = 0.05,
#'   Qtable = Qtable,
#'   covariance_method = "diff_mom_geometric",
#'   covariance_blocks = 10
#' )
#' }
confidence_regions_1D <- function(
    x,
    breakpoints = numeric(),
    delta,
    alpha,
    Qtable = NULL,
    covariance_method = c(
      "pooled",
      "diff_mom_entrywise",
      "diff_mom_geometric"
    ),
    covariance_blocks = NULL,
    geometric_tolerance = 1e-8,
    geometric_max_iterations = 1000L,
    eigenvalue_floor = 1e-8
) {
  covariance_method <- match.arg(covariance_method)

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

  if (length(geometric_tolerance) != 1L ||
      !is.finite(geometric_tolerance) ||
      geometric_tolerance <= 0) {
    log_n_stop("geometric_tolerance must be a positive finite number.")
  }

  if (length(geometric_max_iterations) != 1L ||
      !is.finite(geometric_max_iterations) ||
      geometric_max_iterations < 1 ||
      geometric_max_iterations != as.integer(geometric_max_iterations)) {
    log_n_stop("geometric_max_iterations must be a positive integer.")
  }

  if (length(eigenvalue_floor) != 1L ||
      !is.finite(eigenvalue_floor) ||
      eigenvalue_floor <= 0) {
    log_n_stop("eigenvalue_floor must be a positive finite number.")
  }

  n <- nrow(x)
  signal_dim <- ncol(x)
  n_segments <- length(breakpoints) + 1L

  if (n_segments >= n) {
    log_n_stop("Too many breakpoints in confidence_regions_1D().")
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
    current_rows <- seq.int(
      from = segment_boundaries[i] + 1L,
      to = segment_boundaries[i + 1L]
    )

    means[i, ] <- colMeans(
      x[current_rows, , drop = FALSE]
    )
  }

  covariance <- switch(
    covariance_method,

    pooled = {
      estimate_covariance_pooled(
        x = x,
        means = means,
        lengths = lengths
      )
    },

    diff_mom_entrywise = {
      estimate_covariance_diff_mom(
        x = x,
        number_of_blocks = covariance_blocks,
        aggregation = "entrywise",
        eigenvalue_floor = eigenvalue_floor
      )
    },

    diff_mom_geometric = {
      estimate_covariance_diff_mom(
        x = x,
        number_of_blocks = covariance_blocks,
        aggregation = "geometric",
        geometric_tolerance = geometric_tolerance,
        geometric_max_iterations = geometric_max_iterations,
        eigenvalue_floor = eigenvalue_floor
      )
    }
  )

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
        "The covariance matrix estimated using",
        shQuote(covariance_method),
        "is not positive definite."
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

  if (length(quantile) != 1L ||
      !is.finite(quantile) ||
      quantile <= 0) {
    log_n_stop(
      "Invalid quantile obtained in confidence_regions_1D()."
    )
  }

  minimum_length <- max(2L, ceiling(n * delta))
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

#' Plot a two-dimensional segmented signal in three dimensions
#'
#' Plot observations from a two-dimensional signal against their temporal
#' index. Estimated segment means are represented as piecewise-constant
#' trajectories in the three-dimensional space defined by time, Y1 and Y2.
#'
#' @param x An object of class `posir1DCR` with signal dimension equal to 2.
#' @param segments Integer vector indicating the segments to plot. If `NULL`,
#'   all segments are plotted.
#' @param col Colours assigned to the segments. If `NULL`, colours are generated
#'   automatically.
#' @param show_data Whether to display the observed data points.
#' @param show_means Whether to display the estimated piecewise-constant mean
#'   trajectory.
#' @param show_jumps Whether to connect consecutive segment means at the
#'   breakpoints.
#' @param show_regions Whether to display the confidence-region boundary for
#'   each segment.
#' @param region_location Position along the temporal axis where confidence
#'   regions are drawn. One of `"middle"`, `"start"` or `"end"`.
#' @param region_npoints Number of points used to draw each confidence ellipse.
#' @param data_alpha Transparency of the observations.
#' @param data_size Size of the observations.
#' @param mean_width Width of the estimated mean trajectories.
#' @param jump_width Width of the lines representing changes between consecutive
#'   segment means.
#' @param region_width Width of confidence-region boundaries.
#' @param time_aspect Relative length of the temporal axis in the displayed
#'   three-dimensional scene.
#' @param xlab Label for the temporal axis.
#' @param ylab Label for the first signal component.
#' @param zlab Label for the second signal component.
#' @param main Plot title.
#'
#' @return A plotly object.
#'
#' @export
plot_posir1DCR_3d <- function(
    x,
    segments = NULL,
    col = NULL,
    show_data = TRUE,
    show_means = TRUE,
    show_jumps = TRUE,
    show_regions = FALSE,
    region_location = c("middle", "start", "end"),
    region_npoints = 150,
    data_alpha = 0.6,
    data_size = 2.5,
    mean_width = 5,
    jump_width = 3,
    region_width = 4,
    time_aspect = 2.2,
    xlab = "t",
    ylab = "Y1",
    zlab = "Y2",
    main = "PoSIRd segmented signal"
) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop(
      paste(
        "Package 'plotly' is required for plot_posir1DCR_3d().",
        "Install it with install.packages('plotly')."
      ),
      call. = FALSE
    )
  }

  if (!inherits(x, "posir1DCR") ||
      !validate_posir1DCR(x)) {
    log_n_stop(
      "Invalid object in plot_posir1DCR_3d()."
    )
  }

  if (ncol(x$data) != 2) {
    log_n_stop(
      paste(
        "plot_posir1DCR_3d() is available only",
        "for signal dimension 2."
      )
    )
  }

  region_location <- match.arg(region_location)

  n <- nrow(x$data)
  n_segments <- nrow(x$means)

  if (is.null(segments)) {
    segments <- seq_len(n_segments)
  }

  if (!is.numeric(segments) ||
      length(segments) == 0L ||
      any(!is.finite(segments)) ||
      any(segments != round(segments)) ||
      any(segments < 1L) ||
      any(segments > n_segments) ||
      anyDuplicated(segments)) {
    log_n_stop(
      "Invalid segment indices in plot_posir1DCR_3d()."
    )
  }

  segments <- as.integer(segments)

  logical_options <- c(
    show_data,
    show_means,
    show_jumps,
    show_regions
  )

  if (anyNA(logical_options) ||
      any(!vapply(
        list(show_data, show_means, show_jumps, show_regions),
        function(value) {
          is.logical(value) && length(value) == 1L
        },
        logical(1)
      ))) {
    log_n_stop(
      paste(
        "show_data, show_means, show_jumps and show_regions",
        "must each be a single non-missing logical value."
      )
    )
  }

  if (!is.numeric(region_npoints) ||
      length(region_npoints) != 1L ||
      !is.finite(region_npoints) ||
      region_npoints < 20L ||
      region_npoints != round(region_npoints)) {
    log_n_stop(
      "region_npoints must be an integer greater than or equal to 20."
    )
  }

  region_npoints <- as.integer(region_npoints)

  if (!is.numeric(data_alpha) ||
      length(data_alpha) != 1L ||
      !is.finite(data_alpha) ||
      data_alpha < 0 ||
      data_alpha > 1) {
    log_n_stop(
      "data_alpha must be a number between 0 and 1."
    )
  }

  positive_options <- c(
    data_size = data_size,
    mean_width = mean_width,
    jump_width = jump_width,
    region_width = region_width,
    time_aspect = time_aspect
  )

  if (any(!is.finite(positive_options)) ||
      any(positive_options <= 0)) {
    log_n_stop(
      paste(
        "data_size, mean_width, jump_width, region_width",
        "and time_aspect must be positive finite numbers."
      )
    )
  }

  if (is.null(col)) {
    segment_col <- grDevices::hcl.colors(
      n_segments,
      palette = "Dark 3"
    )
  } else {
    if (!is.character(col) || length(col) == 0L) {
      log_n_stop(
        "col must be a non-empty character vector."
      )
    }

    segment_col <- rep(
      col,
      length.out = n_segments
    )
  }

  boundaries <- c(
    0L,
    x$breakpoints,
    n
  )

  segment_start <- boundaries[-length(boundaries)] + 1L
  segment_end <- boundaries[-1L]
  segment_middle <- (segment_start + segment_end) / 2

  segment_id <- rep(
    seq_len(n_segments),
    times = x$lengths
  )

  selected_rows <- segment_id %in% segments

  point_col <- grDevices::adjustcolor(
    segment_col[segment_id],
    alpha.f = data_alpha
  )

  hover_text <- paste0(
    "t: ", seq_len(n),
    "<br>segment: ", segment_id,
    "<br>Y1: ", signif(x$data[, 1], 5),
    "<br>Y2: ", signif(x$data[, 2], 5)
  )

  p <- plotly::plot_ly()

  # Plot the observations.
  if (show_data) {
    p <- plotly::add_trace(
      p,
      x = seq_len(n)[selected_rows],
      y = x$data[selected_rows, 1],
      z = x$data[selected_rows, 2],
      type = "scatter3d",
      mode = "markers",
      marker = list(
        size = data_size,
        color = point_col[selected_rows]
      ),
      text = hover_text[selected_rows],
      hoverinfo = "text",
      name = "Observations",
      showlegend = FALSE
    )
  }

  # Plot the piecewise-constant estimated means.
  if (show_means) {
    for (segment in segments) {
      mean_hover <- paste0(
        "Segment: ", segment,
        "<br>t: ", segment_start[segment],
        "–", segment_end[segment],
        "<br>estimated mean Y1: ",
        signif(x$means[segment, 1], 5),
        "<br>estimated mean Y2: ",
        signif(x$means[segment, 2], 5)
      )

      p <- plotly::add_trace(
        p,
        x = c(
          segment_start[segment],
          segment_end[segment]
        ),
        y = rep(
          x$means[segment, 1],
          2L
        ),
        z = rep(
          x$means[segment, 2],
          2L
        ),
        type = "scatter3d",
        mode = "lines",
        line = list(
          color = segment_col[segment],
          width = mean_width
        ),
        text = rep(mean_hover, 2L),
        hoverinfo = "text",
        name = paste("Segment", segment),
        showlegend = TRUE
      )

      # Mark the estimated mean at the midpoint of the segment.
      # p <- plotly::add_trace(
      #   p,
      #   x = segment_middle[segment],
      #   y = x$means[segment, 1],
      #   z = x$means[segment, 2],
      #   type = "scatter3d",
      #   mode = "markers",
      #   marker = list(
      #     size = data_size + 3,
      #     color = segment_col[segment],
      #     symbol = "x"
      #   ),
      #   text = mean_hover,
      #   hoverinfo = "text",
      #   showlegend = FALSE
      # )
    }
  }

  # Connect means at consecutive breakpoints.
  if (show_jumps && length(segments) >= 2L) {
    consecutive_pairs <- which(
      segments[-1L] ==
        segments[-length(segments)] + 1L
    )

    for (j in consecutive_pairs) {
      segment_left <- segments[j]
      segment_right <- segments[j + 1L]
      breakpoint <- x$breakpoints[segment_left]

      p <- plotly::add_trace(
        p,
        x = rep(breakpoint + 0.5, 2L),
        y = c(
          x$means[segment_left, 1],
          x$means[segment_right, 1]
        ),
        z = c(
          x$means[segment_left, 2],
          x$means[segment_right, 2]
        ),
        type = "scatter3d",
        mode = "lines",
        line = list(
          color = "#222222",
          width = jump_width,
          dash = "dash"
        ),
        text = rep(
          paste0(
            "Breakpoint: ", breakpoint,
            "<br>segment ", segment_left,
            " to segment ", segment_right
          ),
          2L
        ),
        hoverinfo = "text",
        name = "Mean changes",
        showlegend = j == consecutive_pairs[1L]
      )
    }
  }

  # Plot confidence-region boundaries.
  if (show_regions) {
    invalid_segments <- segments[!x$valid[segments]]

    if (length(invalid_segments) > 0L) {
      log_warn(
        paste(
          "Confidence regions were not drawn for invalid segments:",
          paste(invalid_segments, collapse = ", ")
        )
      )
    }

    theta <- seq(
      0,
      2 * pi,
      length.out = region_npoints
    )

    unit_circle <- cbind(
      cos(theta),
      sin(theta)
    )

    covariance_root <- chol(x$covariance)

    for (segment in segments[x$valid[segments]]) {
      radius_multiplier <- x$quantile /
        sqrt(x$lengths[segment])

      ellipse <- sweep(
        radius_multiplier *
          unit_circle %*% covariance_root,
        MARGIN = 2L,
        STATS = x$means[segment, ],
        FUN = "+"
      )

      region_t <- switch(
        region_location,
        middle = segment_middle[segment],
        start = segment_start[segment],
        end = segment_end[segment]
      )

      region_hover <- paste0(
        "Confidence region for segment ", segment,
        "<br>displayed at t = ", signif(region_t, 5),
        "<br>segment length: ", x$lengths[segment]
      )

      p <- plotly::add_trace(
        p,
        x = rep(region_t, nrow(ellipse)),
        y = ellipse[, 1],
        z = ellipse[, 2],
        type = "scatter3d",
        mode = "lines",
        line = list(
          color = segment_col[segment],
          width = region_width
        ),
        text = rep(region_hover, nrow(ellipse)),
        hoverinfo = "text",
        name = paste("CR segment", segment),
        showlegend = FALSE
      )
    }
  }

  # Use a manual aspect ratio because the range of t is typically much
  # larger than the ranges of Y1 and Y2.
  p <- plotly::layout(
    p,
    title = list(
      text = main
    ),
    scene = list(
      xaxis = list(
        title = list(text = xlab)
      ),
      yaxis = list(
        title = list(text = ylab)
      ),
      zaxis = list(
        title = list(text = zlab)
      ),
      aspectmode = "manual",
      aspectratio = list(
        x = time_aspect,
        y = 1,
        z = 1
      ),
      camera = list(
        eye = list(
          x = 1.5,
          y = 1.5,
          z = 1.2
        )
      )
    ),
    legend = list(
      title = list(
        text = "Segments"
      )
    )
  )

  p
}

#' Plot projected confidence intervals for one signal component
#'
#' Plot one component of a multivariate signal over time, together with the
#' estimated segment means and the coordinate-wise projections of the POSIR
#' confidence regions.
#'
#' @param x A `posir1DCR` object.
#' @param component Component of the signal to plot.
#' @param box A boolean value. If `TRUE`, the projected confidence interval is
#'   displayed as a semi-transparent box. Otherwise, it is displayed as a pair
#'   of dashed lines.
#' @inheritParams base::plot
#' @seealso [posir1DCR()]
#'
#' @return The object `x`, invisibly.
#'
#' @export
#' @importFrom graphics abline lines plot rect
#' @importFrom grDevices adjustcolor
#'
#' @examples
#' \dontrun{
#' plot_posir1DCR_component(
#'   regions,
#'   component = 1,
#'   box = TRUE
#' )
#' }
plot_posir1DCR_component <- function(
    x,
    component = 1,
    box = FALSE,
    ...
) {
  if (!inherits(x, "posir1DCR") ||
      !validate_posir1DCR(x)) {
    log_n_stop(
      "Invalid arguments in plot_posir1DCR_component()."
    )
  }

  signal_dim <- ncol(x$data)

  if (!is.numeric(component) ||
      length(component) != 1 ||
      !is.finite(component) ||
      component != round(component) ||
      component < 1 ||
      component > signal_dim) {
    log_n_stop(
      "component must be a valid signal-component index."
    )
  }

  if (!is.logical(box) ||
      length(box) != 1 ||
      is.na(box)) {
    log_n_stop(
      "box must be a single non-missing logical value."
    )
  }

  component <- as.integer(component)

  n <- nrow(x$data)
  n_segments <- length(x$breakpoints) + 1
  aux_bkp <- c(0, x$breakpoints, n)

  plotmeans <- !is.null(x$means)

  component_name <- colnames(x$data)[component]

  if (is.null(component_name) ||
      is.na(component_name) ||
      component_name == "") {
    component_name <- paste0("Y", component)
  }

  graphics::plot(
    seq_len(n),
    x$data[, component],
    main = paste(
      "Simultaneous confidence intervals for",
      component_name
    ),
    xlab = "t",
    ylab = component_name,
    col = "#33333333",
    pch = 19,
    cex = 0.5,
    ...
  )

  CI_col <- "#FF0000"
  CI_col_alpha <- grDevices::adjustcolor(
    CI_col,
    alpha.f = 0.3
  )

  if (length(x$breakpoints) > 0) {
    for (i in seq_along(x$breakpoints)) {
      graphics::abline(
        v = x$breakpoints[i] + 0.5,
        col = CI_col
      )
    }
  }

  for (i in seq_len(n_segments)) {
    xx <- c(
      aux_bkp[i] + 0.5,
      aux_bkp[i + 1] + 0.5
    )

    if (plotmeans) {
      graphics::lines(
        xx,
        rep(x$means[i, component], 2),
        col = CI_col
      )
    }

    if (x$valid[i]) {
      half_width <- x$quantile *
        sqrt(
          x$covariance[component, component] /
            x$lengths[i]
        )

      lower_bound <- x$means[i, component] -
        half_width

      upper_bound <- x$means[i, component] +
        half_width

      if (box) {
        graphics::rect(
          xx[1],
          lower_bound,
          xx[2],
          upper_bound,
          border = CI_col,
          col = CI_col_alpha
        )
      } else {
        graphics::lines(
          xx,
          rep(lower_bound, 2),
          col = CI_col,
          lty = 2
        )

        graphics::lines(
          xx,
          rep(upper_bound, 2),
          col = CI_col,
          lty = 2
        )
      }
    }
  }

  invisible(x)
}
