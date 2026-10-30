#########################################################
### 27-08-26
### posird vs data thinning/data splitting
# Get Qtable
Qtable_path <- here::here(
  "Simulations",
  "Table_quantiles.txt"
)

Qtable <- read_quantiles(
  Qtable_path
)

source("inst/scripts/posird/utils_competitors.R")


# ============================================================
# Illustrative example:
# PoSIRd vs correct thinning vs misspecified thinning
# ============================================================

set.seed(1535)

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

n <- 1000
n_chg <- 10
delta <- 0.05
alpha <- 0.05
epsilon <- 0.5
a <- 2

rho <- 0.7

Sigma_true <- matrix(
  c(
    1, rho,
    rho, 1
  ),
  nrow = 2,
  byrow = TRUE
)

# Wrong covariance assumed by the thinning procedure
Sigma_wrong <- diag(2)

ps <- c(-a, 0, a)


# ------------------------------------------------------------
# Generate changepoints
# ------------------------------------------------------------

changepoint_data <- generate_changepoints(
  n = n,
  n_chg = n_chg,
  delta = delta
)

true_chg <- changepoint_data$changepoints


# ------------------------------------------------------------
# Generate Gaussian piecewise-constant data
# ------------------------------------------------------------

generated_data <- generate_piecewise_normal_2d(
  n = n,
  changepoints = true_chg,
  Sigma = Sigma_true,
  possible_states = ps
)

x <- generated_data$x
true_means <- generated_data$means


# ------------------------------------------------------------
# Correct data thinning
# ------------------------------------------------------------

thinned_correct <- thin_gaussian_2d(
  x = x,
  Sigma = Sigma_true,
  epsilon = epsilon
)

x1_correct <- thinned_correct$x1
x2_correct <- thinned_correct$x2


# ------------------------------------------------------------
# Misspecified data thinning
#
# The folds are generated using Sigma_wrong instead of Sigma_true.
# ------------------------------------------------------------

thinned_wrong <- thin_gaussian_2d(
  x = x,
  Sigma = Sigma_wrong,
  epsilon = epsilon
)

x1_wrong <- thinned_wrong$x1
x2_wrong <- thinned_wrong$x2


# ------------------------------------------------------------
# Estimate changepoints
#
# Oracle K = true number of changepoints
# ------------------------------------------------------------

est_chg_p <- sort(
  jointseg::doDynamicProgramming(
    x,
    K = n_chg
  )$bkp
)

est_chg_th_correct <- sort(
  jointseg::doDynamicProgramming(
    x1_correct,
    K = n_chg
  )$bkp
)

est_chg_th_wrong <- sort(
  jointseg::doDynamicProgramming(
    x1_wrong,
    K = n_chg
  )$bkp
)


# ------------------------------------------------------------
# Construct confidence regions
# ------------------------------------------------------------

# PoSIRd
regions_p <- confidence_regions_1D(
  x = x,
  breakpoints = est_chg_p,
  delta = delta,
  alpha = alpha,
  Qtable = Qtable,
  covariance_method = "pooled"
)


# Correctly specified thinning
regions_th_correct <- confidence_regions_thinning(
  x2 = x2_correct,
  breakpoints = est_chg_th_correct,
  Sigma = Sigma_true,
  epsilon = epsilon,
  alpha = alpha,
  delta = delta
)


# Misspecified thinning
#
# The procedure believes Sigma_wrong is the covariance,
# so it uses Sigma_wrong also in the confidence-region construction.
regions_th_wrong <- confidence_regions_thinning(
  x2 = x2_wrong,
  breakpoints = est_chg_th_wrong,
  Sigma = Sigma_wrong,
  epsilon = epsilon,
  alpha = alpha,
  delta = delta
)


# ------------------------------------------------------------
# Construct the true signal at every observation
# ------------------------------------------------------------

true_boundaries <- c(
  0,
  true_chg,
  n
)

n_true_segments <- length(true_chg) + 1

true_signal <- matrix(
  NA_real_,
  nrow = n,
  ncol = 2
)

for (i in seq_len(n_true_segments)) {
  current_rows <- seq.int(
    from = true_boundaries[i] + 1L,
    to = true_boundaries[i + 1L]
  )

  true_signal[current_rows, ] <- matrix(
    true_means[i, ],
    nrow = length(current_rows),
    ncol = 2,
    byrow = TRUE
  )
}

colnames(true_signal) <- c("mu1", "mu2")


# ------------------------------------------------------------
# True targets for the three segmentations
# ------------------------------------------------------------

targets_p <- get_true_targets(
  true_signal = true_signal,
  breakpoints = est_chg_p
)

targets_th_correct <- get_true_targets(
  true_signal = true_signal,
  breakpoints = est_chg_th_correct
)

targets_th_wrong <- get_true_targets(
  true_signal = true_signal,
  breakpoints = est_chg_th_wrong
)


# ------------------------------------------------------------
# Segmentwise coverage
# ------------------------------------------------------------

coverage_p <- vapply(
  seq_len(nrow(targets_p)),
  function(i) {
    is_in_confidence_region(
      point = targets_p[i, ],
      region = regions_p,
      segment = i
    )
  },
  logical(1)
)

coverage_th_correct <- vapply(
  seq_len(nrow(targets_th_correct)),
  function(i) {
    is_in_confidence_region(
      point = targets_th_correct[i, ],
      region = regions_th_correct,
      segment = i
    )
  },
  logical(1)
)

coverage_th_wrong <- vapply(
  seq_len(nrow(targets_th_wrong)),
  function(i) {
    is_in_confidence_region(
      point = targets_th_wrong[i, ],
      region = regions_th_wrong,
      segment = i
    )
  },
  logical(1)
)


# ------------------------------------------------------------
# Simultaneous coverage
# ------------------------------------------------------------

simultaneous_coverage <- c(
  PoSIRd = all(
    coverage_p[regions_p$valid]
  ),
  Thinning_correct = all(
    coverage_th_correct[
      regions_th_correct$valid
    ]
  ),
  Thinning_misspecified = all(
    coverage_th_wrong[
      regions_th_wrong$valid
    ]
  )
)

simultaneous_coverage


# ------------------------------------------------------------
# Areas
# ------------------------------------------------------------

areas_p <- confidence_region_areas(
  regions_p
)

areas_th_correct <- confidence_region_areas(
  regions_th_correct
)

areas_th_wrong <- confidence_region_areas(
  regions_th_wrong
)


# ------------------------------------------------------------
# Area summary
# ------------------------------------------------------------

area_summary <- data.frame(
  method = c(
    "PoSIRd",
    "Thinning: correct Sigma",
    "Thinning: misspecified Sigma"
  ),
  mean_area = c(
    mean(
      areas_p[
        regions_p$valid
      ]
    ),
    mean(
      areas_th_correct[
        regions_th_correct$valid
      ]
    ),
    mean(
      areas_th_wrong[
        regions_th_wrong$valid
      ]
    )
  ),
  median_area = c(
    median(
      areas_p[
        regions_p$valid
      ]
    ),
    median(
      areas_th_correct[
        regions_th_correct$valid
      ]
    ),
    median(
      areas_th_wrong[
        regions_th_wrong$valid
      ]
    )
  )
)

area_summary


# ------------------------------------------------------------
# Segmentwise results
# ------------------------------------------------------------

results_p <- data.frame(
  segment = seq_len(nrow(targets_p)),
  length = regions_p$lengths,
  target_mu1 = targets_p[, 1],
  target_mu2 = targets_p[, 2],
  estimated_mu1 = regions_p$means[, 1],
  estimated_mu2 = regions_p$means[, 2],
  valid = regions_p$valid,
  covered = coverage_p,
  area = areas_p
)

results_th_correct <- data.frame(
  segment = seq_len(nrow(targets_th_correct)),
  length = regions_th_correct$lengths,
  target_mu1 = targets_th_correct[, 1],
  target_mu2 = targets_th_correct[, 2],
  estimated_mu1 = regions_th_correct$means[, 1],
  estimated_mu2 = regions_th_correct$means[, 2],
  valid = regions_th_correct$valid,
  covered = coverage_th_correct,
  area = areas_th_correct
)

results_th_wrong <- data.frame(
  segment = seq_len(nrow(targets_th_wrong)),
  length = regions_th_wrong$lengths,
  target_mu1 = targets_th_wrong[, 1],
  target_mu2 = targets_th_wrong[, 2],
  estimated_mu1 = regions_th_wrong$means[, 1],
  estimated_mu2 = regions_th_wrong$means[, 2],
  valid = regions_th_wrong$valid,
  covered = coverage_th_wrong,
  area = areas_th_wrong
)

results_p
results_th_correct
results_th_wrong


# ------------------------------------------------------------
# Plot the three sets of confidence regions
# ------------------------------------------------------------

mu1_grid <- seq(-3, 3, length.out = 300)
mu2_grid <- seq(-3, 3, length.out = 300)

old_par <- par(
  mfrow = c(1, 3),
  mar = c(4, 4, 3, 1)
)

plot(
  regions_p,
  mu1_grid = mu1_grid,
  mu2_grid = mu2_grid,
  main = "PoSIRd",
  show_labels = FALSE
)

points(
  targets_p[, 1],
  targets_p[, 2],
  pch = 4,
  col = "black",
  lwd = 3,
  cex = 1.2
)


plot(
  regions_th_correct,
  mu1_grid = mu1_grid,
  mu2_grid = mu2_grid,
  main = "Thinning: correct Sigma",
  show_labels = FALSE
)

points(
  targets_th_correct[, 1],
  targets_th_correct[, 2],
  pch = 4,
  col = "black",
  lwd = 3,
  cex = 1.2
)


plot(
  regions_th_wrong,
  mu1_grid = mu1_grid,
  mu2_grid = mu2_grid,
  main = "Thinning: wrong Sigma",
  show_labels = FALSE
)

points(
  targets_th_wrong[, 1],
  targets_th_wrong[, 2],
  pch = 4,
  col = "black",
  lwd = 3,
  cex = 1.2
)

par(old_par)
#-----------------------------------------------------------
# MC sims
#-----------------------------------------------------------
run_one_simulation <- function(
    n = 1000,
    n_chg = 10,
    delta = 0.05,
    alpha = 0.05,
    epsilon = 0.5,
    a = 2,
    rho = 0.7,
    Sigma_wrong = diag(2),
    Qtable
) {
  # ------------------------------------------------------------
  # Covariance and possible mean states
  # ------------------------------------------------------------

  Sigma_true <- matrix(
    c(
      1, rho,
      rho, 1
    ),
    nrow = 2,
    byrow = TRUE
  )

  ps <- c(-a, 0, a)

  # ------------------------------------------------------------
  # Generate changepoints
  # ------------------------------------------------------------

  changepoint_data <- generate_changepoints(
    n = n,
    n_chg = n_chg,
    delta = delta
  )

  true_chg <- changepoint_data$changepoints

  # ------------------------------------------------------------
  # Generate Gaussian piecewise-constant data
  # ------------------------------------------------------------

  generated_data <- generate_piecewise_normal_2d(
    n = n,
    changepoints = true_chg,
    Sigma = Sigma_true,
    possible_states = ps
  )

  x <- generated_data$x
  true_means <- generated_data$means

  # ------------------------------------------------------------
  # Correct thinning
  # ------------------------------------------------------------

  thinned_correct <- thin_gaussian_2d(
    x = x,
    Sigma = Sigma_true,
    epsilon = epsilon
  )

  x1_correct <- thinned_correct$x1
  x2_correct <- thinned_correct$x2

  # ------------------------------------------------------------
  # Misspecified thinning
  # ------------------------------------------------------------

  thinned_wrong <- thin_gaussian_2d(
    x = x,
    Sigma = Sigma_wrong,
    epsilon = epsilon
  )

  x1_wrong <- thinned_wrong$x1
  x2_wrong <- thinned_wrong$x2

  # ------------------------------------------------------------
  # Estimate changepoints
  # ------------------------------------------------------------

  est_chg_p <- sort(
    jointseg::doDynamicProgramming(
      x,
      K = n_chg
    )$bkp
  )

  est_chg_th_correct <- sort(
    jointseg::doDynamicProgramming(
      x1_correct,
      K = n_chg
    )$bkp
  )

  est_chg_th_wrong <- sort(
    jointseg::doDynamicProgramming(
      x1_wrong,
      K = n_chg
    )$bkp
  )

  # ------------------------------------------------------------
  # Confidence regions
  # ------------------------------------------------------------

  regions_p <- confidence_regions_1D(
    x = x,
    breakpoints = est_chg_p,
    delta = delta,
    alpha = alpha,
    Qtable = Qtable,
    covariance_method = "pooled"
  )

  regions_th_correct <- confidence_regions_thinning(
    x2 = x2_correct,
    breakpoints = est_chg_th_correct,
    Sigma = Sigma_true,
    epsilon = epsilon,
    alpha = alpha,
    delta = delta
  )

  regions_th_wrong <- confidence_regions_thinning(
    x2 = x2_wrong,
    breakpoints = est_chg_th_wrong,
    Sigma = Sigma_wrong,
    epsilon = epsilon,
    alpha = alpha,
    delta = delta
  )

  # ------------------------------------------------------------
  # True signal
  # ------------------------------------------------------------

  true_boundaries <- c(
    0,
    true_chg,
    n
  )

  n_true_segments <- length(true_chg) + 1

  true_signal <- matrix(
    NA_real_,
    nrow = n,
    ncol = 2
  )

  for (i in seq_len(n_true_segments)) {
    current_rows <- seq.int(
      from = true_boundaries[i] + 1L,
      to = true_boundaries[i + 1L]
    )

    true_signal[current_rows, ] <- matrix(
      true_means[i, ],
      nrow = length(current_rows),
      ncol = 2,
      byrow = TRUE
    )
  }

  colnames(true_signal) <- c("mu1", "mu2")

  # ------------------------------------------------------------
  # True targets for each estimated segmentation
  # ------------------------------------------------------------

  targets_p <- get_true_targets(
    true_signal = true_signal,
    breakpoints = est_chg_p
  )

  targets_th_correct <- get_true_targets(
    true_signal = true_signal,
    breakpoints = est_chg_th_correct
  )

  targets_th_wrong <- get_true_targets(
    true_signal = true_signal,
    breakpoints = est_chg_th_wrong
  )

  # ------------------------------------------------------------
  # Segmentwise coverage
  # ------------------------------------------------------------

  coverage_p <- vapply(
    seq_len(nrow(targets_p)),
    function(i) {
      is_in_confidence_region(
        point = targets_p[i, ],
        region = regions_p,
        segment = i
      )
    },
    logical(1)
  )

  coverage_th_correct <- vapply(
    seq_len(nrow(targets_th_correct)),
    function(i) {
      is_in_confidence_region(
        point = targets_th_correct[i, ],
        region = regions_th_correct,
        segment = i
      )
    },
    logical(1)
  )

  coverage_th_wrong <- vapply(
    seq_len(nrow(targets_th_wrong)),
    function(i) {
      is_in_confidence_region(
        point = targets_th_wrong[i, ],
        region = regions_th_wrong,
        segment = i
      )
    },
    logical(1)
  )

  # ------------------------------------------------------------
  # Simultaneous coverage
  # ------------------------------------------------------------

  c(
    PoSIRd = all(
      coverage_p[regions_p$valid]
    ),
    Thinning_correct = all(
      coverage_th_correct[
        regions_th_correct$valid
      ]
    ),
    Thinning_misspecified = all(
      coverage_th_wrong[
        regions_th_wrong$valid
      ]
    )
  )
}

run_simulation_study <- function(
    B,
    n = 1000,
    n_chg = 10,
    delta = 0.05,
    alpha = 0.05,
    epsilon = 0.5,
    a = 2,
    rho = 0.7,
    Sigma_wrong = diag(2),
    Qtable,
    seed = 1535
) {
  set.seed(seed)

  results <- replicate(
    B,
    run_one_simulation(
      n = n,
      n_chg = n_chg,
      delta = delta,
      alpha = alpha,
      epsilon = epsilon,
      a = a,
      rho = rho,
      Sigma_wrong = Sigma_wrong,
      Qtable = Qtable
    )
  )

  # rows = methods, columns = repetitions
  average_simultaneous_coverage <- rowMeans(results)

  list(
    simultaneous_coverage = results,
    average_simultaneous_coverage =
      average_simultaneous_coverage
  )
}

res <- run_simulation_study(
  B = 1000,
  Qtable = Qtable
)

res$average_simultaneous_coverage
