library(mvtnorm)
library(jointseg)

devtools::load_all()

Qtable_path <- here::here(
  "Simulations",
  "Table_quantiles.txt"
)

Qtable <- read_quantiles(Qtable_path)

# Set parameters
n <- 500
n_chg <- 5
delta <- 0.05
alpha <- 0.05
rho <- 0.5 # correlation among signal dimensions
s <- 1 # variance of signal components
Sigma <- matrix(c(s, rho, rho, s), nrow=2)
possible_states <- c(0, 1, 2, -1, -2)

# Randomly choose changepoints
set.seed(1112)
n_segments <- n_chg + 1
min_length <- ceiling(delta * n)

# Give every segment the minimum required length, then randomly distribute
# the remaining observations among the segments.
remaining_length <- n - n_segments * min_length

segment_lengths <- min_length +
  as.vector(
    rmultinom(
      n = 1,
      size = remaining_length,
      prob = rep(1, n_segments)
    )
  )

# Each changepoint is the last observation of a segment.
true_chg <- cumsum(segment_lengths)[-n_segments]

true_chg
segment_lengths

# Generate data
true_means <- matrix(
  NA_real_,
  nrow = n_segments,
  ncol = 2
)

true_means[1, ] <- sample(
  possible_states,
  size = 2,
  replace = TRUE
)

for (i in 2:n_segments) {
  true_means[i, 1] <- sample(
    setdiff(possible_states, true_means[i - 1, 1]),
    size = 1
  )

  true_means[i, 2] <- sample(
    setdiff(possible_states, true_means[i - 1, 2]),
    size = 1
  )
}

colnames(true_means) <- c("mu1", "mu2")

true_means

set.seed(1418)
boundaries <- c(0, true_chg, n)

x <- matrix(
  NA_real_,
  nrow = n,
  ncol = 2
)

for (i in seq_len(n_segments)) {
  current_rows <- (boundaries[i] + 1):boundaries[i + 1]

  x[current_rows, ] <- rmvnorm(
    n = length(current_rows),
    mean = true_means[i, ],
    sigma = Sigma
  )
}

colnames(x) <- c("Y1", "Y2")

# See posir result with true segmentation
regions <- confidence_regions_1D(
  x = x,
  breakpoints = true_chg,
  delta = delta,
  alpha = alpha,
  Qtable = Qtable
)

mu1_grid <- seq(-2.8, 3, length.out = 300)
mu2_grid <- seq(-2.8, 3, length.out = 300)

plot(
  regions,
  mu1_grid = mu1_grid,
  mu2_grid = mu2_grid,
  main = "Simultaneous PoSIRd confidence regions",
  show_labels = TRUE
)

points(
  true_means[, 1],
  true_means[, 2],
  pch = 4,
  col = "black",
  lwd = 3,
  cex = 1.2
)

plot_posir1DCR_3d(
  regions,
  main = "Two-dimensional signal over time"
)

plot_posir1DCR_3d(
  regions,
  time_aspect = 1.8,
  main = "Simultaneous PoSIRd confidence regions",
  region_location = "middle",
  show_means = TRUE,
  show_jumps = TRUE,
  show_regions = TRUE,
  show_data = TRUE
)

par(mfrow=c(1,1))
plot_posir1DCR_component(
  regions,
  component = 1,
  box = FALSE
)

plot_posir1DCR_component(
  regions,
  component = 2,
  box = FALSE
)
