# PoSIRd semi-synthetic CNV analysis
# Adapted from recursive_realdata.Rmd and recursive_sims.Rmd
# Run from the root of the posir project.

library(jointseg)
library(posir)
library(acnr)
library(here)

source(here::here('inst', 'scripts', 'posird', 'utils_competitors.R'))
Qtable <- posir:::read_quantiles(here::here('Simulations', 'Table_quantiles.txt'))

# ---- 1. Resample real CNV data (reproducible semi-synthetic example) ----
affyDat <- acnr::loadCnRegionData(dataSet = 'GSE11976', tumorFraction = 1)
dat <- subset(affyDat, region %in% c('(0,1)', '(1,1)', '(1,2)'))
set.seed(20250603)
sim <- getCopyNumberDataByResampling(10000, 5, regData = dat, minLength = 200)
is_het <- !is.na(sim$profile$genotype) & sim$profile$genotype == 0.5
Y <- sim$profile$c
AR <- sim$profile$b
DH <- 2 * abs(AR - 0.5)
valid_rows <- is_het & is.finite(Y) & is.finite(DH)
# Original Rmd retains heterozygous observations; this also removes any nonfinite measurements.
if (any(is_het & !valid_rows)) warning('Some heterozygous observations had nonfinite measurements and were removed.')
x <- cbind(CN = Y[valid_rows], DH = DH[valid_rows])
n <- nrow(x)
true_bkp <- as.integer(cumsum(valid_rows)[sim$bkp])
true_bkp <- sort(unique(true_bkp[true_bkp > 0 & true_bkp < n]))

# ---- 2. Parameters ----
delta <- 0.005
alpha <- 0.05
K_high <- 10L
ndelta <- ceiling(n * delta)
stopifnot(n >= 2 * ndelta)

# Same deterministic left-to-right filter as recursive_sims.Rmd.
# It removes candidate boundaries; it does NOT remove data points.
filter_changepoints <- function(bkp, n, ndelta) {
  bkp <- sort(unique(as.integer(bkp)))
  bkp <- bkp[bkp > 0 & bkp < n]
  retained <- integer(0)
  previous <- 0L
  for (b in bkp) {
    if (b - previous >= ndelta && n - b >= ndelta) {
      retained <- c(retained, b)
      previous <- b
    }
  }
  retained
}

check_lengths <- function(bkp, n, ndelta) {
  lengths <- diff(c(0L, bkp, n))
  if (any(lengths < ndelta)) stop('Some segments violate the minimum length constraint.')
  lengths
}

# ---- 3. Initial segmentation and length filtering ----
bkp_initial <- sort(unique(as.integer(doDynamicProgramming(x, K = K_high)$bkp)))
bkp_filtered <- filter_changepoints(bkp_initial, n, ndelta)
lengths_filtered <- check_lengths(bkp_filtered, n, ndelta)

get_regions <- function(bkp) {
  confidence_regions_1D(
    x = x, breakpoints = bkp, delta = delta,
    alpha = alpha, Qtable = Qtable,
    covariance_method = 'diff_mom_entrywise',
    covariance_blocks = 50
  )
}

posird_filtered <- get_regions(bkp_filtered)
stopifnot(all(posird_filtered$valid))

# ---- 4. Confidence region overlap (shared covariance) ----
confidence_regions_overlap <- function(region, i, j) {
  if (!isTRUE(region$valid[i]) || !isTRUE(region$valid[j])) {
    stop('Cannot compare invalid confidence regions.')
  }
  d2 <- as.numeric(mahalanobis(
    matrix(region$means[i, ], nrow = 1),
    center = region$means[j, ], cov = region$covariance
  ))
  radius_sum <- region$quantile * (
    1 / sqrt(region$lengths[i]) + 1 / sqrt(region$lengths[j])
  )
  sqrt(max(0, d2)) <= radius_sum
}

# One-step pruning, consistent with the supplied Rmd: compare every
# adjacent pair using the SAME pre-pruning regions, then remove boundaries.
prune_breakpoints_from_posir <- function(region) {
  bkp <- as.integer(region$breakpoints)
  if (!length(bkp)) return(integer(0))
  overlaps <- vapply(seq_along(bkp), function(i) {
    confidence_regions_overlap(region, i, i + 1L)
  }, logical(1))
  bkp[!overlaps]
}

bkp_pruned <- prune_breakpoints_from_posir(posird_filtered)
lengths_pruned <- check_lengths(bkp_pruned, n, ndelta)
posird_pruned <- get_regions(bkp_pruned)
stopifnot(all(posird_pruned$valid))

# ---- 5. Results ----
summary_table <- data.frame(
  segmentation = c('True', 'Initial', 'Filtered', 'Pruned'),
  n_changepoints = c(length(true_bkp), length(bkp_initial),
                     length(bkp_filtered), length(bkp_pruned)),
  min_segment_length = c(min(diff(c(0, true_bkp, n))),
                         min(diff(c(0, bkp_initial, n))),
                         min(lengths_filtered), min(lengths_pruned))
)
print(summary_table, row.names = FALSE)
cat('n =', n, '; delta =', delta, '; minimum length =', ndelta, '\n')
cat('True:    ', true_bkp, '\n')
cat('Initial: ', bkp_initial, '\n')
cat('Filtered:', bkp_filtered, '\n')
cat('Pruned:  ', bkp_pruned, '\n')
cat('Removed for length:', setdiff(bkp_initial, bkp_filtered), '\n')
cat('Removed for overlap:', setdiff(bkp_filtered, bkp_pruned), '\n')

# ---- 6. Plots suitable for the presentation ----
out_dir <- 'realdata_posird_figures'
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pdf(file.path(out_dir, 'cnv_segmentations.pdf'), width = 12, height = 6)
par(mfrow = c(2, 1), mar = c(3, 4, 1.5, 1))
for (j in 1:2) {
  plot(x[, j], pch = 16, cex = 0.3, xlab = 'Index',
       ylab = colnames(x)[j], col = 'grey40')
  abline(v = true_bkp, col = 'firebrick', lwd = 1.4)
  abline(v = bkp_filtered, col = 'steelblue', lty = 2)
  abline(v = bkp_pruned, col = 'darkgreen', lty = 3, lwd = 1.2)
  if (j == 1) legend('topright',
                     legend = c('True', 'Estimated', 'Pruned'),
                     col = c('firebrick', 'steelblue', 'darkgreen'),
                     lty = c(1, 2, 3), bty = 'n', cex = 0.8)
}
dev.off()

mu1_grid <- seq(min(x[, 1]), max(x[, 1]), length.out = 400)
mu2_grid <- seq(min(x[, 2]), max(x[, 2]), length.out = 400)
for (nm in c('filtered', 'pruned')) {
  obj <- if (nm == 'filtered') posird_filtered else posird_pruned
  pdf(file.path(out_dir, paste0('confidence_regions_', nm, '.pdf')),
      width = 8, height = 7)
  plot(obj, mu1_grid = mu1_grid, mu2_grid = mu2_grid,
       main = paste('PoSIRd confidence regions:', nm),
       show_labels = TRUE, data_cex = 0.5, show_centers = FALSE)
  dev.off()
}

write.csv(summary_table, file.path(out_dir, 'summary.csv'), row.names = FALSE)
saveRDS(list(x = x, true_bkp = true_bkp, initial = bkp_initial,
             filtered = bkp_filtered, pruned = bkp_pruned,
             posird_filtered = posird_filtered,
             posird_pruned = posird_pruned, summary = summary_table),
        file.path(out_dir, 'analysis_results.rds'))

# Common axis limits for both plots
mu1_grid <- seq(min(x[, 1]), max(x[, 1]), length.out = 400)
mu2_grid <- seq(min(x[, 2]), max(x[, 2]), length.out = 400)

# Before pruning: 8 changepoints, 9 segments
pdf("ellipses_filtered.pdf", width = 7, height = 6)

plot(
  posird_filtered,
  mu1_grid = mu1_grid,
  mu2_grid = mu2_grid,
  main = "",
  show_labels = TRUE,
  data_cex = 0.4,
  show_centers = FALSE
)

dev.off()

# After pruning: 5 changepoints, 6 segments
pdf("ellipses_pruned.pdf", width = 7, height = 6)

plot(
  posird_pruned,
  mu1_grid = mu1_grid,
  mu2_grid = mu2_grid,
  main = "",
  show_labels = TRUE,
  data_cex = 0.4,
  show_centers = FALSE
)

dev.off()

xlim <- range(x[, 1]) + c(-0.3, 0.3)
ylim <- range(x[, 2]) + c(-0.3, 0.3)

plot_posir_story(
    posird_filtered,
    stage = "regions",
    xlim = NULL,
    ylim = c(0,1),
    cols = NULL,
    data_cex = 0.8,
    label_cex = 1.2,
    label_angles = NULL,
    lwd = 3,
    asp = 1,
    show_legend = T
)

plot_posir_story(
  posird_pruned,
  stage = "regions",
  xlim = NULL,
  ylim = c(0,1),
  cols = NULL,
  data_cex = 0.8,
  label_cex = 1.2,
  label_angles = NULL,
  lwd = 3,
  asp = 1,
  show_legend = T
)

pdf(
  file.path(out_dir, "cnv_segmentations.pdf"),
  width = 6,
  height = 6,
  pointsize = 11
)

par(
  mfrow = c(2, 1),
  mar = c(3, 3.5, 1, 0.5),
  mgp = c(1.8, 0.5, 0),
  tcl = -0.25
)

for (j in 1:2) {

  plot(
    x[, j],
    pch = 16,
    cex = 0.35,
    xlab = "Index",
    ylab = c("CN", "DH")[j],
    col = adjustcolor("grey40", alpha.f = 0.35),
    cex.axis = 0.9,
    cex.lab = 1.1
  )

  # True breakpoints
  abline(
    v = true_bkp,
    col = "#D55E00",
    lwd = 2.5,
    lty = 1
  )

  # Estimated breakpoints
  abline(
    v = bkp_filtered,
    col = "#0072B2",
    lwd = 2.2,
    lty = 2
  )

  # Pruned breakpoints
  abline(
    v = bkp_pruned,
    col = "#009E73",
    lwd = 2.5,
    lty = 3
  )

  if (j == 1) {
    legend(
      "topright",
      legend = c("True", "Estimated", "Pruned"),
      col = c("#D55E00", "#0072B2", "#009E73"),
      lty = c(1, 2, 3),
      lwd = c(2.5, 2.2, 2.5),
      bty = "n",
      cex = 0.85
    )
  }
}

dev.off()
