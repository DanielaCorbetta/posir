
#' Compute the Monte-Carlo estimations of the quantiles of a discretized process
#'
#' @param Z matrix of simulated trajectories of a discretized random process.
#' @param alphagrid grid of values of the error level \eqn{\alpha}.
#' @param precQ how many decimals for the quantile estimations
#' (used in particular in exporting the table to a file).
#' @seealso [run_simu()], [compute_quantiles()], [extract_error_levels()]
#'
#' @return A table of Quantile Estimations.
#' @export
#'
#' @examples
#' logger::log_info("Running examples for extract_quantiles()",
#'                  namespace = "posir")
#' Ndis <- 100
#' Ntraj <- 1000
#' X <- matrix(rnorm(Ndis * Ntraj), Ntraj, Ndis)
#' W <- t(apply(X, 1, cumsum)) / sqrt(Ndis)
#' # Ntraj trajectories of a discretized Brownian motion on [0, 1],
#' # with discretization times (1:Ndis)/Ndis
#' local_init_testing(do_temp = FALSE, do_log = FALSE, do_multi = TRUE,
#'                    envir=sys.frame(sys.nframe()))
#' extract_quantiles(W, seq(10, 1, -1) / 20, 3)
#' withr::deferred_run(envir=sys.frame(sys.nframe()))
extract_quantiles <- function(Z, alphagrid, precQ) {
  Ntraj <- dim(Z)[1]
  # from above:
  # #' @importFrom foreach foreach
  # #' @importFrom doFuture %dofuture%
  # l <- dim(Z)[2]
  # j <- NULL
  # # Tri en vue du calcul des quantiles empiriques
  # Q <- foreach(
  #   j = 1:l, .combine = "cbind", .multicombine = TRUE,
  #   .options.future = list(seed = TRUE)
  # ) %dofuture% {
  #   sort(Z[, j])
  # }
  # from above:
  # #' @importFrom future.apply future_apply
  # #' Q <- future.apply::future_apply(Z,2,sort)
  Q <- apply(Z,2,sort)
  # Position des quantiles empiriques, par interpolation linéaire
  aux <- Ntraj * (1 - alphagrid)
  Posalpha <- floor(aux)
  Q <- Q[Posalpha, ] * (Posalpha + 1 - aux) + Q[(Posalpha + 1), ] * (aux - Posalpha)
  Q <- round(Q, precQ)
  rownames(Q) <- alphagrid
  colnames(Q) <- colnames(Z)
  return(Q)
}

#' Export quantiles tables to files
#'
#' @param Q quantile table.
#' @param NameF file name (to export to).
#' @param Ntraj number of trajectories
#' used in the Monte-Carlo estimation of the quantiles.
#' @inheritParams simulationDir
#' @seealso [extract_quantiles()], [compute_quantiles()], [simulationDir()]
#'
#' @return NULL
#' @importFrom utils write.table
#'
#' @keywords internal
write_quantiles <- function(Q, NameF, Ndis, Ntraj, d) {
  write("######\n#", NameF)
  write(
    paste("# Quantile Table for the ", toString(d), "D POSIR process, \n",
          "# with alpha in rownames et delta in colnames",
          sep = ""
    ),
    NameF,
    append = TRUE
  )
  write("#\n# Simulation parameters\n#", NameF, append = TRUE)
  write(paste("# Ndiscretisation =", toString(Ndis)), NameF, append = TRUE)
  write(paste("# Ntraj =", toString(Ntraj)), NameF, append = TRUE)
  write("#\n######", NameF, append = TRUE)
  write.table(Q, NameF, append = TRUE)
  return(NULL)
}

### Main function for quantiles ###

#' Computes and save a quantiles table of the POSIR process
#'
#' Computes and save Monte-Carlo estimations
#' of the quantiles of the 1D or 2D POSIR process,
#' only if the quantile file does not already exist.
#'
#' @param filenamebasis file name for the quantile table (without the path).
#' @inheritParams run_simu
#' @inheritParams write_batch
#' @inheritParams simulationDir
#' @inheritParams check_grid
#' @inheritParams batch_filename
#' @inheritParams extract_quantiles
#' @seealso [run_simu()], [extract_quantiles()], [write_quantiles()], [compute_error_levels()]
#'
#' @return A quantile table for the 1D or 2D POSIR process,
#' only if the quantile file does not already exist,
#' else the function returns NULL.
#' @importFrom stats rnorm
#' @export
#'
#' @examples
#' logger::log_info("Running examples for compute_quantiles()",
#'                  namespace = "posir")
#' sim_dir <- local_init_testing(do_multi = TRUE,
#'                               envir=sys.frame(sys.nframe()))
#' compute_quantiles(
#'   5000, 500, 100, seq(10, 1, -1) / 10,
#'   "de 1 à .1", seq(10, 1, -1) / 20, 10^6, sim_dir
#' )
#' future::plan(future::sequential)
#' compute_quantiles(100, 50, 100, seq(10, 1, -1) / 10,
#'   "de 1 à .1", seq(10, 1, -1) / 20, 10^6, sim_dir,
#'   d = 2
#' )
#' withr::deferred_run(envir=sys.frame(sys.nframe()))
compute_quantiles <- function(Ntot, Batchsize, Ndis, deltagrid, gridname,
                              alphagrid, maxmatrixsize,
                              sim_root_dir = "",
                              filenamebasis = "Table_quantiles.txt",
                              NameBaseBatch = "Batch",
                              precQ = floor(log10(Ntot) / 2) + 2,
                              d = 1, ErLev = .001, signal_dim = 1) {
  local_log_init(sim_root_dir)
  NameF <- paste(simulationDir(sim_root_dir, Ndis, gridname,
                               NameDis = "Gaussian", d, signal_dim),
    filenamebasis,
    sep = "/"
  )
  if (!check_table_file(NameF, deltagrid, length(alphagrid))) {
    Z <- run_simu(
      Ntot, Batchsize, Ndis, deltagrid,
      gridname, maxmatrixsize, sim_root_dir, NameBaseBatch,
      rnorm, "Gaussian", TRUE,
      d, ErLev, signal_dim
    )
    log_info("Sorting now the trajectories")
    log_debug("in compute_quantiles().")
    Ntraj <- dim(Z)[1]
    Q <- extract_quantiles(Z, alphagrid, precQ)
    write_quantiles(Q, NameF, Ndis, Ntraj, d)
    log_debug("compute_quantiles(): done.")
    return(Q)
  }
  return(NULL)
}
