
### Preparation of the POSIR simulation ###

#' Define the directory for a POSIR simulation
#'
#' Compute the relative path (directory name) for stocking the batch files
#' in a POSIR simulation.
#'
#' If the directory does not exist, it is created.
#'
#' @param sim_root_dir root directory of the simulations. Both "" and "."
#' stand for the current working directory.
#' @param Ndis discretisation parameter of the POSIR process.
#' @param gridname name of the grid for \eqn{\delta}.
#' @param NameDis name, or short description (of type string),
#' of the distribution of the i.i.d. random variables
#' that are summed in the POSIR process.
#' @param d domain dimension parameter of the POSIR process: 1 for temporal
#' intervals, 2 for spatial rectangles.
#' @param signal_dim dimension of the multivariate signal. Defaults to 1.
#'
#' @return Something similar to
#' "1D_Discretization_1080_Gaussian_grid_de .875 à .25.txt".
#' @seealso [batch_filename()]
#'
#' @examples
#' \dontrun{
#' sim_dir <- local_init_testing(dirNameBase = "posir_test",
#'                               envir=sys.frame(sys.nframe()))
#' simulationDir(sim_dir, 1080, "de .875 à .25", "Gaussian", 1)
#' # some code
#' withr::deferred_run(envir=sys.frame(sys.nframe()))
#' }
#' @keywords internal
simulationDir <- function(sim_root_dir, Ndis, gridname, NameDis, d,
                          signal_dim = 1) {
  NameD <- paste("Discretization", toString(Ndis), NameDis, "grid",
                 gridname,
                 sep = "_"
  )
  NameD <- paste(toString(d), "D_", NameD, sep = "")
  if (signal_dim > 1) {
    NameD <- paste0(NameD, "_signal_dim_", signal_dim)
  }
  if (sim_root_dir != "") {
    if (!dir.exists(sim_root_dir)) {
      log_debug(paste("simulationDir(): creating directory", sim_root_dir))
      dir.create(sim_root_dir)
    }
    if (substring(sim_root_dir, nchar(sim_root_dir)) != "/") {
      NameD <- paste("/", NameD, sep = "")
    }
    NameD <- paste(sim_root_dir, NameD, sep = "")
  }
  if (!dir.exists(NameD)) {
    log_debug(paste("simulationDir(): creating directory", NameD))
    dir.create(NameD)
  }
  return(NameD)
}

### Batches gestion ###

#' Define the name (string) of a batch file in a POSIR simulation
#'
#' Compute the file name (without the path) of a batch in a POSIR simulation.
#'
#' @param j number of the current batch.
#' @param Ltot number of zeros to be added in the file name of batch j=1.
#' @param NameBaseBatch prefix (beginning) of the file name for batches,
#' without the path. See [simulationDir()] for the directory used.
#'
#' @return Something similar to "Batch50_3_0025.txt" or "Batch50_3_2546.txt".
#' @seealso [simulationDir()], [write_batch()], [run_simu()]
#'
#' @examples
#' \dontrun{
#' batch_filename(25, 3, "Batch50_3_")
#' }
#' \dontrun{
#' Batchsize <- 50
#' Nbatch <- 4000
#' Ltot <- floor(log10(Nbatch))
#' batch_filename(2546, Ltot, paste("Batch", toString(Batchsize), "_",
#'   toString(Ltot), "_",
#'   sep = ""
#' ))
#' }
#' @keywords internal
batch_filename <- function(j, Ltot, NameBaseBatch) {
  return(do.call("paste", as.list(c(NameBaseBatch,
                                    rep("0", Ltot - floor(log10(j))),
                                    toString(j), ".txt",
                                    sep = ""
  ))))
}

#' Launch a set a batchs in parallel
#'
#' Simulate and save batches of trajectories of the 1D or 2D POSIR process.
#' If possible, the batches are run in parallel.
#'
#' @param beginpos beginning number of the range of batches to be simulated.
#' @param endpos end number of the range of batches to be simulated.
#' @param Nbatch total number of batches planned in the simulation.
#' @param Batchsize number of trajectories simulated in each batch.
#' @param BaseNameF Base of the file name (with path) for the batches.
#' @inheritParams n_traj_simu
#' @inheritParams write_batch
#' @inheritParams simulationDir
#' @inheritParams check_grid
#' @inheritParams batch_filename
#'
#' @return NULL
#' @importFrom foreach foreach
#' @importFrom doFuture %dofuture%
#' @export
#'
#' @examples
#' logger::log_info("Running examples for run_batches_range()",
#'                  namespace = "posir")
#' sim_dir <- local_init_testing(dirNameBase = "run_batches_example",
#'                               do_multi = TRUE,
#'                               envir=sys.frame(sys.nframe()))
#' BaseNameF <- paste(sim_dir, "/Batch50_1_", sep = "")
#' run_batches_range(
#'   1, 3, 20, 50, 100, seq(10, 1, -1) / 10, 10^6,
#'   BaseNameF, rnorm, TRUE, 2
#' )
#' withr::deferred_run(envir=sys.frame(sys.nframe()))
run_batches_range <- function(beginpos, endpos, Nbatch, Batchsize, Ndis, deltagrid,
                              maxmatrixsize, BaseNameF, rdistrib,
                              is_standard, d, ErLev = .001, signal_dim = 1) {
  if (beginpos != round(beginpos) || endpos != round(endpos) || Batchsize != round(Batchsize) ||
      Nbatch != round(Nbatch) || Ndis != round(Ndis)) {
    log_n_stop("invalid non-integer parameter", "in run_batches().")
  }
  if (beginpos < 1 || endpos < beginpos || Nbatch < endpos || Batchsize < 1 || Ndis < 2) {
    log_n_stop("invalid range or non-positive parameter", "in run_batches().")
  }
  Ltot <- floor(log10(Nbatch))
  T1 <- Sys.time()
  j <- NULL
  foreach(
    j = beginpos:endpos,
    .inorder = FALSE,
    .options.future = list(
      seed = TRUE,
      packages = "posir"
    )
  ) %dofuture% {
    check_table_file(
      batch_filename(j, Ltot, BaseNameF),
      deltagrid, Batchsize
    )
  }

  foreach(
    j = beginpos:endpos,
    .inorder = FALSE,
    .options.future = list(
      seed = TRUE,
      packages = "posir"
    )
  ) %dofuture% {
    write_batch(
      Batchsize, batch_filename(j, Ltot, BaseNameF), Ndis, deltagrid,
      paste("Simulation of batch", toString(j), "among", toString(Nbatch)),
      maxmatrixsize, rdistrib, is_standard, d, ErLev, signal_dim
    )
  } # par tout-à-fait future_compatible (write in parallelling...) FIXME
  T2 <- Sys.time()
  log_info(paste(
    "Elapsed time at the end of the run:",
    toString(round(difftime(T2, T1, units = "secs"))),
    "s.\n"
  ))
  return(NULL)
}

# Writes tables in Files

#' Compute and save a batch of POSIR trajectories
#'
#' Simulates a batch of 1D or 2D POSIR process trajectories,
#' and exports them in a file, only if the batch file does not already exist.
#'
#' @param NameF batch file name.
#' @param message message to be printed if the files does not already exist.
#' @inheritParams n_traj_simu
#' @inheritParams simulationDir
#' @inheritParams check_grid
#' @seealso [run_simu()], [n_traj_simu()], [simulationDir()], [check_grid()], [batch_filename()]
#'
#' @return NULL
#' @importFrom utils write.table
#'
#' @examples
#' \dontrun{
#' logger::log_info("Running examples for write_batch()",
#'                  namespace = "posir")
#' sim_dir <- local_init_testing(dirNameBase = "posir_test",
#'                               envir=sys.frame(sys.nframe()))
#' write_batch(
#'   50, paste(sim_dir,"Batch50_1_01.txt",sep="/"), 100, seq(10, 1, -1) / 10,
#'   "go", 10^6, rnorm, TRUE, 1, .001
#' )
#' withr::deferred_run(envir=sys.frame(sys.nframe()))
#' }
#' @keywords internal
write_batch <- function(n, NameF, Ndis, deltagrid, message, maxmatrixsize,
                        rdistrib, is_standard, d, ErLev,
                        signal_dim = 1) {
  if (!file.exists(NameF)) {
    log_info(message)
    res = n_traj_simu(n, Ndis, deltagrid, maxmatrixsize, rdistrib,
                      is_standard, d, ErLev,
                      signal_dim)
    write.table(res, NameF)
  }
  return(NULL)
}

### Main run function for POSIR simulation ###

#' Launch the simulation of trajectories of the 1D or 2D POSIR process
#'
#' The trajectories are regrouped by batches, which are parallelised.
#'
#' @param Ntot initial number of trajectories planned to be simulated.
#' @inheritParams run_batches_range
#' @inheritParams n_traj_simu
#' @inheritParams write_batch
#' @inheritParams simulationDir
#' @inheritParams check_grid
#' @inheritParams batch_filename
#' @seealso [run_batches_range()], [write_batch()], [compute_quantiles()], [compute_error_levels()]
#'
#' @return A matrix of simulated POSIR process trajectories (in lines).
#' @importFrom foreach foreach %do%
#' @importFrom stats rnorm
#' @export
#'
#' @examples
#' logger::log_info("Running examples for run_simu()",
#'                  namespace = "posir")
#' sim_dir <- local_init_testing(do_multi = TRUE, envir=sys.frame(sys.nframe()))
#' run_simu(500, 50, 100, seq(10, 1, -1) / 10, "de 1 à .1", 10^6, sim_dir)
#' future::plan(future::sequential)
#' run_simu(30, 10, 100, seq(10, 1, -1) / 10, "de 1 à .1", 10^6, sim_dir, d = 2)
#' withr::deferred_run(envir=sys.frame(sys.nframe()))
run_simu <- function(Ntot, Batchsize, Ndis, deltagrid, gridname, maxmatrixsize,
                     sim_root_dir = "", NameBaseBatch = "Batch", rdistrib = rnorm,
                     NameDis = "Gaussian", is_standard = TRUE, d = 1, ErLev = .001,
                     signal_dim = 1) {
  local_log_init(sim_root_dir)
  if (Ntot != round(Ntot) || Batchsize != round(Batchsize)) {
    log_n_stop("invalid non-integer parameter", "in run_simu().")
  }
  if (Ntot < 1 || Batchsize < 1) {
    log_n_stop("invalid non-positive parameter", "in run_simu().")
  }
  Nbatch <- ceiling(Ntot / Batchsize)
  Ltot <- floor(log10(Nbatch))
  NameF <- paste(simulationDir(sim_root_dir, Ndis, gridname, NameDis, d, signal_dim),
    paste(NameBaseBatch, toString(Batchsize), "_",
      toString(Ltot), "_",
      sep = ""
    ),
    sep = "/"
  )
  run_batches_range(
    1, Nbatch, Nbatch, Batchsize, Ndis, deltagrid, maxmatrixsize,
    NameF, rdistrib, is_standard, d, ErLev, signal_dim
  )
  j <- 0
  Z <- as.matrix(foreach(
    j = 1:Nbatch, .export = "batch_filename",
    .combine = rbind, .inorder = FALSE
  ) %do% {
    read.table(batch_filename(j, Ltot, NameF))
  })
  colnames(Z) <- deltagrid
  rownames(Z) <- 1:(Batchsize * Nbatch) # sinon il y a n'importe quoi
  return(Z)
}
