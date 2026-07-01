

### Simulate n trajectories of a POSIR process, version with Rcpp ###

#' Random simulation of trajectories of the POSIR process
#'
#' Simulate n trajectories (as \eqn{\delta} decreases)
#' of the (discretized) 1D or 2D POSIR process.
#'
#' This version calls to C++ written code, with [Rcpp][Rcpp::Rcpp-package] and
#' [RcppArmadillo][RcppArmadillo::RcppArmadillo-package].
#' See [n_traj_simu_1D] and [n_traj_simu_2D] for the R versions.
#'
#' @param n number of generated independent POSIR trajectories.
#' @param deltagrid grid for \eqn{\delta}.
#' @param maxmatrixsize maximum number of i.i.d. random variables
#' to be simultaneously simulated.
#' @param rdistrib distribution (Function) of the i.i.d. random variables
#' that are summed in the POSIR process.
#' @param is_standard whether the POSIR process is to be standardized
#' (in particular if the standard deviation of the distribution is not 1).
#' @inheritParams simulationDir
#' @inheritParams check_grid
#' @seealso [run_simu()], [simulationDir()], [check_grid()], [batch_filename()]
#' @family POSIR process generators
#'
#' @return A matrix with n lines and length(deltagrid) columns.
#' @importFrom foreach foreach
#' @export
#'
#' @examples
#' n_traj_simu(50, 100, seq(10, 1, -1) / 10, 10^6, rnorm, TRUE, 1)
#' n_traj_simu(50, 100, seq(10, 1, -1) / 10, 10^6, rCenteredPareto, FALSE, 2)
n_traj_simu <- function(n, Ndis, deltagrid, maxmatrixsize,
                        rdistrib, is_standard, d, ErLev = .001) {
  if ((d == 1 && n * Ndis > maxmatrixsize) || (d == 2 && Ndis**2 > maxmatrixsize)) {
    log_n_stop("over memory limit", "in n_traj_simu().")
  }
  #gc() # pour libérer la mémoire...
  Intdgrid <- check_grid(Ndis, deltagrid, ErLev)
  if(d==1) {
    Y <- n_traj_simu_1D_C(matrix(rdistrib(n*Ndis),n,Ndis),
                          Intdgrid, is_standard)
    colnames(Y) <- deltagrid
    return(Y)
  }
  if(d==2) {
    Y <- foreach(i=1:n, .combine=rbind) %do% {
      n_traj_simu_2D_C(matrix(rdistrib(Ndis*Ndis),Ndis,Ndis),
                       Intdgrid, is_standard)
    }
    colnames(Y) <- deltagrid
    rownames(Y) <- 1:n
    return(Y)
  }
  # else
  log_n_stop(
    paste("unsupported dimension d=", toString(d), sep = ""),
    "in n_traj_simu()."
  )
  colnames(Y) <- deltagrid
  return(Y)
}

### Simulate n trajectories of a POSIR process, pure R old versions ###

#' Random simulation of trajectories of the 1D POSIR process
#'
#' Simulate n trajectories (as \eqn{\delta} decreases)
#' of the (discretized) 1D POSIR process.
#'
#' This is the R version, see [n_traj_simu] for the C++ version.
#'
#' @inherit n_traj_simu params return
#' @inheritParams simulationDir
#' @inheritParams check_grid
#' @inheritParams n_traj_simu
#' @seealso [run_simu()]
#' @family POSIR process generators
#'
#' @importFrom stats sd
#' @export
#'
#' @examples
#' n_traj_simu_1D(50, 100, seq(10, 1, -1) / 10, 10^6, rnorm, TRUE, .001)
#'
############################### DANIELA NO
n_traj_simu_1D <- function(n, Ndis, deltagrid, maxmatrixsize, rdistrib,
                           is_standard, ErLev) {
  if (n * Ndis > maxmatrixsize) {
    log_n_stop("it's too much !", "in n_traj_simu_1D().")
  }
  grillen <- check_grid(Ndis, deltagrid, ErLev)
  l <- length(deltagrid)
  Y <- matrix(0., n, l)
  colnames(Y) <- deltagrid
  X <- matrix(rdistrib(n * Ndis), n, Ndis)
  # standardisation
  if (!is_standard) {
    sigmas <- apply(X, 1, sd)
    X <- X / sigmas
  }
  X <- as.matrix(cbind(rep(0, n), t(apply(X, 1, cumsum))))
  # Calculs des sups d'accroissements renormalisés Y, en faisant décroître l'écart
  Yprec <- rep(0, n)
  posgrille <- 1
  k <- Ndis
  last_pos <- Ndis + 1
  while (posgrille <= l) {
    dW <- abs(X[, (k + 1):last_pos] - X[, 1:(last_pos - k)]) / sqrt(k)
    Ycur <- apply(cbind(dW, Yprec), 1, max)
    if (k == grillen[posgrille]) {
      Y[, posgrille] <- Ycur
      posgrille <- posgrille + 1
    }
    Yprec <- Ycur
    k <- k - 1
  }
  return(Y)
}

### Simulate n trajectories of a POSIRd process, pure R old versions ###
### Multidimensional extension of POSIR paper/functions

#' Random simulation of trajectories of the dD POSIRd process
#'
#' Simulate n trajectories (as \eqn{\delta} decreases)
#' of the (discretized) POSIRd process.
#'
#' This is the R version, see [n_traj_simu] for the C++ version.
#'
#' @inherit n_traj_simu params return
#' @inheritParams simulationDir
#' @inheritParams check_grid
#' @inheritParams n_traj_simu
#' @param dim dimension of the signal
#' @seealso [run_simu()]
#' @family POSIR process generators
#'
#' @importFrom stats sd
#' @export
#'
#' @examples
#' n_traj_simu_dD(50, 100, 2, seq(10, 1, -1) / 10, 10^6, rnorm, TRUE, .001)
#'
######################## DANIELA YES
n_traj_simu_dD <- function(n, Ndis, dim, deltagrid, maxmatrixsize, rdistrib,
                           is_standard, ErLev) {
  if (n * Ndis * dim > maxmatrixsize) {
    log_n_stop("it's too much !", "in n_traj_simu_dD().")
  }

  grillen <- check_grid(Ndis, deltagrid, ErLev)
  l <- length(deltagrid)

  Y <- matrix(0., n, l)
  colnames(Y) <- deltagrid

  # X[i, t, h] = epsilon_t^(h) for trajectory i, time t, component h
  X <- array(rdistrib(n * Ndis * dim), dim = c(n, Ndis, dim))

  # standardisation
  if (!is_standard) {
    stop("Non-standard multivariate standardisation not implemented yet")
  }

  # Cumulative sums along the time dimension
  S <- array(0, dim = c(n, Ndis + 1, dim))
  for (h in 1:dim) {
    S[, 2:(Ndis + 1), h] <- t(apply(X[, , h], 1, cumsum))
  }

  # Calculs des sups d'accroissements renormalisés Y, en faisant décroître l'écart
  Yprec <- rep(0, n)
  posgrille <- 1
  k <- Ndis
  last_pos <- Ndis + 1

  while (posgrille <= l) {
    # increments[i, a, h] = S_{a+k}^{(h)} - S_a^{(h)}
    increments <- S[, (k + 1):last_pos, , drop = FALSE] -
      S[, 1:(last_pos - k), , drop = FALSE]

    # Euclidean norm over the signal dimension h
    dW <- sqrt(apply(increments^2, c(1, 2), sum)) / sqrt(k)

    # Supremum over all intervals of length k, and previous larger lengths
    Ycur <- apply(cbind(dW, Yprec), 1, max)

    if (k == grillen[posgrille]) {
      Y[, posgrille] <- Ycur
      posgrille <- posgrille + 1
    }
    Yprec <- Ycur
    k <- k - 1
  }
  return(Y)
}



#' Random simulation of trajectories of the 2D POSIR process
#'
#' Simulate n trajectories (as \eqn{\delta} decreases)
#' of the (discretized) 2D POSIR process.
#'
#' This is the R version, see [n_traj_simu] for the C++ version.
#'
#' @inherit n_traj_simu params return
#' @inheritParams simulationDir
#' @inheritParams check_grid
#' @inheritParams n_traj_simu
#' @seealso [run_simu()]
#' @family POSIR process generators
#'
#' @importFrom stats sd
#' @export
#'
#' @examples
#' n_traj_simu_2D(50, 100, seq(10, 1, -1) / 10, 10^6, rnorm, TRUE, .001)
n_traj_simu_2D <- function(n, Ndis, deltagrid, maxmatrixsize, rdistrib,
                           is_standard, ErLev) {
  # TODO éventuel : version par blocs ? probablement inutile
  if (Ndis^2 > maxmatrixsize) {
    log_n_stop("not currently possible", "in n_traj_simu_2D().")
  }
  pos_grid <- check_grid(Ndis, deltagrid, ErLev)
  l <- length(deltagrid)
  res <- matrix(rep(0, n * l), n, l)
  colnames(res) <- deltagrid
  for (a in 1:n) {
    X <- matrix(rdistrib(Ndis^2), Ndis, Ndis)
    # standardisation
    if (!is_standard) X <- X / sd(as.vector(X))
    # sommes cumulatives
    X <- t(apply(X, 1, cumsum))
    X <- apply(X, 2, cumsum)
    # ajouter des 0 au début de ces sommes cumulatives
    X <- as.matrix(rbind(rep(0, Ndis + 1), cbind(rep(0, Ndis), X)))
    prev_pos <- Ndis
    last_pos <- Ndis + 1
    cur_sup <- abs(X[Ndis + 1, Ndis + 1]) / Ndis # cas delta=1
    # ajouter cette dernière valeur en début de résultat (pour delta=1) ? FIXME
    for (i in 1:l) {
      next_pos <- pos_grid[i]
      for (j in seq(prev_pos, next_pos, -1)) {
        # log_debug(paste("  j =", toString(j), "\n"))
        # log_debug("in n_traj_simu_2D()")
        for(k in seq(Ndis, next_pos, -1)) {
          dX <- abs(X[(j + 1):last_pos, (k + 1):last_pos]
                    + X[1:(last_pos - j), 1:(last_pos - k)]
                    - X[(j + 1):last_pos, 1:(last_pos - k)]
                    - X[1:(last_pos - j), (k + 1):last_pos]
                    )/sqrt(j*k)
          cur_sup <- max(cur_sup, dX)
        }
      }
      for (j in seq(Ndis, prev_pos, -1)) {
        # log_debug(paste("  j =", toString(j), "\n"))
        # log_debug("in n_traj_simu_2D()")
        for(k in seq(prev_pos, next_pos, -1)) {
          dX <- abs(X[(j + 1):last_pos, (k + 1):last_pos]
                    + X[1:(last_pos - j), 1:(last_pos - k)]
                    - X[(j + 1):last_pos, 1:(last_pos - k)]
                    - X[1:(last_pos - j), (k + 1):last_pos]
                    )/sqrt(j*k)
          cur_sup <- max(cur_sup, dX)
        }
      }
      res[a, i] <- cur_sup
      prev_pos <- next_pos
    }
  }
  return(res)
}



