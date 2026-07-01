#ifndef C_N_TRAJ_SIMU_H
#define C_N_TRAJ_SIMU_H

#include <string>
#include <sstream>
#include <RcppArmadillo.h>
using namespace Rcpp;

/// Log and raise an exception (stop)
/**
 * The message is logged at level ERROR, then stop(message) is called.
 * An optional additional message (details) can be given,
 * to be logged at level DEBUG.
 *
 * Calls to R functions from package logger.
 *
 * See R version [posir::log_n_stop()] for details.
 *
 * @param message text (string) to be logged and passed to stop().
 * @param m_details optional additional details to be log at DEBUG level.
 *
 */
void log_n_stop_C(std::string message, std::string m_details="");

/// Pointwise multiplication of a matrix by a colvec
/**
 * Multiply the rows of a matrix by the coordinates of a column vector.
 *
 * The number of rows in the matrix has to be the same as the length of the vector.
 *
 * @param M a matrix.
 * @param v a vector.
 * @return A matrix where each row of M has been multiplied
 * by the corresponding coordinate of v.
 */
arma::mat RowMult(arma::mat const & M, arma::vec const & v);

/// Pointwise division of a matrix by a colvec
/**
 * Divide the rows of a matrix by the coordinates of a column vector.
 *
 * The number of rows in the matrix has to be the same as the length of the vector.
 *
 * @param M a matrix.
 * @param v a vector.
 * @return A matrix where each row of M has been divided
 * by the corresponding coordinate of v.
 */
arma::mat RowDiv(arma::mat const & M, arma::vec const & v);

/// Compute a batch of trajectories of the 1D POSIR process
/**
 * The function implements the functional G
 * from paper "Inference post region selection"
 * by Bontemps, Bachoc, and Neuvial, 2024;
 * it takes a matrix whose lines are from simulated white noise trajectories
 * and computes from it the trajectories of the associated 1D POSIR process.
 *
 * @param Z a matrix whose lines are simulated white noise trajectories.
 * @param Intdgrid decreasing vector of the integers \eqn{1\le k\le Ndis}
 * such that \eqn{k/Ndis} is in the grid for \eqn{\delta}.
 * Generally obtained as Ndis*deltagrid,
 * where the vector dgrid is the desired grid for \eqn{\delta}.
 * (Beware: no checking is done in this internal function.)
 * @param is_standard whether the POSIR process is to be standardized
 * (in particular if the standard deviation of rdistrib is not 1).
 */
arma::mat n_traj_simu_1D_C(arma::mat const & Z, arma::ivec const & Intdgrid,
                           bool is_standard);

/// Compute one 2D trajectory of the 2D POSIR process
/**
 * The function implements the 2D version of the functional G
 * from paper "Inference post region selection"
 * by Bontemps, Bachoc, and Neuvial, 2024;
 * it takes a square matrix containing a simulated white noise 2D-sheet,
 * and computes from it one trajectory of the associated 2D POSIR process.
 *
 * @param Z a square matrix containing a simulated white noise 2D-sheet.
 * @inheritParams n_traj_simu_1D_C
 *
 * See n_traj_simu_1D_C() above for further details.
 */
NumericVector n_traj_simu_2D_C(arma::mat const & Z, arma::ivec const & Intdgrid,
                           bool is_standard);

/// Compute a batch of trajectories of the temporal d-dimensional POSIRd process
/**
 * This is the multivariate extension of n_traj_simu_1D_C().
 * Here Z is a cube with dimensions:
 * trajectories x time x signal_dimension.
 *
 * @param Z cube of simulated standard white-noise trajectories.
 * @param Intdgrid decreasing vector of integer minimal interval lengths.
 *
 * @return A matrix with one row per trajectory and one column per delta.
 */
arma::mat n_traj_simu_dD_C(arma::cube const & Z, arma::ivec const & Intdgrid);





#endif // C_N_TRAJ_SIMU_H
