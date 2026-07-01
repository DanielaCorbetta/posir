#include "C_n_traj_simu.h"

// [[Rcpp::depends(RcppArmadillo)]]
/* // [[Rcpp::plugins("cpp11")]] */

/**
 * C++ doc is in .h files.
 *
 * Here is only the R doc processed by Roxygen.
 */

/**
 * R doc for Roxygen
 */
//' Log and raise an exception (stop) from C++
//'
//' A message is logged at level ERROR, then stop(message) is called.
//' An optional additional message (details) can be given,
//' to be logged at level DEBUG.
//'
//' @inheritParams log_n_stop
//' @family log and stop from R or C++
//'
//' @examples
//' \dontrun{
//'   tryCatch(log_n_stop_C("essai de log_n_top_C", "via tryCatch()"),
//'            error= function(e){log_info(paste("erreur :",e))},
//'            finally=function(e){log_info(paste("finally :",e))})
//' }
//' @keywords internal
// [[Rcpp::export]]
void log_n_stop_C(std::string message, std::string m_details) {
  Function log_er("log_error"), log_de("log_debug");
  log_er(message);
  if(!m_details.empty()) log_de(m_details);
  stop(message) ;
  return;
}

arma::mat RowMult(arma::mat const & M, arma::vec const & v) {
  arma::uword k = M.n_rows, l = M.n_cols, i, j ;
  arma::mat res(k,l);
  if(k!=v.n_elem)
    log_n_stop_C("incompatible dimensions", "in posir/Rcpp RowMult().");
  for(i=0; i<k; i++) for(j=0; j<l; j++) res(i,j)=M(i,j)*v[i];
  return res;
}

arma::mat RowDiv(arma::mat const & M, arma::vec const & v) {
  arma::uword i, l=v.n_elem;
  arma::vec vinv(l) ;
  for(i=0; i<l; i++) vinv(i)=1/v(i);
  return RowMult(M, vinv);
}

/**
 * R doc for Roxygen
 */
//' Computation of a batch of trajectories of the 1D POSIR process
//'
//' Implements the functional \eqn{G_\delta(\cdot)}
//' from paper \insertCite{BaBoNe24}{posir}.
//' Takes a matrix whose lines are from simulated white noise trajectories
//' and computes from it the trajectories of the associated 1D POSIR process.
//'
//' @param Z matrix whose lines are simulated white noise trajectories.
//' @param Intdgrid decreasing vector of the integers \eqn{1\le k\le Ndis}
//' such that \eqn{k/Ndis} is in the grid for \eqn{\delta},
//' where Ndis is the rows number of Z.
//' Generally obtained as Ndis*deltagrid,
//' where the vector deltagrid is the desired grid for \eqn{\delta}.
//' (Beware: no checking is done in this internal function.)
//' @inheritParams n_traj_simu
//' @inheritParams simulationDir
//' @inheritParams check_grid
//' @seealso [run_simu()], [batch_filename()]
//' @family POSIR process generators
//'
//' @return A matrix with n lines and length(Intdgrid) columns.
//'
//' @references
//'   \insertAllCited{}
//' @keywords internal
// [[Rcpp::export]]
arma::mat n_traj_simu_1D_C(arma::mat const & Z, arma::ivec const & Intdgrid,
                           bool is_standard) {
  int n=Z.n_rows, Ndis=Z.n_cols, l=Intdgrid.n_elem, k, i, pos=0 ;
  int q=Intdgrid[l-1] ;
  arma::mat Y(n, l), X(n, Ndis+1) ;
  arma::colvec Ycur(n), Yprec(n) ;
  if(!is_standard){
    Ycur = stddev(Z,0,1);
    X.submat(0, 1, n-1, Ndis) = cumsum(RowDiv(Z, Ycur), 1);
  }
  else{
    X.submat(0, 1, n-1, Ndis) = cumsum(Z, 1);
  }
  X.submat(0, 0, n-1, 0).zeros();
  Yprec.zeros();
  for(k=Ndis ; k>=q ; k--){
    Ycur=max(abs(X.submat(0, k, n-1, Ndis)-X.submat(0, 0, n-1, Ndis-k))/sqrt(k),1);
    for(i=0;i<n;i++) {
      Ycur[i]=std::max(Ycur[i], Yprec[i]);
    }
    if(k==Intdgrid[pos]) {
      Y.submat(0, pos, n-1, pos) = Ycur;
      pos++;
    }
    Yprec=Ycur;
  }
  return Y ;
}

/**
 * R doc for Roxygen
 */
//' Computation of a batch of trajectories of the temporal d-dimensional POSIRd process
 //'
 //' Multivariate extension of [n_traj_simu_1D_C].
 //' Takes a cube whose dimensions are:
 //' trajectory x time x signal dimension.
 //'
 //' @param Z cube of simulated standard white-noise trajectories.
 //' @param Intdgrid decreasing vector of integer minimal interval lengths.
 //'
 //' @return A matrix with n lines and length(Intdgrid) columns.
 //'
 //' @keywords internal
 // [[Rcpp::export]]
 arma::mat n_traj_simu_dD_C(arma::cube const & Z,
                            arma::ivec const & Intdgrid) {

   int n = Z.n_rows;
   int Ndis = Z.n_cols;
   int d = Z.n_slices;
   int l = Intdgrid.n_elem;
   int pos = 0;
   int q = Intdgrid[l - 1];

   arma::mat Y(n, l);
   arma::cube X(n, Ndis + 1, d);
   arma::colvec Ycur(n), Yprec(n);

   // X[i, k, h] = cumulative sum up to time k
   // for trajectory i and component h.
   for(int h = 0; h < d; h++) {
     X.slice(h).col(0).zeros();
     X.slice(h).cols(1, Ndis) = cumsum(Z.slice(h), 1);
   }

   Yprec.zeros();

   for(int k = Ndis; k >= q; k--) {

     // sqnorms[i, a] will contain
     // || S_{a+k} - S_a ||_2^2
     // for trajectory i and interval start a.
     arma::mat sqnorms(n, Ndis - k + 1, arma::fill::zeros);

     for(int h = 0; h < d; h++) {
       arma::mat inc =
         X.slice(h).submat(0, k, n - 1, Ndis) -
         X.slice(h).submat(0, 0, n - 1, Ndis - k);

       sqnorms += inc % inc;
     }

     // For each trajectory, take the maximum over intervals of length k.
     Ycur = max(sqrt(sqnorms) / sqrt((double) k), 1);

     // Compare with the previous supremum over longer intervals.
     for(int i = 0; i < n; i++) {
       Ycur[i] = std::max(Ycur[i], Yprec[i]);
     }

     // Save if k corresponds to one requested delta.
     if(k == Intdgrid[pos]) {
       Y.col(pos) = Ycur;
       pos++;
     }

     Yprec = Ycur;
   }

   return Y;
 }

double max2Dwindow(arma::mat const & cummat, int N, int i, int j) {
  arma::mat dY ;
  dY = abs(cummat.submat(i,j,N,N)
             + cummat.submat(0,0,N-i,N-j)
             - cummat.submat(i,0,N,N-j)
             - cummat.submat(0,j,N-i,N)
             )/sqrt(i*j) ;
  dY.reshape((N+1-i)*(N+1-j),1);
  return as_scalar(max(dY)) ;
}

/**
 * R doc for Roxygen
 */
//' Computation one trajectory of the 2D POSIR process
//'
//' Implements the 2D version of the functional \eqn{G_\delta(\cdot)}
//' from paper \insertCite{BaBoNe24}{posir}.
//' Takes a square matrix containing a simulated white noise 2D-sheet,
//' and computes from it one trajectory of the associated 2D POSIR process.
//'
//' @param Z square matrix containing a simulated white noise 2D-sheet.
//' @inheritParams n_traj_simu_1D_C
//' @inheritParams n_traj_simu
//' @inheritParams simulationDir
//' @inheritParams check_grid
//' @seealso [run_simu()], [batch_filename()]
//' @family POSIR process generators
//'
//' @return A vector with length(Intdgrid) elements.
//'
//' @references
//'   \insertAllCited{}
//' @keywords internal
// [[Rcpp::export]]
NumericVector n_traj_simu_2D_C(arma::mat const & Z,
                               arma::ivec const & Intdgrid,
                               bool is_standard) {
  if(Z.n_rows!=Z.n_cols) {
    log_n_stop_C("n_traj_simu_2D_C(): Z is not a square matrix", "");
  }
  int l=Intdgrid.n_elem, k, i, j, Ndis=Z.n_rows;
  int prev_pos=Ndis, next_pos ;
  double s;
  NumericVector res(l);
  arma::mat Y ;
  if(!is_standard) {
    Y = Z;
    Y.reshape(Ndis*Ndis,1);
    s = as_scalar(stddev(Y,0,0)) ;
  }
  Y.set_size(Ndis+1, Ndis+1);
  Y.submat(0, 0, 0, Ndis).zeros();
  Y.submat(1, 0, Ndis, 0).zeros();
  Y.submat(1,1, Ndis, Ndis) = cumsum(cumsum(Z, 0), 1);
  if(!is_standard) {
    Y = Y/s;
  }
  s = std::abs(Y(Ndis,Ndis))/Ndis ;
  for(i=0; i<l; i++) {
    next_pos=Intdgrid[i];
    for(j=prev_pos-1; j>=next_pos; j--){
      for(k=Ndis; k>=next_pos; k--) {
        s = std::max(s, max2Dwindow(Y,Ndis,j,k)) ;
      }
    }
    for(j=Ndis; j>=prev_pos; j--){
      for(k=prev_pos-1; k>=next_pos; k--) {
        s = std::max(s, max2Dwindow(Y,Ndis,j,k)) ;
      }
    }
    res(i) = s ;
    prev_pos=next_pos ;
  }
  return res;
}

