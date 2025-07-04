// matrix_multiply.cpp

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat matrix_multiplication_cpp(const arma::mat& x,
                                    const arma::mat& y,
                                    bool tx = false,
                                    bool ty = false,
                                    unsigned int cores = 1) {
  int n = 0, p = 0;

  if (tx) {
    n = x.n_cols;
    p = y.n_cols;
  } else {
    n = x.n_rows;
    p = ty ? y.n_rows : y.n_cols;
  }

  arma::mat C(n, p, fill::zeros);
  arma::colvec yi;

  if (!tx && !ty) {
    arma::mat xt = x.t();
    for (int i = 0; i < p; ++i) {
      yi = y.col(i);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
      for (int j = 0; j < n; ++j) {
        C(j, i) = dot(xt.col(j), yi);
      }
    }
  } else if (tx) {
    for (int i = 0; i < p; ++i) {
      yi = y.col(i);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
      for (int j = 0; j < n; ++j) {
        C(j, i) = dot(x.col(j), yi);
      }
    }
  } else { // only ty = true
    arma::mat yt = y.t();
    arma::mat xt = x.t();
    for (int i = 0; i < p; ++i) {
      yi = yt.col(i);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
      for (int j = 0; j < n; ++j) {
        C(j, i) = dot(xt.col(j), yi);
      }
    }
  }

  return C;
}

