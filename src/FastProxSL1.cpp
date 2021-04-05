#include <iostream>
#include <RcppArmadillo.h>
#include "FastProxSL1.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' @title FastProxSL1
//' @description
//' A function using RcppArmadillo that computes the proximal operator for SLOPE using Algorithm 3 (FastProxSL1) from Bogdan et al 2015 (SLOPE: adaptive variable selection via convex optimization).
//' @name FastProxSL1
//' @param beta input vector
//' @param lambda vector of lambdas, sorted in decreasing order
//' @examples
//' p = 1e7
//` lambda = 0.5 * (p - (1:p)) + 1
//` y = rnorm(p, mean = p/2, sd = 2)
//' FastProxSL1(y, lambda)
//' 
//' @export
//' 
//' 
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::vec FastProxSL1( arma::vec beta,  arma::vec lambda)
{
  arma::uword p = beta.n_elem;
  
  // collect sign of beta and work with sorted absolutes
  vec beta_vec = vectorise(beta);
  vec beta_sign = sign(beta_vec);
  beta_vec = abs(beta_vec);
  uvec beta_order = sort_index(beta_vec, "descend");
  beta_vec = (beta_vec(beta_order)).eval();
  
  vec s(p);
  vec w(p);
  
  uvec idx_i(p);
  uvec idx_j(p);
  
  uword k = 0;
  
  for (uword i = 0; i < p; i++) {
    
    idx_i(k) = i;
    idx_j(k) = i;
    s(k) = beta_vec(i) - lambda(i);
    w(k) = s(k);
    
    while ((k > 0) && (w(k - 1) <= w(k))) {
      k--;
      idx_j(k) = i;
      s(k) += s(k + 1);
      w(k) = s(k) / (i - idx_i(k) + 1.0);
    }
    k++;
  }
  
  for (uword j = 0; j < k; j++) {
    double d = std::max(w(j), 0.0);
    for (uword i = idx_i(j); i <= idx_j(j); i++) {
      beta_vec(i) = d;
    }
  }
  
  // reset order
  beta_vec(beta_order) = beta_vec;
  
  beta_vec %= beta_sign;

  return beta_vec;
}