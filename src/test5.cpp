#include <iostream>
#include <RcppArmadillo.h>
#include "FastProxSL1.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//` ADMM to solve SLOPE-penalized ordinary least squares regression
//` 
//` @param X matrix of independent  variables
//` @param Y vector of dependent variables
//` @param max_iter maximum number of iterations of the ADMM loop
//` @param rho augmented Lagrangian parameter
//` @export
//` @examples
//` library(SLOPE)
//` fit <- ADMM(as.matrix(bodyfat$x), bodyfat$y, 1e4, 1)

//' @title ADMM
//' @description
//' ADMM to solve SLOPE-penalized ordinary least squares regression.
//' @name ADMM
//' @param X matrix of independent  variables
//' @param Y vector of dependent variables
//' @param max_iter maximum number of iterations of the ADMM loop
//' @param rho augmented Lagrangian parameter
//' @examples
//' library(SLOPE)
//' fit <- ADMM(as.matrix(bodyfat$x), bodyfat$y, 1e4, 1)
//' 
//' @export
//' 
//' 

// [[Rcpp::export]]
arma::vec ADMM(arma::mat X, arma::vec Y, arma::vec lambda, arma::uword max_iter, double rho){
  
  uword p = X.n_cols;
  
  vec beta(p, fill::zeros), z(p, fill::zeros), y(p, fill::zeros);
  mat I(p,p, fill::eye);
  z = beta;
  y = beta;
  
  mat matrix = inv(X.t()*X + rho*I);
  mat y_tilde = X.t()*Y;
  
  
  for(uword i=0; i<max_iter; i++) {
    
    beta = matrix*(y_tilde + rho*(z-y));
    z = beta + y;
    vec z_sign = sign(z);
    z = abs(z);
    uvec z_order = sort_index(z, "descend");
    z = (z(z_order)).eval();
    z = FastProxSL1(z, lambda/rho);
    z(z_order) = z;
    z %= z_sign;
    
    y = y + beta - z;
    
  }
  
  return beta;
}