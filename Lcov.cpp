#include "RcppArmadillo.h"
#include "Rcpp.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


arma::mat vec_log_post_beta_d1_d2_nh(arma::vec y, arma::mat X, arma::vec beta, double prior_var){
  int n = y.size();
  int p = X.n_cols;
  arma::vec q = 2*y - arma::ones(n);
  arma::vec lambda(n, arma::fill::zeros);
  for(int i=0; i<n; ++i){
    lambda(i) = q(i)*arma::normpdf(q(i)*sum(X(i,arma::span(0,p-1))*beta))/arma::normcdf(q(i)*sum(X(i,arma::span(0,p-1))*beta));
  }
  arma::vec del_log_posterior(p, arma::fill::zeros);
  for(int j=0; j<p; ++j){
    del_log_posterior(j) = sum(lambda.t()*X(arma::span(0,n-1), j));
  }
  arma::vec score = del_log_posterior - beta/prior_var;
  
  arma::mat hessian(p, p, arma::fill::zeros);
  for(int i=0; i<n; ++i){
    hessian = hessian + lambda(i)*((sum(X(i,arma::span(0,p-1))*beta)) + lambda(i))*X(i, arma::span(0, p-1)).t()*X(i, arma::span(0,p-1));
  }
  arma::mat I = arma::diagmat(arma::ones(p));
  hessian = hessian + (1/prior_var)*I;
  arma::mat result(p, p+1);
  result(arma::span(0, p-1), 0) = score;
  result(arma::span(0, p -1), arma::span(1, p)) = -hessian;
  return result;
}



arma::vec vec_log_post_beta_laplace_nh(arma::vec y, arma::mat X, double prior_var, int max_it, double epsilon){
  
  // change so that it returns a vector
  int n = y.size();
  arma::mat X_new = join_rows(arma::ones(n, 1), X);
  int p = X.n_cols;
  arma::vec b(p+1, arma::fill::zeros);
  arma::mat H(p+1, p+1, arma::fill::zeros);
  for(int i = 1; i<=max_it; ++i){
    arma::mat res = vec_log_post_beta_d1_d2_nh(y, X_new, b, prior_var);
    arma::vec u = res(arma::span(0, p), 0);
    H = res(arma::span(0, p), arma::span(1, p+1));
    arma::vec err = -inv(H)*u;
    double norm_err = norm(err, 2);
    if(norm_err<epsilon){
      break;
    }
    b = b - inv(H)*u;
  }
  double vec_size = (p+1) + (p+1)*(p+1);
  arma::vec result(vec_size);
  result(arma::span(0, p)) = b;
  arma::mat I = arma::diagmat(arma::ones(p+1));
  arma::mat H_inv = solve(-H, I);
  arma::vec cov_vec = arma::vectorise(H_inv);
  result(arma::span(p+1, p + pow(p+1, 2))) = cov_vec;
  return result;
}

// [[Rcpp::export]]

arma::mat marginal_prv(arma::mat Y, arma::mat X, double prior_var, double epsilon, int max_it){
  int n = Y.n_rows;
  int q = Y.n_cols;
  int p = X.n_cols;
  double vec_size = (p+1) + (p+1)*(p+1);
  arma::mat all_res(vec_size, q);
  //omp_set_dynamic(0);
  //omp_set_num_threads(4);
  //#pragma omp parallel for 
  for(int j=0; j<q; ++j){
    all_res(arma::span(0, vec_size - 1), j) = vec_log_post_beta_laplace_nh(Y(arma::span(0,n-1), j), X, prior_var, max_it, epsilon); 
  }
  return all_res;
}