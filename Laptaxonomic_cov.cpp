
#include "RcppArmadillo.h"
#include "Rcpp.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

arma::mat vec_beta_d1d2(arma::vec y, arma::mat X, arma::vec beta, arma::vec mean_b, arma::mat inv_covb){
  int n = y.size();
  int q = X.n_cols;
  arma::vec qr = 2*y - arma::ones(n);
  arma::vec lambda(n, arma::fill::zeros);
  for(int i=0; i<n; ++i){
    lambda(i) = qr(i)*arma::normpdf(qr(i)*sum(X(i,arma::span(0,q-1))*beta))/arma::normcdf(qr(i)*sum(X(i,arma::span(0,q-1))*beta));
  }
  arma::vec del_log_posterior(q, arma::fill::zeros);
  for(int j=0; j<q; ++j){
    del_log_posterior(j) = sum(lambda.t()*X(arma::span(0,n-1), j));
  }
  arma::vec score = del_log_posterior - inv_covb*(beta - mean_b);
  
  arma::mat hessian(q, q, arma::fill::zeros);
  for(int i=0; i<n; ++i){
    hessian = hessian + lambda(i)*((sum(X(i,arma::span(0,q-1))*beta)) + lambda(i))*X(i, arma::span(0, q-1)).t()*X(i, arma::span(0,q-1));
  }
  hessian = hessian + inv_covb;
  arma::mat result(q, q+1);
  result(arma::span(0, q-1), 0) = score;
  result(arma::span(0, q -1), arma::span(1, q)) = -hessian;
  return result;
}


arma::vec vec_beta_laplace(arma::vec y, arma::mat X, arma::vec mean_b, arma::mat inv_covb, int max_it, double epsilon){
  
  int q = X.n_cols;
  arma::vec b(q, arma::fill::zeros);
  arma::mat H(q, q, arma::fill::zeros);
  for(int i = 1; i<=max_it; ++i){
    arma::mat res = vec_beta_d1d2(y, X, b, mean_b, inv_covb);
    arma::vec u = res(arma::span(0, q-1), 0);
    H = res(arma::span(0, q-1), arma::span(1, q));
    arma::vec err = -inv(H)*u;
    double norm_err = norm(err, 2);
    if(norm_err<epsilon){
      break;
    }
    b = b - inv(H)*u;
  }

  double vec_size = q + q*q;
  arma::vec result(vec_size);
  result(arma::span(0, q-1)) = b;
  arma::mat I = arma::diagmat(arma::ones(q));
  arma::mat H_inv = solve(-H, I);
  arma::vec cov_vec = arma::vectorise(H_inv);
  result(arma::span(q, q-1 + pow(q, 2))) = cov_vec;
  return result;
}



// [[Rcpp::export]]

arma::mat marginal_probit(arma::mat Y, arma::mat X, double gamma, arma::mat mean_b, double prior_var, double epsilon, int max_it){
  int n = Y.n_rows;
  int p = Y.n_cols;
  int q = X.n_cols;
  // compute tau_p and mu_p
  double taup = 2*log(p);
  double mup = sqrt(1+taup)*R::qnorm(gamma / (gamma + p), 0.0, 1.0, 1, 0);
  // arma::mat mean_b of dimension qxp
  mean_b.row(0) += arma::ones<arma::rowvec>(p) * mup;
  // prior_var and inv_covb need to be fixed
  arma::vec var_b(q);
  var_b(0) = taup;
  var_b.subvec(1, q-1).fill(prior_var);
  arma::mat inv_covb = diagmat(1/var_b);

  double vec_size = q + q*q;
  arma::mat all_res(vec_size, p);
  for(int j=0; j<p; ++j){
    all_res(arma::span(0, vec_size - 1), j) = vec_beta_laplace(Y(arma::span(0,n-1), j), X, mean_b(arma::span(0,q-1), j), inv_covb, max_it, epsilon);   

  }
  return all_res;
}

