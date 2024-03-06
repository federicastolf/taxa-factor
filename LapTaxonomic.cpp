
#include "RcppArmadillo.h"
#include "Rcpp.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


arma::vec vec_beta_d1d2(arma::vec y, double beta, double prior_var, double prior_mean){
  int n = y.size();
  arma::vec q = 2*y - arma::ones(n); 
  arma::vec lambda(n, arma::fill::zeros);
  for(int i=0; i<n; ++i){
    lambda(i) = q(i)*arma::normpdf(q(i)*beta)/arma::normcdf(q(i)*beta);
  }
  double del_log_posterior = sum(lambda); 
  double score = del_log_posterior - (beta - prior_mean)/prior_var; 
  
  double hessian=0;
  for(int i=0; i<n; ++i){
    hessian = hessian + lambda(i)*(beta + lambda(i)); 
  }
  hessian = hessian + 1/prior_var;
  arma::vec result(2);
  result(0) = score;
  result(1) = -hessian;
  return result;
}


arma::vec vec_beta_laplace(arma::vec y, double gamma, double alpha_l1j, int p, int max_it, double epsilon){
  double prior_var = 2*log(p);
  double prior_mean =  alpha_l1j + sqrt(1 + prior_var) * R::qnorm(gamma / (gamma + p), 0.0, 1.0, 1, 0);
  double b = 0; 
  double H = 0;
  for(int i = 1; i<=max_it; ++i){
    arma::vec res = vec_beta_d1d2(y, b, prior_var, prior_mean);
    double u = res(0);
    H = res(1);
    double err = -u/H; 
    if(err<epsilon){
      break;
    }
    b = b - u/H;
  }
  double vec_size = 2; 
  arma::vec result(vec_size);
  result(0) = b;
  double H_inv = -1/H; 
  result(1) = H_inv;
  return result;
}



// [[Rcpp::export]]

arma::mat marginal_probit(arma::mat Y, double gamma, arma::vec alpha_l1, double epsilon, int max_it){
  int n = Y.n_rows;
  int p = Y.n_cols;
  double vec_size = 2; 
  arma::mat all_res(2, p);
  for(int j=0; j<p; ++j){
    all_res(arma::span(0, vec_size - 1), j) = vec_beta_laplace(Y(arma::span(0,n-1), j), gamma, alpha_l1(j), p, max_it, epsilon); 
  }
  return all_res;
}

// log lik for gamma for MH
double log_lik_gamma(double gamma, arma::vec beta){
  int p = beta.size();
  double tau2 = 2 * log(p);
  double mu = sqrt(1 + tau2) * R::qnorm(gamma / (gamma + p), 0.0, 1.0, 1, 0);
  if (gamma > 0){
    arma::vec logl = log_normpdf(beta, mu, sqrt(tau2));
    return sum(logl);
  }
  else{
    return R_NegInf;
  }
  
}

// Log-posterior for MH
double logpost_gamma(double gamma, arma::vec beta, double a_gamma, double b_gamma) {
  double log_likelihood = log_lik_gamma(gamma, beta);
  return log_likelihood + R::dgamma(exp(gamma), a_gamma, b_gamma, true) + gamma;
}


// [[Rcpp::export]]

arma::vec emp_bayes(arma::mat Y, arma::vec alpha_l1, double a_gamma, double b_gamma, int max_it, double epsilon, double nmcmc, double burnin, double eps_MH){
  int p = Y.n_cols; 
  arma::vec gamma_samples((nmcmc - burnin));
  // MH initializiation
  double gamma = arma::randg(1)(0);
  
  for(int ii =1; ii<=nmcmc; ++ii){
    //first stage
    arma::mat beta_res = marginal_probit(Y, gamma, alpha_l1, epsilon, max_it);
    arma::vec beta_mean = beta_res.row(0).t();
    arma::vec var_beta = beta_res.row(1).t();
    arma::vec epsilon = arma::randn<arma::vec>(p);
    arma::vec b_samples = beta_mean + var_beta % epsilon;
    // MH step
    double logp = logpost_gamma(gamma, b_samples, a_gamma, b_gamma);
    double gamma_new = arma::randn() * eps_MH + gamma;
    double logp_new = logpost_gamma(gamma_new, b_samples, a_gamma, b_gamma);
    double a_acc = std::min(1.0, exp(logp_new - logp));
    if (R::runif(0.0, 1.0) < a_acc){
      logp = logp_new;
      gamma = gamma_new;
    }
    if(ii> burnin)
    {
      gamma_samples(ii - burnin -1) = gamma;
    }
  }
  //double gamma_hat = mean(gamma_samples);
  return gamma_samples;
}