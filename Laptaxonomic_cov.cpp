
#include "RcppArmadillo.h"
#include "Rcpp.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


arma::mat remove_first_row_col(arma::mat S) {
  int q = S.n_cols;
  arma::mat S1 = S.rows(1, q - 1);
  S1.shed_col(0); // Remove the first column
  return S1;
}

arma::mat gen_mvnrnd_beta(arma::mat M, arma::mat Covs){
  int q = M.n_rows;
  int p = M.n_cols;
  arma::mat X(q,p);
  for(int j=0; j<p; ++j)
  {
    arma::vec m = M(arma::span(0,q-1), j);
    arma::mat S = reshape(Covs(arma::span(0, (q+1)*(q+1) - 1), j), q+1, q+1);
    arma::mat S1 = remove_first_row_col(S);
    X(arma::span(0,q-1),j) = arma::mvnrnd(m, S1);
  }
  return X;
}

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


// [[Rcpp::export]]

List emp_bayes0(arma::mat Y, arma::mat X, double gamma, arma::vec nu0, arma::mat Psi0, double a0, double b0, double epsilon, int max_it, int nmcmc, int burnin){
  int p = Y.n_cols; 
  int q = X.n_cols;
  arma::mat beta0_samples(q-1, (nmcmc - burnin));
  arma::vec sigma0_samples((nmcmc - burnin));
  // initialize
  arma::vec beta0 = arma::zeros(q);
  double sigma0 = 1;
for(int ii =1; ii<=nmcmc; ++ii){
    // obtain approximation of beta^1
    arma::mat beta0mat = arma::repmat(beta0, 1, p);
    arma::mat beta_res = marginal_probit(Y, X, gamma, beta0mat, sigma0, epsilon, max_it);
    arma::mat beta_mean = beta_res(arma::span(1, q-1), arma::span(0, p-1));
    arma::mat cov_mats = beta_res(arma::span(q, q-1 + pow(q, 2)), arma::span(0, p-1));
    //draw beta_j^1

    arma::mat b_samples = gen_mvnrnd_beta(beta_mean, cov_mats);
    // update beta^0
    arma::vec beta_bar = arma::mean(b_samples, 1);
    arma::mat I = diagmat(arma::ones(q-1));
    arma::mat Psi0inv = solve(Psi0, I);    
    arma::mat Sigma1 = arma::diagmat(arma::ones<arma::vec>(q-1) * (p/sigma0)) + Psi0inv;
    arma::mat Sigma1inv = solve(Sigma1, I);   
    arma::vec beta1 = Sigma1inv*(Psi0inv*nu0 + arma::diagmat(arma::ones<arma::vec>(q-1) * (p/sigma0))*beta_bar);
 
    beta0(arma::span(1, q-1)) = arma::mvnrnd(beta1, Sigma1inv);
    //update sigma^0
    double sumbetam = 0;
    for(int j =0; j<p; ++j){
      sumbetam = sumbetam + sum(pow(b_samples(arma::span(0, q-2),j)-beta0(arma::span(1, q-1)),2));
    }
    double sigma0inv = arma::randg(1, distr_param(a0 +((q-1)*p)/2, b0+0.5*sumbetam))(0);
    sigma0 = 1/sigma0inv;
 
  if(ii> burnin)
    {
      beta0_samples(arma::span(0, q-2), (ii - burnin - 1)) = beta0(arma::span(1, q-1));
      sigma0_samples(ii - burnin -1) = sigma0;
    }

    if(ii%20 == 0){
      Rprintf("Iteration = %d \n", ii);
    }
  }

  List result;
  result["beta0"] = beta0_samples;
  result["sigma0"] = sigma0_samples;
  return result;
}




