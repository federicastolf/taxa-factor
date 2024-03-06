
Ltaxa = function(Y, taxonomy, param){
  L = ncol(taxonomy)
  taxa_names = names(taxonomy)
  results = vector("list", L) # each elements is a matrix 2xp_l
  # data aggregated by taxonomic level
  for(l in 1:L){
    Yl = data_tree(Y, taxonomy[,taxa_names[l]])
    if(l==1){
      alphal = rep(0, ncol(Yl))
      gamma_samp = emp_bayes(as.matrix(Yl), alphal, param$a_gamma, param$b_gamma, 
                             param$max_it, param$epsilon, param$Niter, 
                             param$burnin, param$eps_MH)
      # apply Laplace on Yl
      gamma_hat = mean(exp(gamma_samp))
      results[[l]] = marginal_probit(as.matrix(Yl), gamma_hat, alphal, param$epsilon, param$max_it)
    }
    else{
      lb = unique(taxonomy[,taxa_names[l-1]])
      alphal = NULL
      for(i in 1:length(lb)){
        t1 = taxonomy[taxonomy[,l-1] == lb[i],l]
        nchilds = length(unique(t1))
        alphal = c(alphal, rep(results[[l-1]][1,i], nchilds))
      }
      # apply Laplace on Yl
      gamma_samp = emp_bayes(as.matrix(Yl), alphal, param$a_gamma, param$b_gamma, 
                             param$max_it, param$epsilon, param$Niter, 
                             param$burnin, param$eps_MH)
      gamma_hat = mean(exp(gamma_samp))
      results[[l]] = marginal_probit(as.matrix(Yl), gamma_hat, alphal, param$epsilon, param$max_it)
    }
  }
  return(results)
}