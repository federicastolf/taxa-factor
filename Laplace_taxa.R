
Ltaxa_int = function(Y, taxonomy, param){
  L = ncol(taxonomy)
  taxa_names = names(taxonomy)
  results = vector("list", L) # each elements is a matrix 2xp_l
  # data aggregated by taxonomic level
  for(l in 1:L){
    Yl = data_tree(Y, taxonomy[,taxa_names[l]])
    if(l==1){
      alphal = rep(0, ncol(Yl))
      # gamma_samp = emp_bayes(as.matrix(Yl), alphal, param$a_gamma, param$b_gamma, 
      #                        param$max_it, param$epsilon, param$Niter, 
      #                        param$burnin, param$eps_MH)
      # # apply Laplace on Yl
      # gamma_hat = mean(exp(gamma_samp))
      gamma_hat = mean(rowSums(Yl))
      results[[l]] = marginal_probit(as.matrix(Yl), gamma_hat, alphal, param$epsilon, param$max_it)
    }
    else{
      lb = unique(taxonomy[,taxa_names[l-1]])
      #alphal = rep(0, ncol(Yl))
      alphal = NULL
      for(i in 1:length(lb)){
        t1 = taxonomy[taxonomy[,l-1] == lb[i],l]
        nchilds = length(unique(t1))
        alphal = c(alphal, rep(results[[l-1]][1,i], nchilds))
      }
      # apply Laplace on Yl
      # gamma_samp = emp_bayes(as.matrix(Yl), alphal, param$a_gamma, param$b_gamma, 
      #                        param$max_it, param$epsilon, param$Niter, 
      #                        param$burnin, param$eps_MH)
      # gamma_hat = mean(exp(gamma_samp))
      gamma_hat = mean(rowSums(Yl))
      results[[l]] = marginal_probit(as.matrix(Yl), gamma_hat, alphal, param$epsilon, param$max_it)
    }
    cat("layer",l,"\n")
  }
  return(results)
}

# drop observations with only one children
# we initially drop the observations that have only one-children structure from 
# phylum to species (155)
# then starting from l=2 we drop any observation which is only-child
############## N.B need to be fixed ######
Ltaxa_int1 = function(Y, taxonomy, param){
  L = ncol(taxonomy)
  taxa_names = names(taxonomy)
  results = vector("list", L) # each elements is a matrix 2xp_l
  # data aggregated by taxonomic level
  for(l in 1:L){
    if(l==1){
      Yl = data_tree(Y, taxonomy[,taxa_names[l]])
      alphal = rep(0, ncol(Yl))
      # gamma_samp = emp_bayes(as.matrix(Yl), alphal, param$a_gamma, param$b_gamma, 
      #                        param$max_it, param$epsilon, param$Niter, 
      #                        param$burnin, param$eps_MH)
      # gamma_hat = mean(exp(gamma_samp))
      gamma_hat = mean(rowSums(Yl))
      results[[l]] = marginal_probit(as.matrix(Yl), gamma_hat, alphal, param$epsilon, param$max_it)
      taxlab_l1 = colnames(Yl)
    }
    else{
      Yl_all =  dad_tree(Y, taxonomy, l)
      Yl = Yl_all$data
      taxa_update = taxonomy[!(taxonomy[,l] %in% Yl_all$taxa_onechild),]
      # cava tipo taxa drop da taxonomy
      # sistema for l=2
      alphal = NULL
      for(i in 1:length(taxlab_l1)){
        t1 =  taxa_update[taxa_update[,l-1] == taxlab_l1[i],l]
        nchilds = length(unique(t1))
        alphal = c(alphal, rep(results[[l-1]][1,i], nchilds))
      }
      # lb = unique(taxonomy[,taxa_names[l-1]])
      # alphal = NULL
      # # qua c'è l'errore!
      # for(i in 1:length(lb)){
      #   t1 = taxonomy[taxonomy[,l-1] == lb[i],l]
      #   nchilds = length(unique(t1))
      #   alphal = c(alphal, rep(results[[l-1]][1,i], nchilds))
      # }
      # apply Laplace on Yl
      # gamma_samp = emp_bayes(as.matrix(Yl), alphal, param$a_gamma, param$b_gamma, 
      #                        param$max_it, param$epsilon, param$Niter, 
      #                        param$burnin, param$eps_MH)
      # gamma_hat = mean(exp(gamma_samp))
      gamma_hat = mean(rowSums(Yl))
      results[[l]] = marginal_probit(as.matrix(Yl), gamma_hat, alphal, param$epsilon, param$max_it)
      taxlab_l1 = colnames(Yl)
    }
    cat("layer",l,"\n")
  }
  return(results)
}


Ltaxa_int = function(Y, taxonomy, param){
  L = ncol(taxonomy)
  taxa_names = names(taxonomy)
  results = vector("list", L) # each elements is a matrix 2xp_l
  # data aggregated by taxonomic level
  for(l in 1:L){
    Yl = data_tree(Y, taxonomy[,taxa_names[l]])
    if(l==1){
      alphal = rep(0, ncol(Yl))
      # gamma_samp = emp_bayes(as.matrix(Yl), alphal, param$a_gamma, param$b_gamma, 
      #                        param$max_it, param$epsilon, param$Niter, 
      #                        param$burnin, param$eps_MH)
      # # apply Laplace on Yl
      # gamma_hat = mean(exp(gamma_samp))
      gamma_hat = mean(rowSums(Yl))
      results[[l]] = marginal_probit(as.matrix(Yl), gamma_hat, alphal, param$epsilon, param$max_it)
    }
    else{
      lb = unique(taxonomy[,taxa_names[l-1]])
      #alphal = rep(0, ncol(Yl))
      alphal = NULL
      for(i in 1:length(lb)){
        t1 = taxonomy[taxonomy[,l-1] == lb[i],l]
        nchilds = length(unique(t1))
        alphal = c(alphal, rep(results[[l-1]][1,i], nchilds))
      }
      # apply Laplace on Yl
      # gamma_samp = emp_bayes(as.matrix(Yl), alphal, param$a_gamma, param$b_gamma, 
      #                        param$max_it, param$epsilon, param$Niter, 
      #                        param$burnin, param$eps_MH)
      # gamma_hat = mean(exp(gamma_samp))
      gamma_hat = mean(rowSums(Yl))
      results[[l]] = marginal_probit(as.matrix(Yl), gamma_hat, alphal, param$epsilon, param$max_it)
    }
    cat("layer",l,"\n")
  }
  return(results)
}