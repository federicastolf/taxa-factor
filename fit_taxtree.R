
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(pROC)


load("/Users/stolffederica/Library/CloudStorage/Dropbox/dark taxa/code/fungi_data/allData_clim_trait.Rdata")
source("function_fungi.R")
# source("Laplace_taxa.R")
#sourceCpp("LapTaxonomic.cpp") # only-intercept model
sourceCpp("LapTaxonomic_cov.cpp") # probit regregression model
#sourceCpp("bigMVPmarginal.cpp") # bigMVP marginal

####---------------------#### prepare the data ####----------------####
taxonomy = taxonomy[,4:9]

## drop the columns that have only one child structure from phylum to species 
phylan = as.character(taxonomy$phylum)
phyla_cut = rep(0,155)
j=1
for(i in 1:length(phylan)){
  t1 = taxonomy %>% filter(phylum==phylan[i])
  if(nrow(t1)==1){
    phyla_cut[j] = phylan[i]
    j=j+1
  }
}
idxcut = which(taxonomy$phylum %in% phyla_cut)
fungi = otu.table[,-idxcut]
fungi = matrix(as.numeric(fungi >0), nrow(fungi), ncol(fungi))
taxonomy = taxonomy[-idxcut,]

## define model matrix X
Xfungi = data.frame(intercept = rep(1, nrow(fungi)), temperature = meta$temp.mean,
                    temperature2 = meta$temp.mean^2,
                    seqdepth = log(meta$numspikes+meta$numnonspikes),
                    lat = meta$lat/100, seasonality1 = sin(2*pi*meta$yday/365),
                    seasonality2 = cos(2*pi*meta$yday/365))
Xfungi = cbind.data.frame(Xfungi, interaction1 = Xfungi$lat*Xfungi$seasonality1, 
                          interaction2 = Xfungi$lat*Xfungi$seasonality2)

## drop rows with all 0
# idxr = which(rowSums(fungi)==0)
# fungi=fungi[-idxr,]
# # in case you have to do also for covariate! 
# Xfungi = Xfungi[-idxr,]

## drop rows with seqdepth<10000
idxs = which(exp(Xfungi$seqdepth)<10000)
Xfungi = Xfungi[-idxs,]
fungi = fungi[-idxs,]

###--------------------### fit the models ###-----------------------###
q = ncol(Xfungi)
param = list(max_it = 100, epsilon = 0.0001,
             a0 = 1/2, b0 = 2, nu0 = rep(0, q-1), Psi0=diag(q-1), atheta = 3, 
             btheta = 2,  burnin = 0, nmcmc = 10)
# param = list(max_it = 100, epsilon = 0.0001, burnin=500, Niter=2000, 
#              eps_MH = 0.05, a_gamma=10, b_gamma=1/4, prior_var=15)

fitcov = Ltaxa(fungi, as.matrix(Xfungi), taxonomy, param)

#fit_taxonomic = Ltaxa1_int(fungi, taxonomy, param)
fit_taxonomic = Ltaxa_int(fungi, taxonomy, param)


## prove 
res1 = marginal_probit(as.matrix(Yl), as.matrix(Xfungi), gamma_hat, alphal, betal, 
                       param$prior_var, param$epsilon, param$max_it)





###----------------### validate results - AUC ###-------------------###
prob_fit = pnorm(fit_taxonomic[[6]][1,])
auc_all= rep(0, ncol(fungi))
for(i in 1:length(auc_all)){
  auc_all[i] = auc(fungi[,i], rep(prob_fit[i], nrow(fungi)))
}
summary(auc_all)
auc(fungi[1:10,i], rep(prob_fit[i], 10))
