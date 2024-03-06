
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(tidyverse)

load("data_try/tree_fungi.Rdata")
load("data_try/fungi_binary.Rdata")
source("function_fungi.R")
source("Laplace_taxa.R")
sourceCpp("LapTaxonomic.cpp")

tree_fungi = tree_fungi[,2:7]
nunk = nrow(tree_fungi[tree_fungi$Species=="unk",])
lab = rep(0,nunk)
for(i in 1:length(lab)){
  lab[i] = paste0("unk",i) 
}
idxunk = which(tree_fungi$Species=="unk")
tree_fungi$Species[idxunk] = lab

param = list(max_it = 100, epsilon = 0.0001, burnin=500, Niter=2000, 
             eps_MH = 0.03, a_gamma=10, b_gamma=1/4)
# param = list(max_it = 100, epsilon = 0.0001)

fit_taxonomic = Ltaxa(fungi, tree_fungi, param)
