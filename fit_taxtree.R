
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(tidyverse)

# load("data_try/tree_fungi.Rdata")
# load("data_try/fungi_binary.Rdata")
load("/Users/stolffederica/Library/CloudStorage/Dropbox/dark taxa/code/fungi_data/allData_clim_trait.Rdata")
source("function_fungi.R")
source("Laplace_taxa.R")
sourceCpp("LapTaxonomic.cpp")

# tree_fungi = tree_fungi[,2:7]
taxonomy = taxonomy[,4:9]
# nunk = nrow(tree_fungi[tree_fungi$Species=="unk",])
# lab = rep(0,nunk)
# for(i in 1:length(lab)){
#   lab[i] = paste0("unk",i) 
# }
# idxunk = which(tree_fungi$Species=="unk")
# tree_fungi$Species[idxunk] = lab

#--# drop the columns that have only one child structure from phylum to species #--#
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

param = list(max_it = 100, epsilon = 0.0001, burnin=500, Niter=2000, 
             eps_MH = 0.03, a_gamma=10, b_gamma=1/4)
# param = list(max_it = 100, epsilon = 0.0001)

fit_taxonomic = Ltaxa1(fungi, tree_fungi, param)
fit_taxonomic = Ltaxa(fungi, tree_fungi, param)
