
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(pROC)

load("/Users/stolffederica/Library/CloudStorage/Dropbox/dark taxa/code/fungi_data/allData_clim_trait.Rdata")
source("function_fungi.R")
sourceCpp("bigMVPmarginal.cpp")

Y = matrix(as.numeric(otu.table >0), nrow(otu.table), ncol(otu.table))

## define model matrix X
Xfungi = data.frame(intercept = rep(1, nrow(Y)), temperature = meta$temp.mean,
                    temperature2 = meta$temp.mean^2,
                    seqdepth = log(meta$numspikes+meta$numnonspikes),
                    lat = meta$lat/100, seasonality1 = sin(2*pi*meta$yday/365),
                    seasonality2 = cos(2*pi*meta$yday/365))
Xfungi = cbind.data.frame(Xfungi, interaction1 = Xfungi$lat*Xfungi$seasonality1, 
                          interaction2 = Xfungi$lat*Xfungi$seasonality2)


r1 = marginal_prv(Y, as.matrix(Xfungi[,-1]), 15, 0.0001, 100) # not working
r1 = marginal_prv(Y, as.matrix(Xfungi[,-c(1,5,6,7)]), 15, 0.0001, 100) # not working

pre = colSums(Y)
sel.sp = pre>=50
Ycommon = Y[,sel.sp]
r1 = marginal_prv(Ycommon, as.matrix(Xfungi[,-1]), 15, 0.0001, 100) # not working


idxs = which(meta$numspikes+meta$numnonspikes<10000)
r1 = marginal_prv(Ycommon[-idxs,], as.matrix(Xfungi[-idxs,-1]), 15, 0.0001, 100) # working
r1 = marginal_prv(Y[-idxs,], as.matrix(Xfungi[-idxs,-1]), 15, 0.0001, 100) # working

sumr = rowSums(otu.table)
idx0 = which(sumr==0)
r1 = marginal_prv(Ycommon[-idx0,], as.matrix(Xfungi[-idx0,-1]), 15, 0.0001, 100) # working
r1 = marginal_prv(Y[-idx0,], as.matrix(Xfungi[-idx0,-1]), 15, 0.0001, 100) # working

idt = unique(c(idxs, idx0))
r1 = marginal_prv(Ycommon[-idt,], as.matrix(Xfungi[-idt,-1]), 15, 0.0001, 100) # working
