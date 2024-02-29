rm(list=ls())

library(tidyverse)

load("allData_clim_trait.Rdata")
source("function_fungi.R")
fungi = matrix(as.numeric(otu.table >0), nrow(otu.table), ncol(otu.table))

###-------### empirical marginal species occurrence probabilities ###-------####

m1 = apply(fungi, 2, mean)
fungi_prevalence = data.frame(Species = seq(1, NCOL(fungi)), Prevalence = m1)
sz_t = 22
pl1 = ggplot(fungi_prevalence) + 
  geom_histogram(aes(x = Prevalence), binwidth = 0.005) + 
  xlab("Mean prevalence") + 
  ylab("Count") + 
  theme_bw()+
  theme(axis.title.x = element_text(size = sz_t), 
        axis.text.x = element_text(size = sz_t),
        axis.title.y = element_text(size = sz_t), 
        axis.text.y = element_text(size = sz_t))
pl1

# data aggregated by phylum
fungi_phylum = data_tree(fungi, taxonomy$phylum)
mp = apply(fungi_phylum,2,mean)
sort(mp, decreasing = T)[1:10]
# data aggregated by class
fungi_class = data_tree(fungi, taxonomy$class)
m2 = apply(fungi_class,2,mean)
sort(m2, decreasing = T)[1:10]


#-----------------------------------------------------------------------------#
# compute separately for each level of the tree the number of nodes with only 
# one children 

phylan = unique(taxonomy$phylum)
classn = unique(taxonomy$class)
phyla_count = rep(0, length(phylan))

for(i in 1:length(phylan)){
  t1 = taxonomy %>% filter(phylum==phylan[i])
  phyla_count[i] = length(unique(t1$class))
}
table(phyla_count>1)

class_count = rep(0, length(classn))
for(i in 1:length(classn)){
  t1 = taxonomy %>% filter(class==classn[i])
  class_count[i] = length(unique(t1$order))
}
table(class_count>1)

ordern = unique(taxonomy$order)
order_count = rep(0, length(ordern))
for(i in 1:length(ordern)){
  t1 = taxonomy %>% filter(order==ordern[i])
  order_count[i] = length(unique(t1$family))
}
table(order_count>1)

famn = unique(taxonomy$family)
fam_count = rep(0, length(famn))
for(i in 1:length(famn)){
  t1 = taxonomy %>% filter(family==famn[i])
 fam_count[i] = length(unique(t1$genus))
}
table(fam_count>1)

genn = unique(taxonomy$genus)
gen_count = rep(0, length(genn))
for(i in 1:length(genn)){
  t1 = taxonomy %>% filter(genus==genn[i])
  gen_count[i] = length(unique(t1$species))
}
table(gen_count>1)


#-----------------------------------------------------------------------------#
# compute the number of nodes which have only one children for all levels


phylan = as.character(phylan)
c1 = matrix(0, length(phylan), 5)
for(i in 1:length(phylan)){
  t1 = taxonomy %>% filter(phylum==phylan[i])
  c1[i,1] = length(unique(t1$class))
  c1[i,2] = length(unique(t1$order))
  c1[i,3] = length(unique(t1$family))
  c1[i,4] = length(unique(t1$genus))
  c1[i,5] = length(unique(t1$species))
}

sc = apply(c1,1,sum)
table(sc==5)
sc2 = apply(c1[,1:4],1,sum)
table(sc2==4)
sc3 = apply(c1[,1:3],1,sum)
table(sc3==3)

