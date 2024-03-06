require(dplyr)

data_tree = function(Y, label){
  # return the dataset grouped by the taxonomic level choose in 'label'
  data_all = cbind.data.frame(t(Y), "taxa" = label)
  data_taxa = data_all %>% 
    group_by(taxa) %>%
    summarise(across(where(is.numeric), max))
  data_taxa = as.data.frame(t(data_taxa))
  colnames(data_taxa) = data_taxa[1,]
  data_taxa = data_taxa[-1,]
  data_taxa= as.data.frame(sapply(data_taxa, as.numeric))
  data_taxa
}

# same as above, but drop only child
dad_tree = function(Y, taxonomy, l){
  taxonomy = taxonomy[,(l-1):l]
  colnames(taxonomy) = c("vl1","vl")
  data_all = cbind.data.frame(t(Y), "taxa" = taxonomy[,2])
  data_taxa = data_all %>% 
    group_by(taxa) %>%
    summarise(across(where(is.numeric), max))
  ddrop = taxonomy %>% 
    group_by(vl1) %>%
    summarise(ct = vl[which(length(unique(vl))==1)])
  data_taxa = data_taxa[!(data_taxa$taxa %in% ddrop$ct),]
  data_taxa = as.data.frame(t(data_taxa))
  colnames(data_taxa) = data_taxa[1,]
  data_taxa = data_taxa[-1,]
  data_taxa= as.data.frame(sapply(data_taxa, as.numeric))
  data_taxa
}

