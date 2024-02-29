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