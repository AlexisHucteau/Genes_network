library("tidyverse")

setwd("Genes_network/Final_network/")
getwd()

Nodes_features <- read.csv(file = "../Chromatine_network/Results/Node_features.csv")
Genes_DE <- read_tsv(file = "../Genes_Expansion/Networks_of_new_genes/WT_vs_MUT_UP/biogrid/WT_vs_MUT_UP_nodes.tsv")
Genes_DE <- as.array(Genes_DE$gene)

select <- sapply(str_split(Nodes_features$Genes, ";"), function(genes) {
  length(intersect(unlist(genes), Genes_DE)) > 0
})

Nodes_features$DE <- select
write.csv(Nodes_features, file = "Results/Nodes_features.csv")
