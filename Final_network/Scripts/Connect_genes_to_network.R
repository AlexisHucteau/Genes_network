library("tidyverse")

setwd("Genes_network/Final_network/")
getwd()

Nodes_features <- read.csv(file = "../Chromatine_network/Results/Node_features.csv")

########UP DE#########

Genes_DE <- read_tsv(file = "../Genes_Expansion/Networks_of_new_genes/WT_vs_MUT_UP/biogrid/WT_vs_MUT_UP_nodes.tsv")
Genes_network <- as.array(Genes_DE$gene)
Genes_DE_genes <- Genes_DE[which(Genes_DE$DE),1]$gene

select_from_network <- sapply(str_split(Nodes_features$Genes, ";"), function(genes) {
  length(intersect(unlist(genes), Genes_network)) > 0
})

select_DE_gene <- sapply(str_split(Nodes_features$Genes, ";"), function(genes) {
  length(intersect(unlist(genes), Genes_DE_genes)) > 0
})

selected <- data.frame(cbind(select_from_network, select_DE_gene))
selected$selection <- ifelse(selected$select_from_network, ifelse(selected$select_DE_gene, "UP", "Unetwork"), "NinN")
selected[which(selected$select_DE_gene),]

########DOWN DE#########

Genes_DE <- read_tsv(file = "../Genes_Expansion/Networks_of_new_genes/WT_vs_MUT_DOWN/biogrid/WT_vs_MUT_DOWN_nodes.tsv")
Genes_network <- as.array(Genes_DE$gene)
Genes_DE_genes <- Genes_DE[which(Genes_DE$DE),1]$gene

select_from_network <- sapply(str_split(Nodes_features$Genes, ";"), function(genes) {
  length(intersect(unlist(genes), Genes_network)) > 0
})

select_DE_gene <- sapply(str_split(Nodes_features$Genes, ";"), function(genes) {
  length(intersect(unlist(genes), Genes_DE_genes)) > 0
})

selected_DOWN <- data.frame(cbind(select_from_network, select_DE_gene))
selected_DOWN$selection <- ifelse(selected_DOWN$select_from_network, ifelse(selected_DOWN$select_DE_gene, "DOWN", "Dnetwork"), "NinN")
selected_DOWN[which(selected_DOWN$select_DE_gene),]

#########MERGING#############

selection_Genes_UP_DOWN <- data.frame("UP" = selected$selection, "DOWN" = selected_DOWN$selection)
selection_Genes_UP_DOWN$merge <- paste0(selection_Genes_UP_DOWN$UP, selection_Genes_UP_DOWN$DOWN)

final_feature <- sapply(selection_Genes_UP_DOWN$merge, function(merging) {
  if (str_detect(merging, "UP")) {
    return("UP")
  }else if (str_detect(merging, "DOWN")) {
    return("DOWN")
  }else if (str_detect(merging, "Dnetwork")) {
    return("Down_network")
  }else if (str_detect(merging, "Unetwork")) {
    return("UP_network")
  }else {
    return("Not_in_network")
  }
})

Nodes_features$final_feature <- final_feature

write.csv(Nodes_features, file = "Results/Nodes_features.csv", row.names = FALSE)
