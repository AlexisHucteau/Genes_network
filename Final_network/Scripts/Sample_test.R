library("igraph")
library("GenomicRanges")
library("tidyverse")

setwd("Genes_network/Chromatine_network/")
getwd()

anno_450 <- read.csv("/home/alexis/Documents/Data_TCGA/Methylomes/HumanMethylation450_15017482_v1-2.csv", as.is = TRUE, skip = 7)

anno_450 <- anno_450[, c("CHR", "MAPINFO", "Name", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "UCSC_RefGene_Group")]

Sample_test <- sample(anno_450$Name, size = 838)
Sample_test <- anno_450[anno_450$Name %in% Sample_test,]

# PCHiC experiment
load("./DATA/pchic.RData")
pchic <- pchic[, c(1:10)]

List_Promoter <- paste(pchic$baitChr, pchic$baitStart, sep = "_")

# Put all bait and all OE regions into a BED file format and then GRanges object

colnames(pchic)[c(1:5, 6:10)] <- rep(c("chr", "start", "end", "ID", "Name"), 2)
PCHiC_bed <- unique(rbind(pchic[, c(1:3, 5)], pchic[, c(6:8, 10)]))
PCHiC_GRange <- GRanges(
  seqnames = PCHiC_bed$chr,
  IRanges(start = PCHiC_bed$start, end = PCHiC_bed$end),
  Gene_Pchic = PCHiC_bed$Name
)
PCHiC_GRange$ID <- paste(PCHiC_bed$chr, PCHiC_bed$start, sep = "_")
PCHiC_GRange

CpGs_Ranges <- GRanges(
  seqnames = Sample_test$CHR,
  ranges = IRanges(start = Sample_test$MAPINFO, end = Sample_test$MAPINFO + 1),
  genes = Sample_test$UCSC_RefGene_Name
)
CpGs_Ranges

# Overlap peaks with PHiC fragments

overlaps <- findOverlaps(PCHiC_GRange, CpGs_Ranges)

## matching overlapping CpGs IDs (start) with ChIP-seq features

match_hit <- data.frame(mcols(PCHiC_GRange[queryHits(overlaps)]), as.data.frame(mcols(CpGs_Ranges)[subjectHits(overlaps), ]), stringsAsFactors = T)
match_hit$CpGs <- TRUE

colnames(pchic) <- c("chr_bait", "start_bait", "end_bait", "ID_bait", "Name_bait", "chr_oe", "start_oe", "end_oe", "ID_oe", "Name_oe")
pchic$IDbait <- paste(pchic$chr_bait, pchic$start_bait, sep = "_")
pchic$IDoe <- paste(pchic$chr_oe, pchic$start_oe, sep = "_")

CpGs_all_network <- unique(rbind(pchic[which(pchic$IDbait %in% match_hit$ID), ], pchic[which(pchic$IDoe %in% match_hit$ID), ]))
CpGs_all_network <- CpGs_all_network[, c(11, 5, 12, 10)]

CpGs_Promoter_Promoter <- CpGs_all_network[which(CpGs_all_network$IDoe %in% List_Promoter), ]
CpGs_Promoter_Promoter$OePromoter <- TRUE
CpGs_Promoter_NonPromoter <- CpGs_all_network[-which(CpGs_all_network$IDoe %in% List_Promoter), ]
CpGs_Promoter_NonPromoter$OePromoter <- FALSE
CpGs_all_network <- rbind(CpGs_Promoter_Promoter, CpGs_Promoter_NonPromoter)
CpGs_all_network$CpGs_bait <- ifelse(CpGs_all_network$IDbait %in% match_hit$ID, TRUE, FALSE)
CpGs_all_network$CpGs_oe <- ifelse(CpGs_all_network$IDoe %in% match_hit$ID, TRUE, FALSE)
CpGs_all_network$PromoterBait <- TRUE
CpGs_all_network <- CpGs_all_network[, c(1, 2, 6, 8, 3, 4, 7, 5)]
Node_features <- data.frame("ID" = paste0("\"", CpGs_all_network$IDbait, "\""), "Genes" = CpGs_all_network$Name_bait, "Promoter" = CpGs_all_network$PromoterBait, "CpGs" = CpGs_all_network$CpGs_bait)
Node_features <- rbind(Node_features, data.frame("ID" = paste0("\"", CpGs_all_network$IDoe, "\""), "Genes" = CpGs_all_network$Name_oe, "Promoter" = CpGs_all_network$OePromoter, "CpGs" = CpGs_all_network$CpGs_oe))
Node_features <- unique(Node_features)

# For any reason, the bait ID : 19_37726761 and 7_1261811 are duplicated.

CpGs_all_network <- CpGs_all_network[-which(str_detect(CpGs_all_network$IDbait, "19_37726761")), ]
CpGs_all_network <- CpGs_all_network[-which(str_detect(CpGs_all_network$IDbait, "7_1261811")), ]
Node_features <- Node_features[-which(str_detect(Node_features$ID, "\"19_37726761\"")), ]
Node_features <- Node_features[-which(str_detect(Node_features$ID, "\"7_1261811\"")), ]

Node_features$Cytoscape_feature <- "string"

for (i in seq(length(Node_features$CpGs))) {
  if (Node_features[i, 3] & Node_features[i, 4]) {
    Cytoscape_feature <- "Promoter_methylated"
  } else if (Node_features[i, 3] & !Node_features[i, 4]) {
    Cytoscape_feature <- "Promoter"
  } else if (!Node_features[i, 3] & Node_features[i, 4]) {
    Cytoscape_feature <- "Methylated"
  } else {
    Cytoscape_feature <- "Simple_fragment"
  }
  Node_features[i, 5] <- Cytoscape_feature
}










Genes_DE <- read_tsv(file = "../Genes_Expansion/Networks_of_new_genes/WT_vs_MUT_UP/biogrid/WT_vs_MUT_UP_nodes.tsv")
Genes_network <- as.array(Genes_DE$gene)
Genes_DE_genes <- Genes_DE[which(Genes_DE$DE),1]$gene

select_from_network <- sapply(str_split(Node_features$Genes, ";"), function(genes) {
  length(intersect(unlist(genes), Genes_network)) > 0
})

select_DE_gene <- sapply(str_split(Node_features$Genes, ";"), function(genes) {
  length(intersect(unlist(genes), Genes_DE_genes)) > 0
})

selected <- data.frame(cbind(select_from_network, select_DE_gene))
selected$selection <- ifelse(selected$select_from_network, ifelse(selected$select_DE_gene, "UP", "Unetwork"), "NinN")
selected[which(selected$select_DE_gene),]

########DOWN DE#########

Genes_DE <- read_tsv(file = "../Genes_Expansion/Networks_of_new_genes/WT_vs_MUT_DOWN/biogrid/WT_vs_MUT_DOWN_nodes.tsv")
Genes_network <- as.array(Genes_DE$gene)
Genes_DE_genes <- Genes_DE[which(Genes_DE$DE),1]$gene

select_from_network <- sapply(str_split(Node_features$Genes, ";"), function(genes) {
  length(intersect(unlist(genes), Genes_network)) > 0
})

select_DE_gene <- sapply(str_split(Node_features$Genes, ";"), function(genes) {
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

Node_features$final_feature <- final_feature


write.csv(Node_features, file = "Results/Node_features_sample_test.csv", row.names = FALSE)
write.csv(CpGs_all_network[, c(1, 5)], file = "Results/CpGs_chromatine_Network_sample_test.csv", row.names = FALSE)

length(Node_features$ID)
length(Node_features[Node_features$CpGs,1])
length(Node_features[Node_features$Promoter & Node_features$CpGs,1])
length(Node_features[!str_detect(Node_features$final_feature, "Not"), 1])
length(Node_features[!str_detect(Node_features$final_feature, "Not") & Node_features$CpGs, 1])
