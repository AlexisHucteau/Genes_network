###### METHYLATION ANALYSE

##### Choose the comparison between 5:
#1 -> WT vs MUT
#2 -> WT vs HG
#3 -> HG vs MUT
#4 -> WT vs Octyl
#5 -> HG vs Octyl
number <- 1

##### Choose the logFC threshold for the transcriptome analysis
threshold <- 1.5

#### Choose the pvalue threshold for the transcriptome analysis
pvalue <- 0.1

#### Choose the threshold value for methylation analysis (between 0.5 and 1)
methylation_treshold <- 0.65

setwd("~/Genes_network/")
getwd()

load("~/Genes_network/Methylomes/DATA_Methylomes/Methylomes_data.RData")
system("say finished")

library(stringr)

anno_450 <- read.csv("/home/alexis/Documents/Data_TCGA/Methylomes/HumanMethylation450_15017482_v1-2.csv", as.is = TRUE, skip = 7)
anno_450 <- anno_450[, c("CHR", "MAPINFO", "Name", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "UCSC_RefGene_Group")]

threshold_function <- function(up_value) {
  return(-(up_value**2) + up_value - 0.25)
}

Mut <- as.data.frame(RawnoNABeta[, seq(1:13)])
WT <- as.data.frame(RawnoNABeta[, seq(14, 90)])
Mut$Mean <- rowMeans(Mut)
WT$Mean <- rowMeans(WT)
cpgs_Rawdata <- rownames(RawnoNABeta)
Mean_Beta_Values <- data.frame("WT" = WT$Mean, "Mut" = Mut$Mean, "cpgs" = cpgs_Rawdata)

Mean_Beta_Values <- merge(x = Mean_Beta_Values, y = anno_450, by.x = "cpgs", by.y = "Name")

Promoter_Beta_values <- Mean_Beta_Values # [select, ]

true_methylation_treshold <- threshold_function(methylation_treshold)

nb_upmethylation <- sapply((Mean_Beta_Values[, 3] - 0.5) * (Mean_Beta_Values[, 2] - 0.5), function(test_modification) {
  ifelse(test_modification < true_methylation_treshold, TRUE, FALSE)
})
Mean_Beta_Values$up_down_methylation <- nb_upmethylation

Beta_values_change <- Mean_Beta_Values[which(Mean_Beta_Values$up_down_methylation), ]

Beta_values_change <- Beta_values_change[, c("cpgs", "WT", "Mut", "CHR", "MAPINFO", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "UCSC_RefGene_Group")]

######TRANSCRIPTOMIC ANALYSES

rawDataDir <- "Transcriptomes/DATA_Transcriptomes/"

# library(affy)
library(sva)
library(limma)
library(dplyr)
library(tidyr)
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)
# BiocManager::install("pd.hugene.2.0.st")
library(hugene20sttranscriptcluster.db)
library(oligo)

SDRF <- read.csv(file = "Transcriptomes/Samplesheet_Transcriptomes/Samplesheet.csv")
celFiles <- list.celfiles(rawDataDir, full.names = TRUE)
# data <- oligo::read.celfiles(filenames = file.path(rawDataDir, SDRF$Array.Data.File), verbose = FALSE, phenoData = SDRF, checkType = FALSE)

data <- read.celfiles(celFiles)
eset <- oligo::rma(data) # si non normalisé
my_frame <- as.data.frame(exprs(eset))

Annot <- data.frame(ACCNUM = sapply(contents(hugene20sttranscriptclusterACCNUM), paste, collapse = ", "), SYMBOL = sapply(contents(hugene20sttranscriptclusterSYMBOL), paste, collapse = ", "), DESC = sapply(contents(hugene20sttranscriptclusterGENENAME), paste, collapse = ", "))
all <- merge(Annot, my_frame, by.y = 0, by.x = 0, all = T)

TS <- paste(SDRF$Characteristics.genotype., SDRF$Characteristics.treatment., sep = ".")
TS <- factor(TS, levels = c("WT.HG", "WT.octyl", "WT.None", "Mut.None"))
design <- model.matrix(~ 0 + TS)
colnames(design) <- levels(TS)
fit <- lmFit(eset, design)

cont.matrix <- makeContrasts(
  WT.None.vsMut.None = Mut.None - WT.None,
  WT.None.vsWT.HG = WT.HG - WT.None,
  WT.HGvsMut.None = Mut.None - WT.HG,
  WT.NonevsWT.Octyl = WT.octyl - WT.None,
  WT.HGvsWT.octyl = WT.octyl - WT.HG,
  levels = design
)

library(data.table)

# Select the contrast
colnames(cont.matrix)

colnames(cont.matrix)[number]
fit.contrast <- contrasts.fit(fit, cont.matrix)
efit.contrast <- eBayes(fit.contrast)

Genes_values <- list()
Genes_values[["UP"]] <- list()
Genes_values[["DOWN"]] <- list()

# get all the list   #coef parameter = contrast
genes <- Annot$SYMBOL
select <- sapply(genes, function(gene) {
  !is.null(gene)
})
genes <- genes[select]

all_val <- topTable(efit.contrast, coef = number, adjust.method = "BY", n = Inf, genelist = genes)
select_na <- sapply(all_val$ID, function(gene) {
  !(str_detect(gene, "NA"))
})
all_val <- all_val[select_na, ]

# "threshold" for the logFC threshold
# "pvalue" for the pvalue threshold
signif <- all_val[which(all_val$logFC <= -threshold | all_val$logFC >= threshold), ]
signif <- signif[which(signif$P.Value <= pvalue), ]


UP <- signif[which(signif$logFC > 0), ] # UP
UP <- UP[order(UP$logFC, decreasing = T), ]
nrow(UP)
genes_up <- UP$ID
transit <- tstrsplit(as.vector(genes_up), "_///_")
genes_up <- c(transit[[1]]) # ,transit[[4]])
genes_up <- genes_up[which(genes_up != "NA")]
Genes_values[["UP"]][[colnames(cont.matrix)[number]]] <- array(genes_up)

DOWN <- signif[which(signif$logFC < 0), ] # DOWN
DOWN <- DOWN[order(DOWN$logFC, decreasing = T), ]
nrow(DOWN)
genes_down <- DOWN$ID
transit <- tstrsplit(as.vector(genes_down), "_///_")
genes_down <- c(transit[[1]]) # ,transit[[4]])
genes_down <- genes_down[which(genes_down != "NA")]
Genes_values[["DOWN"]][[colnames(cont.matrix)[number]]] <- array(genes_down)

##################GENE EXPENSION


library(igraph)
library(readr)
library(dplyr)

biogrid_data <- read_tsv("Genes_Expansion/Networks/BIOGRID-ORGANISM-Homo_sapiens-3.5.185.tab3.txt") %>%
  filter(`Organism Name Interactor A` == "Homo sapiens" & `Organism Name Interactor B` == "Homo sapiens")
reactome_FI <- read_tsv("Genes_Expansion/Networks/FIsInGene_020720_with_annotations.tsv")

biogrid_network <- graph_from_data_frame(biogrid_data, directed = F)
reactomeFI_network <- graph_from_data_frame(reactome_FI, directed = F)

explore_genes_network <- function(genes_vector) {
  Gene_list <- list()
  selected_nodes <- V(biogrid_network)[V(biogrid_network)$name %in% genes_vector]
  selected_edges <- E(biogrid_network)[from(selected_nodes) | to(selected_nodes)]
  subnet <- subgraph.edges(biogrid_network, selected_edges)
  nodes <- as_tibble(get.data.frame(subnet, what = "vertices"))
  nodes <- data.frame("gene" = nodes, "DE" = nodes$name %in% genes_vector)
  edges <- as_tibble(get.data.frame(subnet, what = "edges"))
  Gene_list[["biogrid"]] <- nodes
  
  selected_nodes <- V(reactomeFI_network)[V(reactomeFI_network)$name %in% genes_vector]
  # First neighborhood is all nodes connected to our selected nodes
  selected_edges <- E(reactomeFI_network)[from(selected_nodes) | to(selected_nodes)]
  subnet <- subgraph.edges(reactomeFI_network, selected_edges)
  nodes <- as_tibble(get.data.frame(subnet, what = "vertices"))
  nodes <- data.frame("gene" = nodes, "DE" = nodes$name %in% genes_vector)
  edges <- as_tibble(get.data.frame(subnet, what = "edges"))
  Gene_list[["reactome"]] <- nodes
  return(Gene_list)
}

Genes_expanded <- list()
Genes_expanded[["UP"]] <- explore_genes_network(Genes_values[["UP"]][[colnames(cont.matrix)[number]]])
Genes_expanded[["DOWN"]] <- explore_genes_network(Genes_values[["DOWN"]][[colnames(cont.matrix)[number]]])


##################Chromatin CPG


library("igraph")
library("GenomicRanges")
library("tidyverse")

load("Chromatine_network/DATA/pchic.RData")
pchic <- pchic[, c(1:10)]

Cpgs_Coordinates <- Beta_values_change

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
  seqnames = Cpgs_Coordinates$CHR,
  ranges = IRanges(start = Cpgs_Coordinates$MAPINFO, end = Cpgs_Coordinates$MAPINFO + 1),
  methylation = ifelse(Cpgs_Coordinates$Mut - Cpgs_Coordinates$WT > 0, "UP", "DOWN"),
  genes = Cpgs_Coordinates$UCSC_RefGene_Name
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

#################FINAL NETWORK

########UP DE#########

Genes_DE <- Genes_expanded[["UP"]][["biogrid"]]
Genes_network <- as.array(Genes_DE$name)
Genes_DE_genes <- Genes_DE[which(Genes_DE$DE),1]

select_from_network <- sapply(str_split(Node_features$Genes, ";"), function(genes) {
  length(intersect(unlist(genes), Genes_network)) > 0
})

select_DE_gene <- sapply(str_split(Node_features$Genes, ";"), function(genes) {
  length(intersect(unlist(genes), Genes_DE_genes)) > 0
})

selected <- data.frame(cbind(select_from_network, select_DE_gene))
selected$selection <- ifelse(selected$select_from_network, ifelse(selected$select_DE_gene, "UP", "Unetwork"), "NinN")

########DOWN DE#########

Genes_DE <- Genes_expanded[["DOWN"]][["biogrid"]]
Genes_network <- as.array(Genes_DE$name)
Genes_DE_genes <- Genes_DE[which(Genes_DE$DE),1]

select_from_network <- sapply(str_split(Node_features$Genes, ";"), function(genes) {
  length(intersect(unlist(genes), Genes_network)) > 0
})

select_DE_gene <- sapply(str_split(Node_features$Genes, ";"), function(genes) {
  length(intersect(unlist(genes), Genes_DE_genes)) > 0
})

selected_DOWN <- data.frame(cbind(select_from_network, select_DE_gene))
selected_DOWN$selection <- ifelse(selected_DOWN$select_from_network, ifelse(selected_DOWN$select_DE_gene, "DOWN", "Dnetwork"), "NinN")

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

write.csv(CpGs_all_network[, c(1, 5)], file = paste0("Final_network/Results/CpGs_chromatine_Network_", colnames(cont.matrix)[number], "_pvalue_", pvalue, "_logFC_", threshold, "_methreshold_", methylation_treshold, ".csv"), row.names = FALSE)
write.csv(Node_features, file = paste0("Final_network/Results/Nodes_features_", colnames(cont.matrix)[number], "_pvalue_", pvalue, "_logFC_", threshold, "_methreshold_", methylation_treshold,  ".csv"), row.names = FALSE)

Number_of_total_CpG <- length(RawnoNABeta[,1])
Number_of_CpG_UP_in_MUT <- length(Mean_Beta_Values[which(Mean_Beta_Values$Mut > 0.5),1])
Number_of_CpG_DOWN_in_MUT <- length(Mean_Beta_Values[which(Mean_Beta_Values$Mut < 0.5),1])
Number_of_CpG_UP_in_WT <- length(Mean_Beta_Values[which(Mean_Beta_Values$WT > 0.5),1])
Number_of_CpG_DOWN_in_WT <- length(Mean_Beta_Values[which(Mean_Beta_Values$WT < 0.5),1])

Number_of_Gene_DE <- length(Genes_values[["UP"]][["WT.None.vsMut.None"]]) + length(Genes_values[["DOWN"]][["WT.None.vsMut.None"]])
Number_of_Gene_expended <- length(Genes_expanded[["UP"]][["biogrid"]][["name"]]) + length(Genes_expanded[["DOWN"]][["biogrid"]][["name"]])

Number_of_extended_genes_in_network <- length(Node_features[!str_detect(Node_features$final_feature, "Not") & Node_features$CpGs,1])
Number_of_extended_genes_in_network <- length(Node_features[!str_detect(Node_features$final_feature, "Not"),1])


tmp <- anno_450[which(anno_450$Name %in% rownames(RawnoNABeta)), ]
Number_of_CpG_in_promoter <- length(tmp[which(str_detect(tmp$UCSC_RefGene_Group, "TSS")),1])
Number_of_CpG_DM <- length(Beta_values_change$cpgs)
Number_of_CpG_UP <- length(Beta_values_change[which(Beta_values_change$Mut > 0.5),1])
Number_of_CpG_DOWN <- length(Beta_values_change[which(Beta_values_change$Mut < 0.5),1])
Number_of_CpG_DM_in_Promoter <- length(Beta_values_change[which(str_detect(Beta_values_change$UCSC_RefGene_Group, "TSS")),1])
CpG_in_gene_neighbor <- length(Node_features[!str_detect(Node_features$final_feature, "Not"),1])
Fragment_with_DM <- length(Node_features[Node_features$CpGs,1])
CpG_in_gene_DE <- length(Node_features[!str_detect(Node_features$final_feature, "_"),1])
CpG_out_gene_DE <- CpG_in_gene_neighbor - CpG_in_gene_DE



length(Node_features$ID)
length(Node_features[Node_features$CpGs,1])
length(Node_features[Node_features$Promoter & Node_features$CpGs,1])
length(Node_features[!str_detect(Node_features$final_feature, "Not"), 1])
length(Node_features[!str_detect(Node_features$final_feature, "Not") & Node_features$CpGs, 1])
