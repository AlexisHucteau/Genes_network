library(ChAMP)
library(stringr)
library(sva)
library(limma)
library(dplyr)
library(tidyr)
library(matrixStats)
library(genefilter)
library(openxlsx)
# BiocManager::install("pd.hugene.2.0.st")
library(hugene20sttranscriptcluster.db)
library(oligo)
library(data.table)
library(igraph)
library(readr)
library("GenomicRanges")
library("tidyverse")

################PRELOADING################

load("~/Genes_network/Methylomes/DATA_Methylomes/Methylomes_data.RData")

anno_450 <- read.csv("/home/alexis/Documents/Data_TCGA/Methylomes/HumanMethylation450_15017482_v1-2.csv", as.is = TRUE, skip = 7)
anno_450 <- anno_450[, c("CHR", "MAPINFO", "Name", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "UCSC_RefGene_Group")]

rawDataDir <- "~/Genes_network/Transcriptomes/DATA_Transcriptomes/"

SDRF <- read.csv(file = "~/Genes_network/Transcriptomes/Samplesheet_Transcriptomes/Samplesheet.csv")
celFiles <- list.celfiles(rawDataDir, full.names = TRUE)
data <- read.celfiles(celFiles)
eset <- oligo::rma(data) # si non normalisé
my_frame <- as.data.frame(exprs(eset))
Annot <- data.frame(ACCNUM = sapply(contents(hugene20sttranscriptclusterACCNUM), paste, collapse = ", "), SYMBOL = sapply(contents(hugene20sttranscriptclusterSYMBOL), paste, collapse = ", "), DESC = sapply(contents(hugene20sttranscriptclusterGENENAME), paste, collapse = ", "))
all <- merge(Annot, my_frame, by.y = 0, by.x = 0, all = T)

biogrid_data <- read_tsv("~/Genes_network/Genes_Expansion/Networks/BIOGRID-ORGANISM-Homo_sapiens-3.5.185.tab3.txt") %>%
  filter(`Organism Name Interactor A` == "Homo sapiens" & `Organism Name Interactor B` == "Homo sapiens")
reactome_FI <- read_tsv("~/Genes_network/Genes_Expansion/Networks/FIsInGene_020720_with_annotations.tsv")

biogrid_network <- graph_from_data_frame(biogrid_data, directed = F)
reactomeFI_network <- graph_from_data_frame(reactome_FI, directed = F)

load("~/Genes_network/Chromatine_network/DATA/pchic.RData")
pchic <- pchic[, c(1:10)]

List_Promoter <- paste(pchic$baitChr, pchic$baitStart, sep = "_")
List_Promoter <- unique(List_Promoter)

colnames(pchic)[c(1:5, 6:10)] <- rep(c("chr", "start", "end", "ID", "Name"), 2)
PCHiC_bed <- unique(rbind(pchic[, c(1:3, 5)], pchic[, c(6:8, 10)]))
PCHiC_GRange <- GRanges(
  seqnames = PCHiC_bed$chr,
  IRanges(start = PCHiC_bed$start, end = PCHiC_bed$end),
  Gene_Pchic = PCHiC_bed$Name
)
PCHiC_GRange$ID <- paste(PCHiC_bed$chr, PCHiC_bed$start, sep = "_")

colnames(pchic) <- c("chr_bait", "start_bait", "end_bait", "ID_bait", "Name_bait", "chr_oe", "start_oe", "end_oe", "ID_oe", "Name_oe")
pchic$IDbait <- paste(pchic$chr_bait, pchic$start_bait, sep = "_")
pchic$IDoe <- paste(pchic$chr_oe, pchic$start_oe, sep = "_")









############Functions##################


#Methylations functions
#Return DATA analysed
Analyze_methylation <- function(DATA, Samplesheet, annotations) {
  TS <- ifelse(str_detect(Samplesheet$Sample_name, "MUT"), "Mut", "WT")
  TS <- factor(TS, levels = c("Mut", "WT"))
  design <- model.matrix(~ 0 + TS)
  colnames(design) <- levels(TS)
  fit <- lmFit(DATA, design)
  cont.matrix <- makeContrasts(
    WTvsMut = Mut - WT,
    levels = design
  )
  fit.contrast <- contrasts.fit(fit, cont.matrix)
  efit.contrast <- eBayes(fit.contrast)
  Methylation_contrast <- data.frame(logFC = efit.contrast[["coefficients"]][, 1], pvalue = efit.contrast[["p.value"]][, 1])
  Methylation_contrast$cpg <- rownames(Methylation_contrast)
  Methylation_contrast <- merge(x = Methylation_contrast, y = annotations, by.x = "cpg", by.y = "Name")
  Mut_column <- str_detect(TS, "Mut")
  WT_column <- str_detect(TS, "WT")
  Mut <- as.data.frame(RawnoNABeta[, Mut_column])
  WT <- as.data.frame(RawnoNABeta[, WT_column])
  Mut$Mean <- rowMeans(Mut)
  WT$Mean <- rowMeans(WT)
  cpgs_Rawdata <- rownames(RawnoNABeta)
  Mean_Beta_Values <- data.frame("WT" = WT$Mean, "Mut" = Mut$Mean, "cpgs" = cpgs_Rawdata)
  Methylation_contrast <- merge(x = Methylation_contrast, y = Mean_Beta_Values, by.x = "cpg", by.y = "cpgs")
  return(Methylation_contrast)
}

#Return DATA filtered
Filter_Differential_Methylation <- function(DATA, pvalue, logFC_treshold) {
  High_contrast <- DATA[DATA$pvalue < pvalue & DATA$logFC**2 > logFC_treshold**2, c(1, 9, 10, 2, 3, 4, 5, 6, 7, 8)]
  return(High_contrast)
}

#Transcription function

#Return a list of genes Differentially expressed
Analyse_transcription <- function(number_of_the_comparison, 
                                  DATA, 
                                  Samplesheet, 
                                  annotation, 
                                  logFC_threshold, 
                                  pvalue) {
  TS <- paste(Samplesheet$Characteristics.genotype., Samplesheet$Characteristics.treatment., sep = ".")
  TS <- factor(TS, levels = c("WT.HG", "WT.octyl", "WT.None", "Mut.None"))
  design <- model.matrix(~ 0 + TS)
  colnames(design) <- levels(TS)
  fit <- lmFit(DATA, design)
  cont.matrix <- makeContrasts(
    WT.None.vsMut.None = Mut.None - WT.None,
    WT.None.vsWT.HG = WT.HG - WT.None,
    WT.HGvsMut.None = Mut.None - WT.HG,
    WT.NonevsWT.Octyl = WT.octyl - WT.None,
    WT.HGvsWT.octyl = WT.octyl - WT.HG,
    levels = design
  )
  fit.contrast <- contrasts.fit(fit, cont.matrix)
  efit.contrast <- eBayes(fit.contrast)
  Genes_values <- list()
  genes <- annotation$SYMBOL
  select <- sapply(genes, function(gene) {
    !is.null(gene)
  })
  genes <- genes[select]
  all_val <- topTable(efit.contrast, coef = number_of_the_comparison, adjust.method = "BY", n = Inf, genelist = genes)
  select_na <- sapply(all_val$ID, function(gene) {
    !(str_detect(gene, "NA"))
  })
  all_val <- all_val[select_na, ]
  signif <- all_val[which(all_val$logFC <= -logFC_threshold | all_val$logFC >= logFC_threshold), ]
  signif <- signif[which(signif$P.Value <= pvalue), ]
  UP <- signif[which(signif$logFC > 0), ] # UP
  UP <- UP[order(UP$logFC, decreasing = T), ]
  genes_up <- UP$ID
  transit <- tstrsplit(as.vector(genes_up), "_///_")
  genes_up <- c(transit[[1]]) # ,transit[[4]])
  genes_up <- genes_up[which(genes_up != "NA")]
  Genes_values[["UP"]] <- array(genes_up)
  DOWN <- signif[which(signif$logFC < 0), ] # DOWN
  DOWN <- DOWN[order(DOWN$logFC, decreasing = T), ]
  genes_down <- DOWN$ID
  transit <- tstrsplit(as.vector(genes_down), "_///_")
  genes_down <- c(transit[[1]]) # ,transit[[4]])
  genes_down <- genes_down[which(genes_down != "NA")]
  Genes_values[["DOWN"]] <- array(genes_down)
  return(Genes_values)
}

#Gene_expansion_function

#Return a vector of genes from a network
explore_genes_network <- function(genes_vector, network) {
  selected_nodes <- V(network)[V(network)$name %in% genes_vector]
  selected_edges <- E(network)[from(selected_nodes) | to(selected_nodes)]
  subnet <- subgraph.edges(network, selected_edges)
  nodes <- as_tibble(get.data.frame(subnet, what = "vertices"))
  nodes <- data.frame("gene" = nodes, "DE" = nodes$name %in% genes_vector)
  edges <- as_tibble(get.data.frame(subnet, what = "edges"))
  return(nodes)
}

#Connect_CpG_to_chromatine_network

#Return 
find_overlap_between_data_and_chromatine_network <- function(DATA, network) {
  DATA <- DATA[str_detect(DATA$cpg, "cg"),]
  CpGs_Ranges <- GRanges(
    seqnames = DATA$CHR,
    ranges = IRanges(start = DATA$MAPINFO, end = DATA$MAPINFO + 1),
    genes = DATA$UCSC_RefGene_Name,
    position = DATA$UCSC_RefGene_Group
  )
  overlaps <- findOverlaps(network, CpGs_Ranges)
  match_hit <- data.frame(mcols(network[queryHits(overlaps)]), as.data.frame(mcols(CpGs_Ranges)[subjectHits(overlaps), ]), stringsAsFactors = T)
  match_hit$CpGs <- TRUE
  return(match_hit)
}

#Return a network
Create_network <- function(DATA, network, list_of_promoter) {
  CpGs_all_network <- unique(rbind(network[which(network$IDbait %in% DATA$ID), ], network[which(network$IDoe %in% DATA$ID), ]))
  CpGs_all_network <- CpGs_all_network[, c(11, 5, 12, 10)]
  CpGs_Promoter_Promoter <- CpGs_all_network[which(CpGs_all_network$IDoe %in% list_of_promoter), ]
  CpGs_Promoter_Promoter$OePromoter <- TRUE
  CpGs_Promoter_NonPromoter <- CpGs_all_network[-which(CpGs_all_network$IDoe %in% list_of_promoter), ]
  CpGs_Promoter_NonPromoter$OePromoter <- FALSE
  CpGs_all_network <- rbind(CpGs_Promoter_Promoter, CpGs_Promoter_NonPromoter)
  CpGs_all_network$CpGs_bait <- ifelse(CpGs_all_network$IDbait %in% DATA$ID, TRUE, FALSE)
  CpGs_all_network$CpGs_oe <- ifelse(CpGs_all_network$IDoe %in% DATA$ID, TRUE, FALSE)
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
  res <- list()
  res[["network"]] <- CpGs_all_network
  res[["features"]] <- Node_features
  return(res)
  
}

#ADD transcription feature

#create feature from vector of genes
Transcription_feature <- function(DATA, network_features) {
  Genes_network <- as.array(DATA$name)
  Genes_DE_genes <- DATA[which(DATA$DE),1]
  
  select_from_network <- sapply(str_split(network_features$Genes, ";"), function(genes) {
    length(intersect(unlist(genes), Genes_network)) > 0
  })
  
  select_DE_gene <- sapply(str_split(network_features$Genes, ";"), function(genes) {
    length(intersect(unlist(genes), Genes_DE_genes)) > 0
  })
  
  selected <- data.frame(cbind(select_from_network, select_DE_gene))
  selected$selection <- ifelse(selected$select_from_network, ifelse(selected$select_DE_gene, "UP", "Unetwork"), "NinN")
  return(selected)
}

#Add features from 2 sens expression genes
Transcription_features <- function(gene_list, network_features) {
  selected_UP <- Transcription_feature(gene_list[["UP"]], network_features)
  selected_DOWN <- Transcription_feature(gene_list[["DOWN"]], network_features)
  selection_Genes_UP_DOWN <- data.frame("UP" = selected_UP$selection, "DOWN" = selected_DOWN$selection)
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
  
  network_features$final_feature <- final_feature
  return(network_features)
}


Filter_network <- function(network) {
  network <- network[network$CpGs_bait | network$CpGs_oe | network$OePromoter,]
  network <- network[-which(network$Name_oe == "."),]
  return(network)
}



From_data_to_network <- function(DATA_MET, 
                                 Samplesheet_MET, 
                                 annotations_MET, 
                                 pvalue_met = 0.01, 
                                 logFC_treshold = 0.1,
                                 number_of_the_comparison = number, 
                                 DATA_RNAseq, 
                                 Samplesheet_RNAseq, 
                                 annotation_RNAseq, 
                                 logFC_threshold_RNAseq, 
                                 pvalue_RNAseq, 
                                 chromatine_BED_file = PCHiC_GRange, 
                                 chromatine_network = pchic, 
                                 List_of_promoter = List_Promoter, 
                                 gene_gene_network = reactomeFI_network) {
  Analysed_DATA_Met <- Analyze_methylation(DATA_MET, Samplesheet_MET, annotations_MET)
  print("Analyze_methylation: DONE")
  Filtered_DATA_Met <- Filter_Differential_Methylation(Analysed_DATA_Met, pvalue_met, logFC_treshold)
  print("Filter_Differential_Methylation: DONE")
  Genes_Differentially_expressed <- Analyse_transcription(number_of_the_comparison = number_of_the_comparison, 
                                                          DATA = DATA_RNAseq, 
                                                          Samplesheet = Samplesheet_RNAseq, 
                                                          annotation = annotation_RNAseq, 
                                                          logFC_threshold = logFC_threshold_RNAseq, 
                                                          pvalue = pvalue_RNAseq)
  print("Analyse_transcription: DONE")
  Genes_network <- list()
  Genes_network[["UP"]] <- explore_genes_network(Genes_Differentially_expressed[["UP"]], gene_gene_network)
  Genes_network[["DOWN"]] <- explore_genes_network(Genes_Differentially_expressed[["DOWN"]], gene_gene_network)
  print("explore_genes_network: DONE")
  Overlap_bed_file <- find_overlap_between_data_and_chromatine_network(Filtered_DATA_Met, chromatine_BED_file)
  print("find_overlap_between_data_and_chromatine_network: DONE")
  print("Creation of the network in progress... May be long to compute!")
  Network <- Create_network(Overlap_bed_file, chromatine_network, List_of_promoter)
  print("Create_network: DONE")
  Network[["features"]] <- Transcription_features(Genes_network, Network[["features"]])
  Network[["network"]] <- Filter_network(Network[["network"]])
  print("DONE!")
  return(Network)
}


###############################################

##### Choose the comparison between 5:
Comparison <- data.frame("WT vs MUT" = "WT_vs_MUT", "WT vs HG" = "WT_vs_HG", "HG vs MUT" = "HG_vs_MUT", "WT vs Octyl" = "WT_vs_Octyl", "HG vs Octyl" = "HG_vs_Octyl")
number <- 1

##### Choose the logFC threshold for the transcriptome analysis
threshold <- 1.5

#### Choose the pvalue threshold for the transcriptome analysis
pvalue <- 0.1

#### Choose the pvalue threshold for the methylation analysis
pvalue_met <- 0.0001

logFC_treshold_met <- 0.3

Network <- From_data_to_network(DATA_MET = RawnoNABeta, 
                                Samplesheet_MET = Samplesheet, 
                                annotations_MET = anno_450, 
                                pvalue_met = pvalue_met, 
                                logFC_treshold = logFC_treshold_met,
                                number_of_the_comparison = number, 
                                DATA_RNAseq = eset, 
                                Samplesheet_RNAseq = SDRF, 
                                annotation_RNAseq = Annot, 
                                logFC_threshold_RNAseq = threshold, 
                                pvalue_RNAseq = pvalue, 
                                chromatine_BED_file = PCHiC_GRange, 
                                chromatine_network = pchic, 
                                List_of_promoter = List_Promoter, 
                                gene_gene_network = reactomeFI_network)

Analysed_methyl <- Analyze_methylation(RawnoNABeta, Samplesheet, anno_450)
Filtered_methyl <- Filter_Differential_Methylation(Analysed_methyl, pvalue_met, logFC_treshold_met)

filtered_network <- Filter_network(Network[["network"]][,c(1,5, 2, 3, 4, 6, 7, 8)])


write.csv(Network[["network"]][, c(1, 5)], file = paste0("~/Genes_network/Final_network/Results/CpGs_chromatine_Network_", 
                                                         Comparison[1,number], 
                                                         "_pvalue_", pvalue, 
                                                         "_logFC_", threshold, 
                                                         "_pval_met_", pvalue_met, 
                                                         "_logFCmet_", logFC_treshold_met,  
                                                         ".csv"), row.names = FALSE)
write.csv(Network[["features"]], file = paste0("~/Genes_network/Final_network/Results/Nodes_features_", 
                                               Comparison[1,number], 
                                               "_pvalue_", pvalue, 
                                               "_logFC_", threshold, 
                                               "_pval_met_", pvalue_met, 
                                               "_logFCmet_", logFC_treshold_met,  
                                               ".csv"), row.names = FALSE)
