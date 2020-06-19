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
  print(paste0("The number of CpGs differentially methylated is ", length(High_contrast[,1])))
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
  print(paste0("The number of genes up regulated is ", length(Genes_values[["UP"]]))  )
  print(paste0("The number of genes down regulated is ", length(Genes_values[["DOWN"]])))
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

Filter_network <- function(network) {
  network <- network[network$CpGs_bait | network$CpGs_oe | network$OePromoter,]
  network <- network[-which(network$Name_oe == "."),]
  return(network)
}

Sample_function <- function(methylation_annotations, number_of_sample=100, size_of_sample) {
  cpt <- 1
  List_of_test <- list()
  for (i in seq(number_of_sample)) {
    tmp <- sample(methylation_annotations$Name, size = size_of_sample, replace = FALSE)
    List_of_test[[i]] <- methylation_annotations[methylation_annotations$Name %in% tmp,]
    List_of_test[[i]] <- List_of_test[[i]][str_detect(List_of_test[[i]]$Name, "cg"),]
    if (cpt%%10==0){
      print(cpt)
    }
    cpt <- cpt + 1
  }
  return (List_of_test)
}

find_overlap_sample_test_function <- function(list_of_Sample) {
  cpt <- 1
  match_hit <- list()
  for (sample in list_of_Sample) {
    CpGs_Ranges <- GRanges(
      seqnames = sample$CHR,
      ranges = IRanges(start = sample$MAPINFO, end = sample$MAPINFO + 1),
      genes = sample$UCSC_RefGene_Name,
      position = sample$UCSC_RefGene_Group
    )
    overlaps <- findOverlaps(PCHiC_GRange, CpGs_Ranges)
    match_hit[[cpt]] <- data.frame(mcols(PCHiC_GRange[queryHits(overlaps)]), as.data.frame(mcols(CpGs_Ranges)[subjectHits(overlaps), ]), stringsAsFactors = T)
    match_hit[[cpt]]$CpGs <- TRUE
    if (cpt%%10==0){
      print(cpt)
    }
    cpt <- cpt + 1
  }
  return(match_hit)
}

find_promoter_distribution <- function(overlaps) {
  cpt <- 1
  promoter_distribution <- list()
  for (o in overlaps) {
    promoter_distribution[[cpt]] <- length(o[str_detect(o$position, "TSS"),1])
    if (cpt%%10==0){
      print(cpt)
    }
    cpt <- cpt + 1
  }
  return(promoter_distribution)
}

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
  res <- list()
  res[["network"]] <- CpGs_all_network
  res[["features"]] <- Node_features
  return(res)
}

Create_samples_network <- function(matchit, pchic, list_of_promoter) {
  res <- list()
  cpt <- 1
  for (match_hit in matchit) {
    res[[cpt]] <- Create_network(match_hit, pchic, list_of_promoter)
    if (cpt%%10==0){
      print(cpt)
    }
    cpt <- cpt + 1
  } 
  return(res)
}

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
Transcription_features <- function(gene_list, network) {
  cpt <- 1
  res <- list()
  for (sample in network) {
    selected_UP <- Transcription_feature(gene_list[["UP"]], sample$features)
    selected_DOWN <- Transcription_feature(gene_list[["DOWN"]], sample$features)
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
    res[[cpt]] <- final_feature
    print(paste0("The number of genes from the DE network connected to the final network is ",  length(res[[cpt]][which(!str_detect(res[[cpt]], "Not"))])))
    res[[cpt]] <- length(res[[cpt]][which(!str_detect(res[[cpt]], "Not"))])
    cpt <- cpt + 1
  }
  return(res)
}

Compute_promoter_in_network <- function(features) {
  cpt <- 1
  res <- list()
  for (sample in features) {
    tmp <- sample$features
    res[[cpt]] <- length(tmp[tmp$Promoter,1])
    print(paste0("The number of genes in the network final is ",  res[[cpt]]))
    if (cpt%%10==0){
      print(cpt)
    }
    cpt <- cpt + 1
  }
  return(res)
}


From_Sampletest_to_network <- function(methylation_annotations, 
                                       number_of_sample=100, 
                                       size_of_sample,
                                       number_of_the_comparison = number,                                  
                                       DATA_RNAseq, 
                                       Samplesheet_RNAseq, 
                                       annotation_RNAseq, 
                                       logFC_threshold_RNAseq, 
                                       pvalue_RNAseq, 
                                       chromatine_BED_file = PCHiC_GRange, 
                                       List_of_promoter = List_Promoter, 
                                       gene_gene_network = reactomeFI_network) {
  Filtered_DATA_Met <- Sample_function(methylation_annotations, number_of_sample, size_of_sample)
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
  Overlaps <- find_overlap_sample_test_function(Filtered_DATA_Met)
  Networks <- Create_samples_network(Overlaps, pchic, List_of_promoter)
  Promoter_distribution <- find_promoter_distribution(Overlaps)
  res <- list()
  res[["Promoter_distribution"]] <- Promoter_distribution
  res[["Networks"]] <- Networks
  # res[["Gene_connection_down"]] <- Gene_connection_down
  return(res)
}

result <- From_Sampletest_to_network(methylation_annotations = anno_450, 
                                     number_of_sample=100, 
                                     size_of_sample = 256,
                                     number_of_the_comparison = 1,                                  
                                     DATA_RNAseq = eset, 
                                     Samplesheet_RNAseq = SDRF, 
                                     annotation_RNAseq = Annot, 
                                     logFC_threshold_RNAseq = 1.5, 
                                     pvalue_RNAseq = 0.1, 
                                     chromatine_BED_file = PCHiC_GRange, 
                                     List_of_promoter = List_Promoter, 
                                     gene_gene_network = reactomeFI_network)


Genes_Differentially_expressed <- Analyse_transcription(number_of_the_comparison = 1, 
                                                        DATA = eset, 
                                                        Samplesheet = SDRF, 
                                                        annotation = Annot, 
                                                        logFC_threshold = 1.5, 
                                                        pvalue = 0.1)
print("Analyse_transcription: DONE")
Genes_network <- list()
Genes_network[["UP"]] <- explore_genes_network(Genes_Differentially_expressed[["UP"]], reactomeFI_network)
Genes_network[["DOWN"]] <- explore_genes_network(Genes_Differentially_expressed[["DOWN"]], reactomeFI_network)
promoters <- Compute_promoter_in_network(result[["Networks"]])
Features <- Transcription_features(Genes_network, result[["Networks"]])

mean(unlist(result[["Promoter_distribution"]]))
mean(unlist(Features))

