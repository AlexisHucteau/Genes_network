library("igraph")
library("GenomicRanges")
library("tidyverse")

setwd("~/Genes_network/Chromatine_network/")
getwd()

anno_450 <- read.csv("/home/alexis/Documents/Data_TCGA/Methylomes/HumanMethylation450_15017482_v1-2.csv", as.is = TRUE, skip = 7)

anno_450 <- anno_450[, c("CHR", "MAPINFO", "Name", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "UCSC_RefGene_Group")]

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



Genes_DE_UP <- read_tsv(file = "../Genes_Expansion/Networks_of_new_genes/WT_vs_MUT_UP/biogrid/WT_vs_MUT_UP_nodes.tsv")
Genes_network_UP <- as.array(Genes_DE_UP$gene)
Genes_DE_DOWN <- read_tsv(file = "../Genes_Expansion/Networks_of_new_genes/WT_vs_MUT_DOWN/biogrid/WT_vs_MUT_DOWN_nodes.tsv")
Genes_network_DOWN <- as.array(Genes_DE_DOWN$gene)


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


find_connection_to_gene <- function(overlaps, vector_of_genes) {
  cpt <- 1
  gene_connection <- list()
  for (o in overlaps) {
    select_from_network <- sapply(str_split(o$gene, ";"), function(genes) {
      length(intersect(unlist(genes), vector_of_genes)) > 0
    })
    gene_connection[[cpt]] <- length(o[select_from_network,1])
    if (cpt%%10==0){
      print(cpt)
    }
    cpt <- cpt + 1
  }
  return(gene_connection)
}

Sample_full_analyze <- function(annotations, number_of_test, size_of_sample, vector_of_genes_up, vector_of_genes_down) {
  Sample_list <- Sample_function(annotations, size_of_sample = size_of_sample, number_of_sample = number_of_test)
  Overlaps <- find_overlap_sample_test_function(Sample_list)
  Promoter_distribution <- find_promoter_distribution(Overlaps)
  Gene_connection_up <- find_connection_to_gene(Overlaps, vector_of_genes_up)
  Gene_connection_down <- find_connection_to_gene(Overlaps, vector_of_genes_down)
  res <- list()
  res[["Promoter_distribution"]] <- Promoter_distribution
  res[["Gene_connection_up"]] <- Gene_connection_up
  res[["Gene_connection_down"]] <- Gene_connection_down
  return(res)
}

result <- Sample_full_analyze(anno_450, size_of_sample = 13432, number_of_test = 100, vector_of_genes_up = Genes_network_UP, vector_of_genes_down = Genes_network_DOWN)

mean(unlist(result[["Promoter_distribution"]]))
mean(unlist(result[["Gene_connection_up"]]))
mean(unlist(result[["Gene_connection_down"]]))
