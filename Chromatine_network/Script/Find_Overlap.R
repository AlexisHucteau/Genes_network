library("igraph")
library("GenomicRanges")
library("tidyverse")

# Importation of the data

setwd("Genes_network/Chromatine_network/")
getwd()
# PCHiC experiment
load("./DATA/pchic.RData")
pchic <- pchic[, c(1:10)]
# Coordinates of CpGs differently methylated
Cpgs_Coordinates <- read.csv("./DATA/CpGs_mapinfo.csv")
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

# For any reason, the bait ID : 19_37726761 and 7_1261811 are duplicated.

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

write.csv(Node_features, file = "Results/Node_features.csv", row.names = FALSE)
write.csv(CpGs_all_network[, c(1, 5)], file = "Results/CpGs_chromatine_Network.csv", row.names = FALSE)
