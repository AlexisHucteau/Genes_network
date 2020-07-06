#########Metylation part

Samplesheet <- read.csv("~/Documents/Methylations_IDH1m/Wiehle/E-MTAB-6059.raw.1/Phenotype.csv")
Phenotype1 <- ifelse(str_detect(Samplesheet$Characteristics.phenotype., "wild"), "WT", "Mut")
Phenotype2 <- ifelse(str_detect(Samplesheet$Factor.Value.compound., "hydro"), "HG", "None")
Phenotype <- as.factor(paste(Phenotype1, Phenotype2, sep = "."))

DATA_met <- Import_Methylation_DATA(pathfile = "~/Documents/Methylations_IDH1m/Wiehle/E-MTAB-6059.raw.1/", Illumina_type_array = "EPIC", Samplesheet = Samplesheet, Phenotype = Phenotype)


fit_design <- Fit_and_design(comparison = Phenotype, DATA = DATA_met[["BMIQ"]])
Comparison_matrix <- makeContrasts(
  WTvsMut = Mut.None - WT.None,
  levels = fit_design[["design"]]
)

Analysed_methylation <- Analyze_methylation(DATA = DATA_met[["BMIQ"]], annotations = DATA_met[["anno"]], comparison = Phenotype, comparison_matrix = Comparison_matrix, fit_design = fit_design, Condition_1 = "WT.None", Condition_2 = "Mut.None")

Filtered_methylation <- Filter_Differential_Methylation(Analysed_methylation, 0.001, 0.3)

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

matchit_met <- find_overlap_between_data_and_chromatine_network(Filtered_methylation, network = PCHiC_GRange)

Met_network <- Create_network(matchit_met, network = pchic, list_of_promoter = List_Promoter)

#Transcription part

DATA_trscr <- Import_Transcription_DATA(DATA_dir = "~/Genes_network/Transcriptomes/DATA_Transcriptomes/", Sample_file = "~/Genes_network/Transcriptomes/Samplesheet_Transcriptomes/Samplesheet.csv")

Phenotype_Trscr <- paste(DATA_trscr[["Samplesheet"]]$Characteristics.genotype., DATA_trscr[["Samplesheet"]]$Characteristics.treatment., sep = ".")

fit_design_trans <- Fit_and_design(comparison = Phenotype_Trscr, DATA = DATA_trscr[["DATA"]])

Annot <- data.frame(ACCNUM = sapply(contents(hugene20sttranscriptclusterACCNUM), paste, collapse = ", "),
                    SYMBOL = sapply(contents(hugene20sttranscriptclusterSYMBOL), paste, collapse = ", "),
                    DESC = sapply(contents(hugene20sttranscriptclusterGENENAME), paste, collapse = ", "),
                    ENTREZID = sapply(contents(hugene20sttranscriptclusterENTREZID), paste, collapse = ", "))

Genes_DE <- Analyse_transcription(number_of_the_comparison = 1, DATA = DATA_trscr[["DATA"]], Samplesheet = DATA_trscr[["Samplesheet"]], annotation = Annot, logFC_threshold = 0.75, pvalue = 0.1)

biogrid_data <- read_tsv("~/Genes_network/Genes_Expansion/Networks/BIOGRID-ORGANISM-Homo_sapiens-3.5.185.tab3.txt") %>%
  filter(`Organism Name Interactor A` == "Homo sapiens" & `Organism Name Interactor B` == "Homo sapiens")
reactome_FI <- read_tsv("~/Genes_network/Genes_Expansion/Networks/FIsInGene_020720_with_annotations.tsv")

biogrid_network <- graph_from_data_frame(biogrid_data, directed = F)
reactomeFI_network <- graph_from_data_frame(reactome_FI, directed = F)

Up_genes <- explore_genes_network(Genes_DE[["UP"]], biogrid_network)
Down_genes <- explore_genes_network(Genes_DE[["DOWN"]], biogrid_network)

library("biomaRt")

Up_genes <- add_genes_coordinates(Up_genes, Annot)
Down_genes <- add_genes_coordinates(Down_genes, Annot)

matchit_UP <- find_overlap_between_genes_and_chromatine_network(Up_genes, PCHiC_GRange)
matchit_DOWN <- find_overlap_between_genes_and_chromatine_network(Down_genes, PCHiC_GRange)

matchit <- rbind(matchit_UP, matchit_DOWN)
matchit <- unique(matchit)

Network_genes <- Create_network_genes(matchit, pchic)
Network_genes <- Network_genes[,c(11, 12)]

Network <- Combine_networks(Network_genes, Met_network[["network"]], matchit_UP, matchit_DOWN, matchit_met, List_Promoter)

nodes_attributes <- add_nodes_attributes(Network, pchic, matchit, matchit_met) 
  
