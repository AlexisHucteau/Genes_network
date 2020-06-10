# JOBIM_Network
Connect Methylomes and Transcriptomes

In IDH1m cells, the DNA is up methylated. By looking at CpGs we can look at the specific localization of thoses methylations. The impact of methylations on genes is knwown to have a down regulation of their transcription. By connecting transcriptomics and methylomics data, we can see the direct impact of methylation due to the mutation.

To connect the two type of data, we first have to analyze both separetly. Then we have to connect the specific CpGs found to genes that are closed to. That means, looking at 3D conformation of the chromatine.

## **STEP 1**

### Analyze of the transcriptomes

### **STEP 1.1**

##### *Data 1*

RNAseq from HL60 cell lines, a promyeloblast cell line. Somes are WT, others are IDH1 mutated, others are WT but treated with 2HG (the direct result of reaction by IDH1m) and others are WT treated by Octyl (a small molecule involoved in the transport of HG into the cells)

**DATA :** Transcriptomes/DATA_Transcriptomes/DATA

**Samplesheet :** Transcriptomes/Samplesheet_Transcriptomes

**Scripts :** Transcriptomes/Scripts_Transcriptomes/Transcriptomic_Analyze.Rmd

**Results :** Transcriptomes/Results_Transcriptomes

**Gene_Ontology :** Transcriptomes/Results_Transcriptomes/Gene_Ontology_enrichment

### **STEP 1.2**

The list of genes that are differently expressed is poor for some comparision. To extend thoses lists we have to find the genes that are linked functionnaly to them. Thanks to networks like Biogrid or Reactome, we expend the data.

##### *Data 2*

List of genes for each comparison (WT/Mut, WT/WT+HG, WT+HG/MUT, WT/WT+Octyl, WT+Octyl/WT+HG) and in each direction (up or down regulation)

**DATA :** Transcriptomes/DATA_Transcriptomes/Results_Transcriptomes/[Comparison]/Genes_values/

With logFC threshold = 1.5 and P-value = 0.1

- 4 genes WT vs MUT UP regulation
- 3 genes WT vs MUT DOWN regulation


- 14 genes WT vs WT + HG UP regulation
- 5 genes WT vs WT + HG DOWN regulation


- 11 genes WT + HG vs MUT UP regulation
- 21 genes WT + HG vs MUT DOWN regulation


- 18 genes WT vs WT + Octyl UP regulation
- 2 genes WT vs WT + Octyl DOWN regulation

With logFC threshold = 0.75 and P-value = 0.1

- 8 genes WT + HG vs WT + Octyl UP regulation
- 6 genes WT + HG vs WT + Octyl DOWN regulation


#### *Data 3*

Networks from Biogrid and Reactome

**DATA :**
- Genes_Expension/Networks/BIOGRID-ORGANISM-Homo_sapiens-3.5.185.tab3.txt
- Genes_Expension/Networks/FIsInGene_020720_with_annotations.tsv

**Results :** Genes_Expension/Networks_of_new_genes

- 313 genes WT vs MUT UP regulation
- 11 genes WT vs MUT DOWN regulation


- 4 genes WT vs WT + HG UP regulation
- 3 genes WT vs WT + HG DOWN regulation


- 255 genes WT + HG vs MUT UP regulation
- 45 genes WT + HG vs MUT DOWN regulation


- 278 genes WT vs WT + Octyl UP regulation
- 0 genes WT vs WT + Octyl DOWN regulation

With logFC threshold = 0.75 and P-value = 0.1

- 24 genes WT + HG vs WT + Octyl UP regulation
- 75 genes WT + HG vs WT + Octyl DOWN regulation

**Genes Networks Visualization :** Genes_Expension/Genes_Networks_Visualizations/Genes_Networks_Visualization.cys

**Gene Ontology Visualization :** Genes_Expension/Gene_Ontology/Visualizations/Gene_Ontology.cys

**Gene Ontology Table :** Genes_Expension/Gene_Ontology/Tables/

**Scripts :** Genes_Expension/Scripts_Expension/Expend_genes.Rmd

## **STEP 2**

### Analyze of the methylomes

##### *Data 2*

**DATA :** Methylomes/DATA_Methylomes/DATA

**Samplesheet :** Methylomes/Samplesheet_Methylomes

**Scripts :** Methylomes/Scripts_Methylomes

**Results :** Methylomes/Results_Methylomes
