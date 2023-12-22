
# We will use the following files:

# data_clinical_patient.txt, data_mrna_seq_v2_rsem.txt, data_mutations.txt and data_cna.txt

clinical = read.delim("data_clinical_patient.txt")

rnaseq = read.delim("data_mrna_seq_v2_rsem.txt")


# in this assignment we will delete the genes for which there's more than one Hugo Symbol
# These are typically genes with no Hugo Symbol ("" as an entry) or pseudogenes.

# This is more for simplicity.If you keep your analysis would still be correct so no worries.

keep = !duplicated(rnaseq[,1])

rnaseq = rnaseq[keep,]

# set rownames of rnaseq to hugo symbols

rownames(rnaseq)  = rnaseq[,1]

# Read CNA Data

cna = read.delim('data_cna.txt')

# find ERBB2 in cna

erbb2_indx = which(cna[,1] == 'ERBB2')

# Plot histogram to visualize explore the data.

hist(as.numeric(cna[erbb2_indx,-c(1,2)]))

# match patients in rnaseq to patients in cna.

rna_cna_id = which(is.element(colnames(rnaseq[,-c(1,2)]), colnames(cna[,-c(1,2)])))

# select only the rna cases which have cna data.

rna_cna_sub = rnaseq[,2+rna_cna_id]

# check all patients in rna_can_sub are in cna

no_pats_in_rna_cna_sub_and_cna = sum(is.element(colnames(rnaseq[,2+rna_cna_id]), colnames(cna[,-c(1,2)]))) 

# sanity check.This will print an error if the result is not the same.

sanity_check = no_pats_in_rna_cna_sub_and_cna == dim(rna_cna_sub)[2]

# Pre-allocate memory for ERBB2

meta_erbb2 = matrix(0,length(rna_cna_id),1)

for (i in 1:length(rna_cna_id)){
  # access the colnames of i
  col_i = colnames(rna_cna_sub)[i]
  # get the index in cna for the same patient
  col_cna = which(colnames(cna)==col_i)
  # store if they're amplified.
  meta_erbb2[i,] = 1*(cna[erbb2_indx,col_cna]>0)
  
}

# This are some checks you can do to make sure your code worked.
# There's some more systematic checks you can do. See unit testing.


# simple checks to make sure. 

col_i = colnames(rna_cna_sub)[1]

col_cna = which(colnames(cna)==col_i)

# sanity check

(cna[erbb2_indx,col_cna]>0) == meta_erbb2[1,1]

# see now if a positive meta_erbb2 is amplified.

pos_example = which(meta_erbb2==1)[1]


col_i = colnames(rna_cna_sub)[pos_example]

col_cna = which(colnames(cna)==col_i)

# sanity check

(cna[erbb2_indx,col_cna]>0) == meta_erbb2[pos_example,1]

# botch checks should print true.

# add a title to the metadata.

colnames(meta_erbb2) = 'ERBB2Amp'

# transform into integers

rna_cna_sub = round(rna_cna_sub)

# Install DESeq2.


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install DeSeq2

BiocManager::install("DESeq2")

library(DESeq2)

# Build DESeq Object
# Based on the Session 15 R-tutorial: TCGA Differential Expression Analysis DESEQ2

dds <- DESeqDataSetFromMatrix(countData = round(rna_cna_sub),
                              colData = meta_erbb2,
                              design = ~ ERBB2Amp)
# Filter

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Normalize

dds <- DESeq(dds)

# Transform the data to visualize
rld <- vst(dds, blind=FALSE)

# Do Principal Components Analysis
plotPCA(rld, intgroup = "ERBB2Amp")

# Get Results

res <- results(dds)

# Summary
summary(res)
rownames(res) = rnaseq[keep,1]

# Top 10 Differentially Expressed Genes Ranked by Fold Change
signif = which(res$padj<0.05)
deg = res[signif,]
deg_sorted <- deg[order(deg$log2FoldChange, decreasing = TRUE), ]
top10_deg <- head(deg_sorted, 10)
print(top10_deg)

# KEGG Pathway Enrichment and Visualisation
# Based on Session 17 R-tutorial: Pathway Analysis
entrez_ids = rnaseq[keep,2]

entrez_all = entrez_ids[signif]
entrez_up = entrez_all[signif[deg[,2]>0.]]
entrez_down = entrez_all[signif[deg[,2]<0.]]

if (!require("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
if (!require("pathview", quietly = TRUE))
  BiocManager::install("pathview")

library(clusterProfiler)
library(enrichplot)
library(pathview)

enrich_paths <- enrichKEGG(gene = entrez_all,
                           organism = 'hsa',
                           pvalueCutoff = 0.05)
dotplot(enrich_paths) + ggtitle("KEGG Pathway Enrichment - enrichKEGG")

pathview(gene.data  = entrez_all,
         pathway.id = "hsa04110",
         species    = "hsa",
         limit      = list(gene=max(abs(entrez_all)), cpd=1))

# Gene Expression Cluster
pca_result <- plotPCA(rld, intgroup = "ERBB2Amp", returnData = TRUE)
pca_result$group <- factor(pca_result$group)

# DBSCAN Cluster
if (!require("dbscan", quietly = TRUE))
  install.packages("dbscan")
library(dbscan)
library(ggplot2)
clusters <- dbscan(pca_result[, c("PC1", "PC2")], eps = 3, MinPts = 10)
pca_result$cluster <- as.factor(clusters$cluster)

ggplot(pca_result, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(100 * attr(pca_result, "percentVar")[1]), "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pca_result, "percentVar")[2]), "% variance")) +
  ggtitle("PCA - DBSCAN Cluster") +
  theme_minimal()
