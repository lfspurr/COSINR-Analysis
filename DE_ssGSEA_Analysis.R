# "Differential expression analyses"

library(data.table)
library(tidyverse)
library(qvalue)
library(DESeq2)
library(fgsea)
library(GSVA)

source("Enrichment.R")


# Perform ssGSEA
pathways <- gmtPathways("h.all.v7.4.symbols.gmt")
ssgsea <- gsva(ssgsea_mat, pathways, method = "ssgsea", kcdf = "Poisson")


# Perform DE analysis
run_deseq <- function(clinical, counts, var, ref_level, gene_names, min_counts = 10, paired = T) {
  if(paired) design <- as.formula(paste0("~ Patient + ", var)) else design <- as.formula(paste0("~ ", var))
  de <- DESeqDataSetFromMatrix(countData = counts,
                               colData = clinical,
                               design = design)
  keep <- rowSums(counts(de)) >= min_counts
  de <- de[keep,]
  
  if(ref_level != levels(de@colData[,var])[1]) de@colData[,var] <- relevel(de@colData[,var], ref = ref_level)
  
  de <- DESeq(de)
  res_name <- resultsNames(de)
  res_name <- res_name[!grepl("Patient|Intercept", res_name)]
  
  resLFC <- lfcShrink(de, coef = res_name, type = "apeglm", format = "DataFrame")
  
  resLFC <- as.data.frame(resLFC)
  resLFC$ensembl_gene_id <- rownames(resLFC)
  resLFC <- left_join(resLFC, gene_names)
  cbind(resLFC, stat = results(de)$stat)
}
