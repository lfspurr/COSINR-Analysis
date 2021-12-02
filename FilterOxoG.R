# Resources:
# https://github.com/cpwardell/bin/blob/master/metalfox.py
# https://academic.oup.com/nar/article/41/6/e67/2902364

library(data.table)
library(tidyverse)
library(NMF)
library(BSgenome.Hsapiens.UCSC.hg38)
devtools::load_all("MutationalPatterns/")

# First, I will compute the FoxoG score as follows for all C>A/G>T mutations:
# FoxoG = F1R2 / (F1R2 + F2R1) or F2R1 / (F1R2 + F2R1)
# Numerator is F2R1 if C>A, F1R2 if G>T 
maf <- fread("COSINR_Mutations_Raw.maf") %>% 
  mutate(poss_artifact = ifelse((Reference_Allele == "C" & Tumor_Seq_Allele2 == "A") | 
                                  (Reference_Allele == "G" & Tumor_Seq_Allele2 == "T"), T, F),
         foxog = ifelse((Reference_Allele == "C" & Tumor_Seq_Allele2 == "A"), F2R1 / (F1R2 + F2R1),
                        ifelse((Reference_Allele == "G" & Tumor_Seq_Allele2 == "T"), F1R2 / (F1R2 + F2R1), NA)),
         artifact_score = -10 + (100 / 3) * foxog) 

# Then filter any mutations out where tumor_lod < -10 + (100 / 3) * FoxoG
oxog_filtered_muts <- maf %>% filter(TLOD < artifact_score & poss_artifact)

maf <- maf %>%
  filter(TLOD > artifact_score | !poss_artifact)

# QC the artifacts to make sure they look like OxoG (SBS45) mutations
arifact_qc <- oxog_filtered_muts %>% 
  filter(Variant_Type == "SNP", 
         !(Chromosome %in% c("__UNKNOWN__")))

artifact_GR <- GRanges(
  seqnames = arifact_qc$Chromosome,
  ranges = IRanges(arifact_qc$Start_Position, arifact_qc$End_Position),
  ref = arifact_qc$Reference_Allele,
  alt = arifact_qc$Tumor_Seq_Allele2,
  sample = arifact_qc$Tumor_Sample_Barcode,
)
artifact_GR@seqinfo@genome <- rep("hg38", length(artifact_GR@seqinfo@genome))
artifact_mat <- mut_matrix(artifact_GR, "BSgenome.Hsapiens.UCSC.hg38")

signatures <- get_known_signatures(incl_poss_artifacts = T)
most_sim <- t(cos_sim_matrix(artifact_mat, signatures)) %>% as.data.frame() %>% arrange(desc(My_sample)) %>% head(1)

plot_96_profile(artifact_mat) + labs(title = rownames(most_sim[1]), subtitle = paste("Cosine similarity =", round(most_sim$My_sample[1], 2)))
