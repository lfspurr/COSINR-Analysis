# "Exome signature extraction"

library(data.table)
library(tidyverse)
library(NMF)
library(BSgenome.Hsapiens.UCSC.hg38)
devtools::load_all("~/Documents/Research/Pitroda_RT-IO/MutationalPatterns/")

# keep only samples with at least 20 mutations
maf <- fread("~/Documents/Research/Pitroda_RT-IO/Data/COSINR_Mutations_Final_092121.oxog.OncoKb.CCF.maf") %>% 
  group_by(Tumor_Sample_Barcode) %>%
  mutate(mut_count = n()) %>%
  filter(mut_count >= 20) %>% 
  ungroup()

snps <- maf %>% 
  filter(Variant_Type == "SNP", 
         !(Chromosome %in% c("__UNKNOWN__")))

mutations_GR <- GRanges(
  seqnames = snps$Chromosome,
  ranges = IRanges(snps$Start_Position, snps$End_Position),
  ref = snps$Reference_Allele,
  alt = snps$Tumor_Seq_Allele2,
  sample = snps$Tumor_Sample_Barcode,
)
mutations_GR@seqinfo@genome <- rep("hg38", length(mutations_GR@seqinfo@genome))
mutations_GR <- split(mutations_GR, as.factor(mutations_GR$sample))

type_occurrences <- mut_type_occurrences(mutations_GR, "BSgenome.Hsapiens.UCSC.hg38")

mut_mat <- mut_matrix(mutations_GR, "BSgenome.Hsapiens.UCSC.hg38")


# Fit to all signatures
signatures <- get_known_signatures(muttype = "snv", incl_poss_artifacts = F)
fit <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.02)


candidate_sigs <- paste("SBS", c("1", "2", "3", "4", "5", "6", "7a", "7b", "7c", "8", "9", "10c", "10d", "13", "15", "17a", "17b", "18", "28", "29", "30", "40", "44", "86", "87", "89"), sep = "")
signatures <- signatures[,candidate_sigs[candidate_sigs %in% colnames(signatures)]]
# https://townsend-lab-yale.github.io/cancereffectsizeR/articles/cosmic_cancer_type_note.html

strict_refit <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.02)
fit_res_strict <- strict_refit$fit_res$contribution
fit_res_strict <- fit_res_strict[which(rowSums(fit_res_strict ) > 0),]

contri_boots <- fit_to_signatures_bootstrapped(mut_mat,
                                               signatures,
                                               n_boots = 50,
                                               method = "strict",
                                               max_delta = 0.02,
)
