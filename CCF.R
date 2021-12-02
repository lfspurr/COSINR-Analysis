library(data.table)
library(tidyverse)

facets_summary <- fread("COSINR_FACETS_Purity_Ploidy.txt")
maf_annotated <- fread("COSINR_Mutations_Final.maf")
facets_gene_level <- fread("COSINR_FACETS_GeneLevelCNVCalls.txt")

# Functions adapted from: https://github.com/mskcc/facets-suite/blob/master/R/ccf-annotate-maf.R, based on PMID 25877892
# Estimate most likely CCF given observed VAF, purity and local ploidy
estimate_ccf <- function(purity,
                        total_copies,
                        mutant_copies,
                        t_alt_count,
                        t_depth) {
  
  
  ccfs = seq(0.001, 1, 0.001)
  expected_vaf  = function(ccf, purity, total_copies) {
    purity * ccf * mutant_copies / (2 * (1 - purity) + purity * total_copies)
  }
  
  probs = sapply(ccfs, function(c) {
    stats::dbinom(t_alt_count, t_depth, expected_vaf(c, purity, total_copies))
  })
  probs = probs / sum(probs)
  
  ccf_max = which.max(probs)
  if (identical(ccf_max, integer(0))) ccf_max = NA
  ccf_half_max = which(probs > max(probs) / 2)
  ccf_lower = max(ccf_half_max[1] - 1, 1) # closest ccf value before half-max range (within 0-1 range)
  ccf_upper = min(ccf_half_max[length(ccf_half_max)] + 1, length(ccfs)) # closest ccf value after half-max range (within 0-1 range)
  if (is.na(purity)) ccf.upper = NA 
  ccf_max = ccf_max / length(ccfs)
  ccf_lower = ccf_lower / length(ccfs)
  ccf_upper = ccf_upper / length(ccfs)
  prob95 = sum(probs[950:1000])
  prob90 = sum(probs[900:1000])
  
  data.frame(ccf_max = ccf_max, ccf_lower = ccf_lower, ccf_upper = ccf_upper, prob95 = prob95, prob90 = prob90)
}

# Estimate mutant copy number, given observed VAF, purity, and local ploidy
expected_mutant_copies <- function(t_var_freq,
                                  total_copies,
                                  purity) {
  
  if (is.na(total_copies)) {
    NA_real_
  } else {
    if (total_copies == 0) total_copies = 1
    mu = t_var_freq * (1 / purity) * (purity * total_copies + (1 - purity) * 2)
    alt_copies = ifelse(mu < 1, 1, abs(mu))
    round(alt_copies)
  }
}

maf_ccf <- cbind(maf_annotated, bind_rows(lapply(1:nrow(maf_annotated), function(i) estimate_ccf(maf_annotated$purity[i],
                                           maf_annotated$tcn[i],
                                           mutant_copies = expected_mutant_copies(maf_annotated$Allele_Fraction[i], maf_annotated$tcn[i], maf_annotated$purity[i]),
                                           maf_annotated$t_alt_count[i],
                                           maf_annotated$t_alt_count[i] + maf_annotated$t_ref_count[i]))))

