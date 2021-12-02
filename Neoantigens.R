library(biomaRt)
library(tidyverse)
library(data.table)
library(ggpubr)

# write files in proper format for netMHC
mutations <- fread("COSINR_Mutations.maf")

non_syn_snv <- mutations %>% filter(Variant_Type == "SNP",
                                    HGVSp_Short != "",
                                    Variant_Classification != "Silent")

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "uswest")
peptides <- getSequence(id = unique(non_syn_snv$Annotation_Transcript_NoVers),
            type = "ensembl_transcript_id",
            seqType = "peptide", 
            mart = mart)

non_syn_snv <- left_join(non_syn_snv, peptides)
non_syn_snv <- non_syn_snv %>% mutate(peptide_17mer = substr(peptide, AA_Position - 8, AA_Position + 8))

peptide_9mers <- bind_rows(lapply(1:nrow(non_syn_snv), function(i) {
  list_of_9mers <- c()
  for(j in 1:9) {
    list_of_9mers <- c(list_of_9mers, substr(peptide_9mers$peptide_17mer[i], j, j + 8))
  }
  data.frame(id = peptide_9mers$sample_id[i],
             peptide_sequence = list_of_9mers)
}))

peptide_9mers <- peptide_9mers %>% filter(peptide_sequence != "",
                                          nchar(peptide_sequence) == 9)

# import HLA results from PolySolver
hla_result_path <- "HLA_Results/winners/"
hla_result_files <- list.files(path = hla_result_path)
hla_results <- bind_rows(lapply(hla_result_files, function(x) fread(paste0(hla_result_path, x), header = F) %>%
                                  mutate(sample = str_extract(x, "SP-[0-9]+")))) %>%
  dplyr::select(sample, everything())
names(hla_results) <- c("sample_id", "hla_type", "allele1", "allele2")


fpkm <- fread("COSINR_FPKM.txt")

# check to make sure all HLA alleles are expressed
hla <- fpkm %>% filter(external_gene_name %in% c("HLA-A", "HLA-B", "HLA-C"))

# Use these as inputs to netMHC
sample_list_netmhc <- peptide_9mers %>% 
  dplyr::select(sample_id) %>% 
  distinct() %>%
  as.data.frame() %>%
  drop_na(sample_id)

hla_netmhc <- split(hla_netmhc, hla_netmhc$sample)
peptides_netmhc <- split(peptides_netmhc, peptides_netmhc$sample)

# after netMHC is run (see run_netMHC.sh)
netmhc_results_path <- "netMHCpan-4.1/Results/"
netmhc_results_files <- list.files(path = netmhc_results_path, pattern = "xls")
netmhc_results <- bind_rows(lapply(netmhc_results_files, function(x) {
  temp = fread(paste0(netmhc_results_path, x)) %>% as.data.frame()
  n = unlist(temp[2,c(1:2,4:7)])
  alleles = temp[1,][temp[1,] != ""]
  temp = temp[-1:-2,]
  
  temp = bind_rows(lapply(1:((ncol(temp) - 5) / 4), function(y) {
    data.frame(setNames(temp[,c(1:2, (4*y):(4*y+3))], c(1:6)))
  }))

  names(temp) = n
  temp$sample_id = substr(x, 1, 6)
  temp$HLA = sort(rep(alleles, nrow(temp) / length(alleles)))
  temp
  }))

netmhc_results_pass <- netmhc_results %>% dplyr::select(sample_id, HLA, Peptide, EL_Rank) %>%
  left_join(peptide_9mers) %>%
  dplyr::rename(percent_rank = EL_Rank,
                hla_allele_formatted = HLA,
                peptide_sequence = Peptide) %>%
  mutate(percent_rank = as.numeric(percent_rank)) %>%
  left_join(non_syn_snv) %>%
  left_join(transcript_names) %>%
  left_join(fpkm) %>%
  distinct() %>% filter(FPKM > 0, percent_rank < 2)

neoant_neopep_counts <- netmhc_results_pass %>% group_by(sample_id) %>%
  dplyr::summarize(n_neopeptides = length(unique(peptide_sequence)),
                   n_neoantigens = length(unique(HGVSg))) %>%
  ungroup()
