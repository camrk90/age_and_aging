#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu

library(tidyverse)

rna_meta_full<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_metadata_wELA_18Jul25.rds")
rna_counts<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_counts_18Jul25.rds")

#Subset to control and truncate df to relevant variables
rna_meta<- rna_meta_full %>%
  filter(Stimulation == "H2O") %>%
  select(c(animal_ID, sex, trapped_age, trapping_ID, Sample_ID, Seq_batch))

#Generate within and between age
rna_meta<- rna_meta %>%
  group_by(animal_ID) %>%
  mutate(mean_age = mean(trapped_age)) %>%
  mutate(within_age = trapped_age - mean_age) %>%
  arrange(animal_ID, trapped_age)

#Make batch factor
rna_meta$Seq_batch<- gsub("batch", "", rna_meta$Seq_batch)
rna_meta$Seq_batch<- factor(rna_meta$Seq_batch, levels = c("1", "2", "3"))

#Subset counts to RIDs in metadata
rna_counts<- rna_counts[, colnames(rna_counts) %in% rna_meta$Sample_ID]

rna_pca<- prcomp(cor(t(rna_counts), use="pairwise.complete.obs"))

saveRDS(rna_pca, "/home/ckelsey4/Cayo_meth/rna_seq/rna_pca.rds")



