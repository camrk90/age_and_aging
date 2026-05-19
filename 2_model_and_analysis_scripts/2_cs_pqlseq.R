#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mem=50G 
#SBATCH --array=1-21

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

#### THIS SCRIPT RUNS PQLseq FOR THE CROSS SECTIONAL DATA ####

library(tidyverse)
library(PQLseq2)
setwd("/scratch/ckelsey4/Cayo_meth/cross_models")

#Generate function--------------------------------------------------------------
run_pqlseq<- function(pheno, covariates, type){
  
  mod<- pqlseq2(Y = meth, x = pheno, 
                K = kinship, W = covariates, 
                lib_size = cov, model="BMM")

  
  mod<- mod %>%
    filter(converged == TRUE) %>%
    mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
    relocate(fdr, .after = pvalue) %>%
    dplyr::select(-c(converged, elapsed_time))
  
  colnames(mod)<- c("outcome", "n", paste(names(mod[,3:length(mod)]), type, sep = "_"))
  
  return(mod)
  
}

#Import metadata----------------------------------------------------------------
blood_metadata<- read.table("/scratch/ckelsey4/Cayo_meth/blood_metadata_full.txt", sep = "\t", header = T)
long_data<- read.table("/scratch/ckelsey4/Cayo_meth/long_data_adjusted.txt")

long_ids<- unique(long_data$monkey_id)

overlap_lids<- blood_metadata[blood_metadata$monkey_id %in% long_ids,]

overlap_lids<- overlap_lids %>%
  group_by(monkey_id) %>%
  sample_n(1)

lids_to_remove<- long_data[!long_data$lid_pid %in% overlap_lids$lid_pid,]
blood_metadata<- blood_metadata[!blood_metadata$lid_pid %in% lids_to_remove$lid_pid,]

blood_metadata<- blood_metadata %>%
  filter(age_at_sampling > 1)

rm(lids_to_remove);rm(long_data);rm(long_ids)

#Import kinship matrix----------------------------------------------------------
kinship<- readRDS("/scratch/ckelsey4/Cayo_meth/full_kin_matrix")

#Subset and rearrange kinship rows and cols to match metadata
kinship<- kinship[blood_metadata$lid_pid, blood_metadata$lid_pid]

#Import m/cov rds------------------------------------------------------------
# load region lists that have been filtered for 5x coverage in 90% of samples
regions_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_filtered_cs.rds")
regions_m<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_m_filtered_cs.rds")

#Filter metadata to lids in regions list
blood_metadata<- blood_metadata[blood_metadata$lid_pid %in% colnames(regions_cov[[1]]),]

regions_cov<- lapply(names(regions_cov), function(x){
  regions_cov<- subset(regions_cov[[x]], select=blood_metadata$lid_pid)
  return(regions_cov)
})

regions_m<- lapply(names(regions_m), function(x){
  regions_m<- subset(regions_m[[x]], select=blood_metadata$lid_pid)
  return(regions_m)
})

names(regions_cov)<- 1:21 #turn all chroms into integers (X = 21)
names(regions_m)<- 1:21 #turn all chroms into integers (X = 21)

#Check metadata lids match the lids (cols) of a random chromosome
if (all.equal(blood_metadata$lid_pid, colnames(regions_cov[[runif(1, 1, 21)]]))) {
  
  #Model Vectors for lme4-------------------------------------------------------
  cov<- regions_cov[[SAMP]]
  meth<- regions_m[[SAMP]]
  
  ###################################
  #####        Run PQLseq       #####
  ###################################
  #Run PQLseq-------------------------------------------------------------------
  #Generate model matrix
  cs_matrix<- model.matrix(~ age_at_sampling + individual_sex + university, data = blood_metadata)
  
  vars <- c("age_at_sampling", "individual_sexM")
  
  cs_model <- lapply(setNames(vars, vars), function(i) {
    
    cs_phenotype <- cs_matrix[, i]
    cs_covariates <- cs_matrix[, setdiff(colnames(cs_matrix), i)]
    
    run_pqlseq(cs_phenotype, cs_covariates)
    
  })
  
  #Save pqlseq model
  saveRDS(cs_model, paste("cs_pqlseq2_age", SAMP, sep = "_"))
  
} else {
  
  print("blood_metadata lids did not match cov matrix lids")
  
}


