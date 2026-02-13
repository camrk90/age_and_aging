#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mem=10G 
#SBATCH --array=1-21

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(tidyverse)
library(PQLseq2)
setwd("/scratch/ckelsey4/Cayo_meth/glmer_model_compare")

#Import metadata----------------------------------------------------------------
long_data<- read.table("/scratch/ckelsey4/Cayo_meth/long_data_adjusted.txt")

long_data<- long_data %>%
  group_by(monkey_id) %>%
  mutate(n = n()) %>%
  ungroup()

long_data<- long_data %>%
  filter(age_at_sampling > 1) %>%
  filter(n > 1) %>%
  dplyr::rename(perc_unique = unique) %>%
  drop_na() %>%
  arrange(lid_pid)

#Import kinship matrix----------------------------------------------------------
kinship<- readRDS("/scratch/ckelsey4/Cayo_meth/full_kin_matrix")

#Subset and rearrange kinship rows and cols to match metadata
kinship<- kinship[long_data$lid_pid, long_data$lid_pid]

#Import m/cov rds------------------------------------------------------------
# load region lists that have been filtered for 5x coverage in 90% of samples
regions_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_filtered.rds")
regions_m<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_m_filtered.rds")

#Filter metadata to lids in regions list
long_data<- long_data[long_data$lid_pid %in% colnames(regions_cov[[1]]),]

regions_cov<- lapply(names(regions_cov), function(x){
  regions_cov<- subset(regions_cov[[x]], select=long_data$lid_pid)
  return(regions_cov)
})

regions_m<- lapply(names(regions_m), function(x){
  regions_m<- subset(regions_m[[x]], select=long_data$lid_pid)
  return(regions_m)
})

names(regions_cov)<- 1:21 #turn all chroms into integers (X = 21)
names(regions_m)<- 1:21 #turn all chroms into integers (X = 21)

#Check metadata lids match the lids (cols) of a random chromosome
if (all.equal(long_data$lid_pid, colnames(regions_cov[[runif(1, 1, 21)]]))) {
  
  #Model Vectors for lme4---------------------------------------------------------
  cov<- regions_cov[[SAMP]]
  meth<- regions_m[[SAMP]]
  
  ###################################
  #####        Run PQLseq       #####
  ###################################
  #Run PQLseq for within_age-------------------------------------------------------------
  #Generate model matrix
  predictor_matrix<- model.matrix(~ within.age + mean.age + individual_sex + perc_unique, data = long_data)
  w.age_phenotype<- predictor_matrix[, 2]
  w.age_covariates<- predictor_matrix[, 3:5]
  
  #Run pqlseq model
  w.age_pqlseq2_model<- pqlseq2(Y = meth, x = w.age_phenotype, 
                                K = kinship, W = w.age_covariates, 
                                lib_size = cov, model="BMM")
  
  #Run PQLseq for chronological age---------------------------------------------
  #Generate model matrix
  predictor_matrix<- model.matrix(~ age_at_sampling + individual_sex + perc_unique, data = long_data)
  agechron_phenotype<- predictor_matrix[, 2]
  agechron_covariates<- as.matrix(predictor_matrix[, 3:4])
  
  #Run pqlseq model
  agechron_pqlseq2_model<- pqlseq2(Y = meth, x = agechron_phenotype, 
                                   K = kinship, W = agechron_covariates, 
                                   lib_size = cov, model="BMM")
  
  #Save pqlseq models
  saveRDS(w.age_pqlseq2_model, paste("within", "age", "no", "singles", SAMP, sep = "_"))
  saveRDS(agechron_pqlseq2_model, paste("agechron", "no", "singles", SAMP, sep = "_"))
  
} else {
  
  print("long_data lids did not match cov matrix lids")
  
}



