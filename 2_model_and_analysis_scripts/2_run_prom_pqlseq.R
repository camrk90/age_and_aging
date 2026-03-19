#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mem=50G 
#SBATCH --array=1-21

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(tidyverse)
library(PQLseq2)
setwd("/scratch/ckelsey4/Cayo_meth/glmer_model_compare")

#Generate function--------------------------------------------------------------
run_pqlseq<- function(pheno, covariates, type){
  
  mod_df<- pqlseq2(Y = meth, x = pheno, 
                   K = kinship, W = covariates, 
                   lib_size = cov, model="BMM")
  
  mod_df<- mod_df %>%
    filter(converged == TRUE) %>%
    mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
    relocate(fdr, .after = pvalue) %>%
    dplyr::select(-c(converged, elapsed_time))
  
  colnames(mod_df)<- c("outcome", "n", paste(names(mod_df[,3:length(mod_df)]), type, sep = "_"))
  
  return(mod_df)
  
}

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

#Load promoters-----------------------------------------------------------------
prom_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/prom_cov_filtered")
prom_cov<- prom_cov[c(1:21)]
prom_m<- readRDS("/scratch/ckelsey4/Cayo_meth/prom_m_filtered")
prom_m<- prom_m[c(1:21)]

#Filter metadata to lids in regions list
long_data<- long_data[long_data$lid_pid %in% colnames(prom_cov[[1]]),]

prom_cov<- lapply(names(prom_cov), function(x){
  prom_cov<- subset(prom_cov[[x]], select=long_data$lid_pid)
  return(prom_cov)
})

prom_m<- lapply(names(prom_m), function(x){
  prom_m<- subset(prom_m[[x]], select=long_data$lid_pid)
  return(prom_m)
})

names(prom_cov)<- 1:21 #turn all chroms into integers (X = 21)
names(prom_m)<- 1:21 #turn all chroms into integers (X = 21)

#Check metadata lids match the lids (cols) of a random chromosome
if (all.equal(long_data$lid_pid, colnames(prom_cov[[runif(1, 1, 21)]]))) {
  
  #Model Vectors for lme4-------------------------------------------------------
  cov<- prom_cov[[SAMP]]
  meth<- prom_m[[SAMP]]
  
  ###################################
  #####           Eq.1          #####
  ###################################
  #Run PQLseq for chronological age---------------------------------------------
  #Generate model matrix
  eq1_matrix<- model.matrix(~ age_at_sampling + individual_sex + perc_unique, data = long_data)
  eq1_phenotype<- eq1_matrix[, 2]
  eq1_covariates<- as.matrix(eq1_matrix[, 3:4])
  
  #Run pqlseq model
  agechron_pqlseq2_model<- run_pqlseq(eq1_phenotype, eq1_covariates, "chron_age")
  
  #Save pqlseq model
  saveRDS(agechron_pqlseq2_model, paste("prom", "pqlseq2", "agechron", SAMP, sep = "_"))
  
  ###################################
  #####           Eq.2          #####
  ###################################
  #Run PQLseq for within_age----------------------------------------------------
  #Generate model matrix
  eq2_matrix<- model.matrix(~ within.age + mean.age + individual_sex + perc_unique, data = long_data)
  eq2_phenotype<- eq2_matrix[, 2]
  eq2_covariates<- eq2_matrix[, 3:5]
  
  #Run pqlseq model
  w.age_pqlseq2_model<- run_pqlseq(eq2_phenotype, eq2_covariates, "eq2_w_age")
  
  #Save pqlseq model
  saveRDS(w.age_pqlseq2_model, paste("prom", "pqlseq2", "within", "age", SAMP, sep = "_"))
  
  #Run PQLseq for mean_age------------------------------------------------------
  #Generate model matrix
  eq2_m_phenotype<- eq2_matrix[, 3]
  eq2_m_covariates<- eq2_matrix[, c(2,4:5)]
  
  #Run pqlseq model
  m.age_pqlseq2_model<- run_pqlseq(eq2_m_phenotype, eq2_m_covariates, "eq2_m_age")
  
  #Save pqlseq model
  saveRDS(m.age_pqlseq2_model, paste("prom", "pqlseq2", "mean", "age", SAMP, sep = "_"))
  
  ###################################
  #####           Eq.3          #####
  ###################################
  #Run PQLseq for eq3-----------------------------------------------------------
  #Generate model matrix
  eq3_matrix<- model.matrix(~ age_at_sampling + mean.age + individual_sex + perc_unique, data = long_data)
  eq3_phenotype<- eq3_matrix[, 2]
  eq3_covariates<- eq3_matrix[, 3:5]
  
  #Run pqlseq model
  eq3_pqlseq2_model<-  run_pqlseq(eq3_phenotype, eq3_covariates, "eq3_age")
  
  #Save pqlseq model
  saveRDS(eq3_pqlseq2_model, paste("prom", "pqlseq2", "eq3", SAMP, sep = "_"))
  
  #Run PQLseq for chronological age---------------------------------------------
  eq3_m_phenotype<- eq3_matrix[, 3]
  eq3_m_covariates<- eq3_matrix[, c(2,4:5)]
  
  #Run pqlseq model
  eq3_m_pqlseq2_model<- run_pqlseq(eq3_m_phenotype, eq3_m_covariates, "eq3_age_m")
  
  #Save pqlseq model
  saveRDS(eq3_m_pqlseq2_model, paste("prom", "pqlseq2", "eq3", "m", SAMP, sep = "_"))
  
} else {
  
  print("long_data lids did not match cov matrix lids")
  
}




