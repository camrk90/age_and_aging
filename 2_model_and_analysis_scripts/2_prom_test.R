#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mem=50G 
#SBATCH --array=1-21

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(tidyverse)
library(PQLseq2)
setwd("/home/ckelsey4/age_and_aging/models_out/")

#Generate function--------------------------------------------------------------
run_pqlseq<- function(pheno, covariates){
  
  mod_df<- pqlseq2(Y = meth, x = pheno, 
                   K = kinship, W = covariates, 
                   lib_size = cov, model="BMM")
  
  mod_df<- mod_df %>%
    filter(converged == TRUE) %>%
    mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
    relocate(fdr, .after = pvalue) %>%
    dplyr::select(-c(converged, elapsed_time))
  
  return(mod_df)
  
}

#Import metadata----------------------------------------------------------------
long_data<- read.table("/scratch/ckelsey4/Cayo_meth/long_data_adjusted.txt")

#Load promoters-----------------------------------------------------------------
prom_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/prom_cov_filtered2.rds")
prom_m<- readRDS("/scratch/ckelsey4/Cayo_meth/prom_m_filtered2.rds")

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

#Import kinship matrix----------------------------------------------------------
kinship<- readRDS("/scratch/ckelsey4/Cayo_meth/full_kin_matrix")

#Subset and rearrange kinship rows and cols to match metadata
kinship<- kinship[long_data$lid_pid, long_data$lid_pid]

#Check metadata lids match the lids (cols) of a random chromosome
if (all.equal(long_data$lid_pid, colnames(prom_cov[[runif(1, 1, 21)]]))) {
  
  #Model Vectors for lme4-------------------------------------------------------
  cov<- prom_cov[[SAMP]]
  meth<- prom_m[[SAMP]]
  
  ###################################
  #####       Eq.3 No Uni       #####
  ###################################
  #Generate model matrix
  eq3_matrix<- model.matrix(~ age_at_sampling + mean.age + individual_sex + perc_unique, data = long_data)
  
  eq3_phenotype<- eq3_matrix[, "age_at_sampling"]
  eq3_covariates<- eq3_matrix[, setdiff(colnames(eq3_matrix), "age_at_sampling")]
  
  eq3_model<- run_pqlseq(eq3_phenotype,eq3_covariates)
  
  #Save pqlseq model
  saveRDS(eq3_model, paste("eq3_no_uni", SAMP, sep = "_"))
  
  rm(eq3_matrix);rm(eq3_model);rm(eq3_phenotype);rm(eq3_covariates)
  
  ###################################
  #####         Eq.3 Uni        #####
  ###################################
  #Generate model matrix
  eq3_matrix<- model.matrix(~ age_at_sampling + mean.age + individual_sex + university, data = long_data)
  
  eq3_phenotype<- eq3_matrix[, "age_at_sampling"]
  eq3_covariates<- eq3_matrix[, setdiff(colnames(eq3_matrix), "age_at_sampling")]
  
  eq3_model<- run_pqlseq(eq3_phenotype,eq3_covariates)
  
  #Save pqlseq model
  saveRDS(eq3_model, paste("eq3_uni", SAMP, sep = "_"))
  
} else {
  
  print("long_data lids did not match cov matrix lids")
  
}




