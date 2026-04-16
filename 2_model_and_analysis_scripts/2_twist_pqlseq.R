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

#Import methylation data--------------------------------------------------------
cov<- readRDS("/home/ckelsey4/Cayo_meth/twist_cov_list")
meth<- readRDS("/home/ckelsey4/Cayo_meth/twist_meth_list")

#Import metadata----------------------------------------------------------------
long_data<- read.table("/home/ckelsey4/Cayo_meth/twist_metadata.txt", sep = "\t", header=T)

#Import kinship matrix----------------------------------------------------------
kinship<- readRDS("/home/ckelsey4/Cayo_meth/dnam_kin_matrix.rds")

#Subset and rearrange kinship rows and cols to match metadata
kinship<- kinship[long_data$vantage_id, long_data$vantage_id]

#Check metadata lids match the lids (cols) of a random chromosome
if (all.equal(long_data$vantage_id, colnames(cov[[runif(1, 1, 21)]]))) {
  
  #Model Vectors for lme4-------------------------------------------------------
  cov<- cov[[SAMP]]
  meth<- meth[[SAMP]]
  
  ###################################
  #####           Eq.1          #####
  ###################################
  #Run PQLseq for chronological age---------------------------------------------
  #Generate model matrix
  eq1_matrix<- model.matrix(~ age + subject_sex + plate_id, data = long_data)
  eq1_phenotype<- eq1_matrix[, 2]
  eq1_covariates<- as.matrix(eq1_matrix[, 3:10])

  #Run pqlseq model
  agechron_pqlseq2_model<- run_pqlseq(eq1_phenotype, eq1_covariates, "chron_age")

  #Save pqlseq model
  saveRDS(agechron_pqlseq2_model, paste("twist", "pqlseq2", "agechron", SAMP, sep = "_"))
  
  ###################################
  #####           Eq.2          #####
  ###################################
  #Run PQLseq for within_age----------------------------------------------------
  #Generate model matrix
  eq2_matrix<- model.matrix(~ within_age + mean_age + subject_sex + plate_id, data = long_data)
  eq2_phenotype<- eq2_matrix[, 2]
  eq2_covariates<- eq2_matrix[, 3:11]
  
  #Run pqlseq model
  w.age_pqlseq2_model<- run_pqlseq(eq2_phenotype, eq2_covariates, "eq2_w_age")
  
  #Save pqlseq model
  saveRDS(w.age_pqlseq2_model, paste("twist", "pqlseq2", "within", "age", SAMP, sep = "_"))
  
  #Run PQLseq for mean_age------------------------------------------------------
  #Generate model matrix
  eq2_m_phenotype<- eq2_matrix[, 3]
  eq2_m_covariates<- eq2_matrix[, c(2,4:11)]
  
  #Run pqlseq model
  m.age_pqlseq2_model<- run_pqlseq(eq2_m_phenotype, eq2_m_covariates, "eq2_m_age")
  
  #Save pqlseq model
  saveRDS(m.age_pqlseq2_model, paste("twist", "pqlseq2", "mean", "age", SAMP, sep = "_"))
  
  ###################################
  #####           Eq.3          #####
  ###################################
  #Run PQLseq for eq3-----------------------------------------------------------
  #Generate model matrix
  eq3_matrix<- model.matrix(~ age + mean_age + subject_sex + plate_id, data = long_data)
  eq3_phenotype<- eq3_matrix[, 2]
  eq3_covariates<- eq3_matrix[, 3:11]
  
  #Run pqlseq model
  eq3_pqlseq2_model<-  run_pqlseq(eq3_phenotype, eq3_covariates, "eq3_age")
  
  #Save pqlseq model
  saveRDS(eq3_pqlseq2_model, paste("twist", "pqlseq2", "eq3", SAMP, sep = "_"))
  
  #Run PQLseq for chronological age---------------------------------------------
  eq3_m_phenotype<- eq3_matrix[, 3]
  eq3_m_covariates<- eq3_matrix[, c(2,4:11)]
  
  #Run pqlseq model
  eq3_m_pqlseq2_model<- run_pqlseq(eq3_m_phenotype, eq3_m_covariates, "eq3_age_m")
  
  #Save pqlseq model
  saveRDS(eq3_m_pqlseq2_model, paste("twist", "pqlseq2", "eq3", "m", SAMP, sep = "_"))
  
} else {
  
  print("long_data lids did not match cov matrix lids")
  
}



