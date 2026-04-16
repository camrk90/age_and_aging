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
  #Eq.1-------------------------------------------------------------------------
  #Generate model matrix
  eq1_matrix<- model.matrix(~ age_at_sampling + individual_sex + perc_unique, data = long_data)
    
  vars <- c("age_at_sampling")
  
  eq1_model <- lapply(setNames(vars, vars), function(i) {
    
    eq1_phenotype<- eq1_matrix[, i]
    eq1_covariates<- eq1_matrix[, setdiff(colnames(eq1_matrix), i)]
    
    run_pqlseq(eq1_phenotype, eq1_covariates)
    
  })
  
  #Save pqlseq model
  saveRDS(eq1_model, paste("dnam_prom_eq1_model", SAMP, sep = "_"))
  
  ###################################
  #####           Eq.2          #####
  ###################################
  #Generate model matrix
  eq2_matrix<- model.matrix(~ within.age + mean.age + individual_sex + perc_unique, data = long_data)
  
  vars <- c("within.age", "mean.age")
  
  eq2_model <- lapply(setNames(vars, vars), function(i) {
    
    eq2_phenotype <- eq2_matrix[, i]
    eq2_covariates <- eq2_matrix[, setdiff(colnames(eq2_matrix), i)]
    
    run_pqlseq(eq2_phenotype, eq2_covariates)
    
  })
  
  #Save pqlseq model
  saveRDS(eq2_model, paste("dnam_prom_eq2_model", SAMP, sep = "_"))
  
  ###################################
  #####           Eq.3          #####
  ###################################
  #Run PQLseq for eq3-----------------------------------------------------------
  #Generate model matrix
 eq3_matrix<- model.matrix(~ age_at_sampling + mean.age + individual_sex + perc_unique, data = long_data)
  
  vars <- c("age_at_sampling", "mean.age")
  
 eq3_model <- lapply(setNames(vars, vars), function(i) {
    
   eq3_phenotype<- eq3_matrix[, i]
   eq3_covariates<- eq3_matrix[, setdiff(colnames(eq3_matrix), i)]
    
    run_pqlseq(eq3_phenotype,eq3_covariates)
    
  })
  
  #Save pqlseq model
  saveRDS(eq3_model, paste("dnam_prom_eq3_model", SAMP, sep = "_"))
  
} else {
  
  print("long_data lids did not match cov matrix lids")
  
}




