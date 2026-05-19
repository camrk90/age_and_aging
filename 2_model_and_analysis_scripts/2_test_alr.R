#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mem=50G 
#SBATCH --array=1-21

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(tidyverse)
library(PQLseq2)
setwd("/home/ckelsey4/age_and_aging/models_out/ALR")

#Generate function--------------------------------------------------------------
run_pqlseq<- function(pheno, covariates, type){
  
  mod<- pqlseq2(Y = meth, x = pheno, 
                   K = kinship, W = covariates, 
                   lib_size = cov, model="BMM",
                   verbose = T)
  
  mod_df<- mod@estimates
  
  aics <- sapply(rownames(cov), function(i) {
    mod@others[[i]][["AIC"]]
  })
  
  #mod_df<- cbind(mod_df, aics)
  
  mod_df<- mod_df %>%
    filter(converged == TRUE) %>%
    mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
    relocate(fdr, .after = pvalue) %>%
    dplyr::select(-c(converged, elapsed_time))
  
  colnames(mod_df)<- c("outcome", "n", paste(names(mod_df[,3:length(mod_df)]), type, sep = "_"))
  
  return(list(mod_df=mod_df, aics=aics))
  
}

#Import metadata----------------------------------------------------------------
long_data<- read.table("/scratch/ckelsey4/Cayo_meth/long_data_adjusted.txt")

long_data<- long_data %>%
  group_by(monkey_id) %>%
  mutate(n = n(),
         alr = max(age_at_sampling)) %>%
  ungroup() %>% 
  relocate(alr, .after = age_at_sampling)

#Import kinship matrix----------------------------------------------------------
kinship<- readRDS("/scratch/ckelsey4/Cayo_meth/full_kin_matrix")

#Subset and rearrange kinship rows and cols to match metadata
kinship<- kinship[long_data$lid_pid, long_data$lid_pid]

#Import m/cov rds------------------------------------------------------------
# load region lists that have been filtered for 5x coverage in 90% of samples
regions_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_filtered2.rds")
regions_m<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_m_filtered2.rds")

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
  
#Subset M and Cov to chromosome-------------------------------------------------
cov<- regions_cov[[SAMP]][1:10,]
meth<- regions_m[[SAMP]][1:10,]

###################################
#####        Run PQLseq       #####
###################################
#Run PQLseq for within_age------------------------------------------------------
#Generate model matrix
predictor_matrix<- model.matrix(~ age_at_sampling + alr + individual_sex + university, data = long_data)
phenotype<- predictor_matrix[, 2]
covariates<- predictor_matrix[, 3:5]

#Run pqlseq model
alr_model<- run_pqlseq(phenotype, covariates, "alr")

#Save pqlseq models
saveRDS(alr_model, paste("alr", 'pqlseq', "model", SAMP, sep = "_"))

  
  