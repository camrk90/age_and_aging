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
regions_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_filtered.rds")
regions_m<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_m_filtered.rds")

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
  
#Separate m/cov into chromosomes based on array number--------------------------
cov<- regions_cov[[SAMP]]
meth<- regions_m[[SAMP]]

###################################
#####        Run PQLseq       #####
###################################
#Run PQLseq-------------------------------------------------------------------
#Generate model matrix
cs_matrix<- model.matrix(~ age_at_sampling + individual_sex + university, data = blood_metadata)

vars <- c("age_at_sampling", "individual_sexM")

cs_model<- lapply(setNames(vars, vars), function(i) {
  
  pheno <- cs_matrix[, i]
  covariates <- cs_matrix[, setdiff(colnames(cs_matrix), i)]
  
  rr_list<- vector("list", nrow(cov))
 # err_list<- vector("list", nrow(cov))
  
  for (r in 1:nrow(cov)) {
    
    tryCatch({
      
      rr<- pqlseq2(Y = meth[20:25,], x = pheno, 
                          K = kinship, W = covariates, 
                          lib_size = cov[20:25,], model="BMM")
      
      rr<- rr %>%
        filter(converged == TRUE) %>%
        mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
        relocate(fdr, .after = pvalue) %>%
        dplyr::select(-c(converged, elapsed_time))
      
      colnames(rr)<- c("outcome", "n", paste(names(rr[,3:length(rr)]), i, sep = "_"))
      
      rr$note <- NA_character_
      
      rr_list[[r]]<- rr
      
    }, error = function(e) {
      
      message(e$message)
      
      err<- as.data.frame(matrix(NA, nrow = 1, ncol = 10))
      
      colnames(err)<- c("outcome", "n", 
                        paste(c("intercept", "se_intercept","beta", "se_beta", "pvalue", "fdr", "h2", "sigma2"), 
                              i,
                        sep = "_"))
      
      err$note<- e$message
      
      rr_list[[r]]<<- err
      
    })
    
    print(paste("Region", r, "done", sep = " "))
    
  } 
  
  # combine into data frames
  rr_df <- dplyr::bind_rows(rr_list)
  
  # return both as a list
  list(results = rr_df)
  
})

#Save pqlseq model
saveRDS(cs_model, paste("cs_pqlseq2_age", SAMP, sep = "_"))
