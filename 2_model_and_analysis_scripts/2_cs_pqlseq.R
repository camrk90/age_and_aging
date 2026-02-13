#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mem=10G 
#SBATCH --array=1-21

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

#### THIS SCRIPT RUNS PQLseq FOR THE CROSS SECTIONAL DATA ####

library(tidyverse)
library(PQLseq2)
setwd("/scratch/ckelsey4/Cayo_meth/cross_models")

#Define model function----
run_cs_model<- function(meta){
  
  #Subset and rearrange kinship rows and cols to match metadata
  kin<- kin[meta$lid_pid, meta$lid_pid]
  
  #Check kinship matrix dims match metadata
  all.equal(nrow(meta), ncol(kin))
  all.equal(nrow(meta), nrow(kin))
  
  #Filter metadata to lids in regions list
  meta<- meta[meta$lid_pid %in% colnames(cov),]
  
  cov<- subset(cov, select=meta$lid_pid)
  meth<- subset(meth, select=meta$lid_pid)
  
  #Check metadata lids match the lids (cols) of a random chromosome
  if (all.equal(meta$lid_pid, colnames(cov)) == T) {
    
    ###################################
    #####        Run PQLseq       #####
    ###################################
    #Run PQLseq for within_age-------------------------------------------------------------
    #Generate model matrix
    predictor_matrix<- model.matrix(~ age_at_sampling + individual_sex + unique, data = meta)
    age_pheno<- predictor_matrix[, 2]
    age_cov<- predictor_matrix[, 3:4]
    
    #Run pqlseq model
    age_model<- pqlseq2(Y = meth, x = age_pheno, 
                        K = kin, W = age_cov, 
                        lib_size = cov, model="BMM")
    
  } else {
    
    print("metadata lids did not match cov matrix lids")
    
  }
  
}

#Import metadata----------------------------------------------------------------
blood_metadata<- read.table("/scratch/ckelsey4/Cayo_meth/blood_metadata_full.txt")
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

overlap_lids<- overlap_lids %>%
  filter(age_at_sampling > 1)

rm(lids_to_remove);rm(long_data);rm(long_ids)

#Import kinship matrix
kin<- readRDS("/scratch/ckelsey4/Cayo_meth/full_kin_matrix")

#Import m/cov rds
regions_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_filtered.rds")
cov<- regions_cov[[SAMP]]

regions_m<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_m_filtered.rds")
meth<- regions_m[[SAMP]]

#Run model with ids that 100% overlap longitudinal set
short<- run_cs_model(meta=overlap_lids)

#Save short model
saveRDS(short, paste("cs", SAMP, "short.rds", sep = "_"))

#Run model with all cross-sectional ids
long<- run_cs_model(meta=blood_metadata)

#Save long model
saveRDS(long, paste("cs", SAMP, "long.rds", sep = "_"))

