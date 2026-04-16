#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mem=50G 
#SBATCH --array=1-21

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(tidyverse)
library(lme4)
library(broom)
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

#Model Vectors for lme4---------------------------------------------------------
SAMP<- 1
cov<- t(regions_cov[[SAMP]])[, 1:50]
meth<- t(regions_m[[SAMP]])[, 1:50]
age<- long_data$age_at_sampling
age.w<- long_data$within.age
age.m<- long_data$mean.age
ids<- long_data$monkey_id
sex<- long_data$individual_sex
p_unique<- long_data$perc_unique

vars<- c("eq2_int", "eq3_int", "eq2_slope", "eq3_slope")

#Intercepts---------------------------------------------------------------------
run_glm<- function(eq) {
  
  form <- switch(eq, 
                 "eq2_int"  = as.formula("cbind(meth[,i],cov[,i]) ~ age.w + age.m + sex + (1|ids)"), 
                 "eq3_int"  = as.formula("meth[,i]/cov[,i] ~ age + age.m + sex + (1|ids)"),
                 "eq2_slope"  = as.formula("cbind(meth[,i],cov[,i]) ~ age.w + age.m + sex + (1 + age.w|ids)"), 
                 "eq3_slope"  = as.formula("meth[,i]/cov[,i] ~ age + age.m + sex + (1 + age|ids)"),
                 stop("Invalid eq value"))
  
  results_list <- vector("list", ncol(meth))
  
  #Run intercept glm
  for (i in 1:ncol(meth)) {
    
    res <- tryCatch({
      
      rfx <- glmer(form, family = binomial(link = "logit"))
      
      rfx_sum <- summary(rfx)
      cfs <- as.data.frame(rfx_sum[["coefficients"]])
      cfs$term <- rownames(cfs)
      colnames(cfs) <- c("estimate", "se", "z", "pval", "term")
      
      rfx_row <- pivot_wider(cfs,names_from = term,
                             values_from = c(estimate, se, z, pval))
      
      if(length(rfx@optinfo[["conv"]][["lme4"]]) != 0){
        rfx_row$note<- rfx@optinfo[["conv"]][["lme4"]][["messages"]][[1]]
      }else{
        rfx_row$note<- "Converged"
      }
      
      rfx_row<- as.data.frame(rfx_row)
      rfx_row
      
    }, error = function(e) {
      
      message(paste("Region", i, "failed:", e$message))
      
      # Create NA row
      na_row <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(rfx_row)))
      colnames(na_row) <- colnames(rfx_row)
      na_row$note<- e$message
      na_row
      
    }, warning = function(w) {
      
      message(w$message)
      
      rfx_sum <- summary(rfx)
      cfs <- as.data.frame(rfx_sum[["coefficients"]])
      cfs$term <- rownames(cfs)
      colnames(cfs) <- c("estimate", "se", "z", "pval", "term")
      
      rfx_row <- pivot_wider(cfs,names_from = term,
                             values_from = c(estimate, se, z, pval))
     
      rfx_row$note<- w$message
      rfx_row<- as.data.frame(rfx_row)
      rfx_row
      
    })
    
    results_list[[i]]<- res
    
    #print(paste(i, "of", ncol(rna), "done"))
    
  }
  
  # bind once at the end
  df<- do.call(rbind, results_list)
  
  df$region<- colnames(meth)
  
  return(results_list)
}

out<- run_glm("eq2_int")
out2<- run_glm("eq2_slope")

saveRDS(out, file.path("/home/ckelsey4/rna_data", paste(type, "_glmer", sep = "")))



