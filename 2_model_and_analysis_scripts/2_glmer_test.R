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
cov<- t(regions_cov[[SAMP]])
meth<- t(regions_m[[SAMP]])
age<- long_data$age_at_sampling
age.w<- long_data$within.age
age.m<- long_data$mean.age
ids<- long_data$monkey_id
sex<- long_data$individual_sex

###################################
#####   Run Intercept GLMM    #####
###################################
#Run GLMMM----------------------------------------------------------------------
#Generate blank objects for glms
eq2_df_int<- data.frame(matrix(nrow=0, ncol=21))
colnames(eq2_df_int)<- c("estimate_(Intercept)", "estimate_age.w", "estimate_age.m", "estimate_sexM",
                 "se_(Intercept)", "se_age.w", "se_age.m", "se_sexM", 
                 "z_(Intercept)", "z_age.w", "z_age.m", "z_sexM",
                 "pval_(Intercept)", "pval_age.w", "pval_age.m", "pval_sexM",
                 "aic", "bic", "loglik", "converged")

#Run slope/intercept glm
for(i in 1:ncol(meth)){
  
  rfx=glmer(cbind(meth[, i], cov[, i]) ~ age.w + age.m + sex + (1|ids), 
            family = binomial(link = "logit"))
  
  slopes<- as.data.frame(ranef(rfx))
  rfx_sum<- summary(rfx)
  cfs<- rfx_sum[["coefficients"]]
  cfs<- as.data.frame(cfs)
  cfs$term<- rownames(cfs)
  #cfs<- cfs[1:4, 1:5] #this is for when the batch variable gets incorporated
  colnames(cfs)<- c("estimate", "se", "z", "pval", "term")
  rfx_row<- pivot_wider(cfs, names_from = term, values_from = c(estimate, se, z, pval))
  rfx_row$aic<- rfx_sum[["AICtab"]][["AIC"]]
  rfx_row$bic<- rfx_sum[["AICtab"]][["BIC"]]
  rfx_row$loglik<- rfx_sum[["logLik"]]
  
  if(length(rfx@optinfo[["conv"]][["lme4"]]) == 2){
    rfx_row$converged<- "FALSE"
  }else{
    rfx_row$converged<- "TRUE"
  }
  
  eq2_df_int<- rbind(eq2_df_int, rfx_row)

}

eq2_df_int$region<- colnames(meth)

#Save outputs to RDS
saveRDS(eq2_df_int, paste("glmer", "eq2", "int", SAMP, sep = "_"))

#EQ3----------------------------------------------------------------------------
#Generate blank objects for glms
eq3_df_int<- data.frame(matrix(nrow=0, ncol=21))
colnames(eq3_df_int)<- c("estimate_(Intercept)", "estimate_age", "estimate_age.m", "estimate_sexM",
                     "se_(Intercept)", "se_age", "se_age.m", "se_sexM", 
                     "z_(Intercept)", "z_age", "z_age.m", "z_sexM",
                     "pval_(Intercept)", "pval_age", "pval_age.m", "pval_sexM",
                     "aic", "bic", "loglik", "converged")

#Run slope/intercept glm
for(i in 1:ncol(meth)){
  
  rfx=glmer(cbind(meth[, i], cov[, i]) ~ age + age.m + sex + (1|ids), 
            family = binomial(link = "logit"))
  
  slopes<- as.data.frame(ranef(rfx))
  rfx_sum<- summary(rfx)
  cfs<- rfx_sum[["coefficients"]]
  cfs<- as.data.frame(cfs)
  cfs$term<- rownames(cfs)
  #cfs<- cfs[1:4, 1:5] #this is for when the batch variable gets incorporated
  colnames(cfs)<- c("estimate", "se", "z", "pval", "term")
  rfx_row<- pivot_wider(cfs, names_from = term, values_from = c(estimate, se, z, pval))
  rfx_row$aic<- rfx_sum[["AICtab"]][["AIC"]]
  rfx_row$bic<- rfx_sum[["AICtab"]][["BIC"]]
  rfx_row$loglik<- rfx_sum[["logLik"]]
  
  if(length(rfx@optinfo[["conv"]][["lme4"]]) == 2){
    rfx_row$converged<- "FALSE"
  }else{
    rfx_row$converged<- "TRUE"
  }
  
  eq3_df_int<- rbind(eq3_df_int, rfx_row)

}

#Save outputs to RDS
saveRDS(eq3_df_int, paste("glmer", "eq3", "int", SAMP, sep = "_"))

###################################
#####     Run Slope GLMM      #####
###################################
#Run GLMMM----------------------------------------------------------------------
#Generate blank objects for glms
eq2_df_slopes<- data.frame(matrix(nrow=0, ncol=21))
colnames(eq2_df_slopes)<- c("estimate_(Intercept)", "estimate_age.w", "estimate_age.m", "estimate_sexM",
                     "se_(Intercept)", "se_age.w", "se_age.m", "se_sexM", 
                     "z_(Intercept)", "z_age.w", "z_age.m", "z_sexM",
                     "pval_(Intercept)", "pval_age.w", "pval_age.m", "pval_sexM",
                     "aic", "bic", "loglik", "converged")

#Run slope/intercept glm
for(i in 1:ncol(meth)){
  
  rfx=glmer(cbind(meth[, i], cov[, i]) ~ age.w + age.m + sex + (1 + age.w|ids), 
            family = binomial(link = "logit"))
  
  slopes<- as.data.frame(ranef(rfx))
  rfx_sum<- summary(rfx)
  cfs<- rfx_sum[["coefficients"]]
  cfs<- as.data.frame(cfs)
  cfs$term<- rownames(cfs)
  #cfs<- cfs[1:4, 1:5] #this is for when the batch variable gets incorporated
  colnames(cfs)<- c("estimate", "se", "z", "pval", "term")
  rfx_row<- pivot_wider(cfs, names_from = term, values_from = c(estimate, se, z, pval))
  rfx_row$aic<- rfx_sum[["AICtab"]][["AIC"]]
  rfx_row$bic<- rfx_sum[["AICtab"]][["BIC"]]
  rfx_row$loglik<- rfx_sum[["logLik"]]
  
  if(length(rfx@optinfo[["conv"]][["lme4"]]) == 2){
    rfx_row$converged<- "FALSE"
  }else{
    rfx_row$converged<- "TRUE"
  }
  
  eq2_df_slopes<- rbind(eq2_df_slopes, rfx_row)
  
}

eq2_df_slopes$region<- colnames(meth)

#Save outputs to RDS
saveRDS(eq2_df_slopes, paste("glmer", "eq2", "slopes", SAMP, sep = "_"))

#EQ3----------------------------------------------------------------------------
#Generate blank objects for glms
eq3_df_slopes<- data.frame(matrix(nrow=0, ncol=21))
colnames(eq3_df_slopes)<- c("estimate_(Intercept)", "estimate_age", "estimate_age.m", "estimate_sexM",
                     "se_(Intercept)", "se_age", "se_age.m", "se_sexM", 
                     "z_(Intercept)", "z_age", "z_age.m", "z_sexM",
                     "pval_(Intercept)", "pval_age", "pval_age.m", "pval_sexM",
                     "aic", "bic", "loglik", "converged")

#Run slope/intercept glm
for(i in 1:ncol(meth)){
  
  rfx=glmer(cbind(meth[, i], cov[, i]) ~ age + age.m + sex + (1 + age|ids), 
            family = binomial(link = "logit"))
  
  slopes<- as.data.frame(ranef(rfx))
  rfx_sum<- summary(rfx)
  cfs<- rfx_sum[["coefficients"]]
  cfs<- as.data.frame(cfs)
  cfs$term<- rownames(cfs)
  #cfs<- cfs[1:4, 1:5] #this is for when the batch variable gets incorporated
  colnames(cfs)<- c("estimate", "se", "z", "pval", "term")
  rfx_row<- pivot_wider(cfs, names_from = term, values_from = c(estimate, se, z, pval))
  rfx_row$aic<- rfx_sum[["AICtab"]][["AIC"]]
  rfx_row$bic<- rfx_sum[["AICtab"]][["BIC"]]
  rfx_row$loglik<- rfx_sum[["logLik"]]
  
  if(length(rfx@optinfo[["conv"]][["lme4"]]) == 2){
    rfx_row$converged<- "FALSE"
  }else{
    rfx_row$converged<- "TRUE"
  }
  
  eq3_df_slopes<- rbind(eq3_df_slopes, rfx_row)
  
}

#Save outputs to RDS
saveRDS(eq3_df_slopes, paste("glmer", "eq3", "slopes", SAMP, sep = "_"))



