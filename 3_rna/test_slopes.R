#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mem=30G 
#SBATCH --array=1-4

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(tidyverse)
library(lme4)
library(limma)
library(edgeR)

#Load data----------------------------------------------------------------------
base_meta<- read.table("/home/ckelsey4/rna_data/base_meta.txt")
rna_counts<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_counts_9Jan26.rds")

base_meta<- base_meta %>%
  arrange(Sample_ID) %>%
  mutate(y = 1)

rna_counts<- rna_counts[, base_meta$Sample_ID]

rna_counts<- t(rna_counts)

age<- base_meta$trapped_age
age.w<- base_meta$within_age
age.m<- base_meta$mean_age
ids<- base_meta$animal_ID
sex<- base_meta$sex
p_gene_counts<- base_meta$p_gene_counts

vars<- c("eq2_int", "eq3_int", "eq2_slope", "eq3_slope")
type<- vars[SAMP]

#Intercepts---------------------------------------------------------------------
run_glm<- function(rna, eq) {
  
  form <- switch(eq, 
                 "eq2_int"  = as.formula("rna[,i] ~ age.w + age.m + sex + p_gene_counts + (1|ids)"), 
                 "eq3_int"  = as.formula("rna[,i] ~ age + age.m + sex + p_gene_counts + (1|ids)"),
                 "eq2_slope"  = as.formula("rna[,i] ~ age.w + age.m + sex + p_gene_counts + (1 + age.w|ids)"), 
                 "eq3_slope"  = as.formula("rna[,i] ~ age + age.m + sex + p_gene_counts + (1 + age|ids)"),
                 stop("Invalid eq value"))
  
  results_list <- vector("list", ncol(rna))
  
  #Run intercept glm
  for (i in 1:ncol(rna)) {
    
    res <- tryCatch({
      
      rfx <- glmer(form, family = poisson(link = "log"))
      
      rfx_sum <- summary(rfx)
      cfs <- as.data.frame(rfx_sum[["coefficients"]])
      cfs$term <- rownames(cfs)
      colnames(cfs) <- c("estimate", "se", "z", "pval", "term")
      
      rfx_row <- pivot_wider(cfs,names_from = term,
                             values_from = c(estimate, se, z, pval))
      
      if(length(rfx@optinfo[["conv"]][["lme4"]]) == 2){
        rfx_row$converged<- "FALSE"
      }else{
        rfx_row$converged<- "TRUE"
      }
      rfx_row<- as.data.frame(rfx_row)
      rfx_row
      
    }, error = function(e) {
      
      message(paste("Gene", i, "failed:", e$message))
      
      # Create NA row
      na_row <- as.data.frame(matrix(NA, nrow = 1, ncol = 21))
      colnames(na_row) <- colnames(rfx_row)
      na_row$converged <- e$message
      na_row
      
    })
    
    results_list[[i]]<- res
    
    print(paste(i, "of", ncol(rna), "done"))
    
  }
  
  # bind once at the end
  df<- do.call(rbind, results_list)
  
  df$gene<- colnames(rna)
  
  return(df)
}

out<- run_glm(rna_counts, type)

saveRDS(out, file.path("/home/ckelsey4/rna_data", paste(type, "_glmer", sep = "")))

  



