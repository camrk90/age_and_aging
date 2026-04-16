#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mem=50G 
#SBATCH --array=1-3

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(tidyverse)
library(lme4)
library(limma)
library(edgeR)
library(EMMREML)

#Load data
base_meta<- read.table("/home/ckelsey4/rna_data/base_meta.txt")
rna_counts<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_counts_9Jan26.rds")
rna_kin<- readRDS("/home/ckelsey4/rna_data/rna_kin_matrix.rds")

#Normalize RNA Count Data-------------------------------------------------------
base_meta<- base_meta %>%
  arrange(Sample_ID) %>%
  mutate(y = 1)

rna_counts<- rna_counts[, base_meta$Sample_ID]

#Generate normalized counts
if (all.equal(base_meta$Sample_ID, colnames(rna_counts)) == T) {
  
  rna_norm<- voom(calcNormFactors(DGEList(counts=rna_counts)), plot=FALSE)
  rna_norm<- rna_norm[["E"]]
  colnames(rna_norm)<- colnames(rna_counts)
  rownames(rna_norm)<- rownames(rna_counts)
  rna_norm<- t(rna_norm)
  
} else {
  print("Metadata Sample IDs and RNA Count cols do not match")
}

#Run EMMA-----------------------------------------------------------------------
#Generate emma function
run_emma<- function(eq, re) {
  
  # Define fixed effects based on eq
  fixed_part <- switch(eq,
                       "eq2"  = "within_age + mean_age + sex + p_gene_counts",
                       "eq3"  = "trapped_age + mean_age + sex + p_gene_counts",
                       "chron"= "trapped_age + sex + p_gene_counts",
                       stop("Invalid eq value"))
  
  # Determine slope variable (if needed)
  slope_var <- switch(eq,
                      "eq2"  = "within_age", "eq3"  = "trapped_age",
                      "chron"= "trapped_age")
  
  # Define random effects based on re
  random_part <- if (re == "intercept") {
    "(1|animal_ID)"
  } else if (re == "slopes") {
    paste0("(1 + ", slope_var, "|animal_ID)")
  } else {
    stop("Invalid re value")
  }
  
  # Create model matrix
  mat <- model.matrix(as.formula(paste("~", fixed_part)), data = base_meta)
  
  # Create full model formula
  re_eq <- paste("y ~", fixed_part, "+", random_part)
  
  #Generates random effects matrix
  re_mat <- lFormula(eval(re_eq), base_meta)
  re_matZ <- as.matrix(t(re_mat$reTrms$Zt))
  
  rna_kin<- rna_kin[colnames(re_matZ), colnames(re_matZ)]
  
  df<- data.frame(matrix(nrow=0, ncol=16))
  
  for (i in 1:ncol(rna_norm)) {
    
    em<- emmreml(rna_norm[, i], mat, re_matZ, rna_kin, varbetahat=T,varuhat=T, PEVuhat=T, test=T)
    
    vars <- rownames(em$betahat)[1:4]
    beta_vals   <- as.numeric(em$betahat[1:4, 1])
    chi_sq_vals <- as.numeric(em$Xsqtestbeta[1:4, 1])
    pvals       <- as.numeric(em$pvalbeta[1:4, "none"])
    
    # return a one-row data.frame (same shape as serial version)
    em_row<- data.frame(outcome = colnames(rna_norm)[i],
                        t(setNames(beta_vals, paste0("beta_", vars))),
                        t(setNames(chi_sq_vals, paste0("chi_square_", vars))),
                        t(setNames(pvals, paste0("pvalue_", vars))),
                        Vu = em$Vu, Ve = em$Ve, loglik = em$loglik,
                        stringsAsFactors = FALSE)
    
    df<- rbind(df, em_row)
    
  }
  
  return(df)
  
}

params <- list(
  c("eq2", "intercept"),
  c("eq2", "slopes"),
  c("eq3", "intercept"),
  c("eq3", "slopes"),
  c("chron", "intercept"),
  c("chron", "slopes")
)

model_name<- params[[SAMP]][1]
type_name<- params[[SAMP]][2]

result <- run_emma(model_name, type_name)

saveRDS(result,
        file.path("/home/ckelsey4/rna_data",
                  paste0("rna_", model_name, "_",
                         ifelse(type_name == "intercept", "int", "slopes"))))



