#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mem=50G 

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(tidyverse)
library(lme4)
library(limma)
library(edgeR)
library(EMMREML)

#Load data
base_meta<- read.table("/home/ckelsey4/rna_data/base_meta.txt")
cell_counts<- readRDS("/scratch/ckelsey4/Cayo_meth/cell_counts/lymphocyte_proportions.rds")
rna_counts<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_counts_9Jan26.rds")
rna_kin<- readRDS("/home/ckelsey4/rna_data/rna_kin_matrix.rds")

#Normalize RNA Count Data-------------------------------------------------------
base_meta<- base_meta %>%
  arrange(Sample_ID) %>%
  mutate(y = 1)

cell_counts<- cell_counts %>% rename(animal_ID = monkey_id)
cell_counts<- cell_counts %>% rename(trapping_ID = trapping_id)

base_meta<- left_join(base_meta, cell_counts, by = c("animal_ID", "trapping_ID"))

base_meta<- base_meta %>%
  drop_na()

base_meta<- base_meta %>%
  distinct(trapping_ID, .keep_all = T)

#Run EMMA for EQ3---------------------------------------------------------------
run_emma<- function(df){
  
  df<- base_meta
  rna_counts<- rna_counts[, df$Sample_ID]
  
  #Generate normalized counts
  if (all.equal(df$Sample_ID, colnames(rna_counts)) == T) {
    
    rna_norm<- voom(calcNormFactors(DGEList(counts=rna_counts)), plot=FALSE)
    rna_norm<- rna_norm[["E"]]
    colnames(rna_norm)<- colnames(rna_counts)
    rownames(rna_norm)<- rownames(rna_counts)
    rna_norm<- t(rna_norm)
    
  } else {
    print("Metadata Sample IDs and RNA Count cols do not match")
  }
  
  rna_norm<- rna_norm[rownames(rna_norm) %in% df$Sample_ID,]
  
  # Create model matrix
  mat <- model.matrix(~ trapped_age + sex + mean_age + p_gene_counts, data = df)
  re_eq <- "y ~ trapped_age + sex + mean_age + p_gene_counts + (1|animal_ID)"
  
  #Generates random effects matrix
  re_mat <- lFormula(eval(re_eq), df)
  re_matZ <- as.matrix(t(re_mat$reTrms$Zt))
  
  rna_kin<- rna_kin[colnames(re_matZ), colnames(re_matZ)]
  
  df<- data.frame(matrix(nrow=0, ncol=4*(ncol(mat))))
  
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
    
    print(paste("gene", i, "out of", ncol(rna_norm), "done", sep = " "))
    
  }
  return(df)
}

rna_no_cell_counts<- run_emma(base_meta)

saveRDS(df, "rna_interaction_out.rds")

#Run EMMA for EQ1---------------------------------------------------------------
run_eq1_emma<- function(df, sex){
  
  if (sex == "M") {
    df<- df %>%
      filter(sex == "M")
  } else if (sex == "F") {
    df<- df %>%
      filter(sex == "F")
  } else if (sex == "both") {
    df=df
  }
  
  rna_counts<- rna_counts[, df$Sample_ID]
  
  #Generate normalized counts
  if (all.equal(df$Sample_ID, colnames(rna_counts)) == T) {
    
    rna_norm<- voom(calcNormFactors(DGEList(counts=rna_counts)), plot=FALSE)
    rna_norm<- rna_norm[["E"]]
    colnames(rna_norm)<- colnames(rna_counts)
    rownames(rna_norm)<- rownames(rna_counts)
    rna_norm<- t(rna_norm)
    
  } else {
    print("Metadata Sample IDs and RNA Count cols do not match")
  }
  
  rna_norm<- rna_norm[rownames(rna_norm) %in% df$Sample_ID,]
  
  # Create model matrix
  if (sex == "M" | sex == "F") {
    mat <- model.matrix(~ trapped_age + p_gene_counts, data = df)
    re_eq <- "y ~ trapped_age + p_gene_counts + (1|animal_ID)"
    
  } else if (sex == "both") {
    mat <- model.matrix(~ trapped_age:sex + p_gene_counts, data = df)
    re_eq <- "y ~ trapped_age:sex + p_gene_counts + (1|animal_ID)"
  }
  
  #Generates random effects matrix
  re_mat <- lFormula(eval(re_eq), df)
  re_matZ <- as.matrix(t(re_mat$reTrms$Zt))
  
  rna_kin<- rna_kin[colnames(re_matZ), colnames(re_matZ)]
  
  df<- data.frame(matrix(nrow=0, ncol=4*(ncol(mat))))
  
  for (i in 1:ncol(rna_norm)) {
    
    em<- emmreml(rna_norm[, i], mat, re_matZ, rna_kin, varbetahat=T,varuhat=T, PEVuhat=T, test=T)
    
    vars <- rownames(em$betahat)[1:3]
    beta_vals   <- as.numeric(em$betahat[1:3, 1])
    chi_sq_vals <- as.numeric(em$Xsqtestbeta[1:3, 1])
    pvals       <- as.numeric(em$pvalbeta[1:3, "none"])
    
    # return a one-row data.frame (same shape as serial version)
    em_row<- data.frame(outcome = colnames(rna_norm)[i],
                        t(setNames(beta_vals, paste0("beta_", vars))),
                        t(setNames(chi_sq_vals, paste0("chi_square_", vars))),
                        t(setNames(pvals, paste0("pvalue_", vars))),
                        Vu = em$Vu, Ve = em$Ve, loglik = em$loglik,
                        stringsAsFactors = FALSE)
    
    df<- rbind(df, em_row)
    
    print(paste("gene", i, "out of", ncol(rna_norm), "done", sep = " "))
    
  }
  return(df)
  
}

rna_eq1_male<- run_eq1_emma(base_meta, "M")
rna_eq1_female<- run_eq1_emma(base_meta, "F")
rna_eq1_nested<- run_eq1_emma(base_meta, "both")

rna_eq1_sex<- list(rna_m = rna_eq1_male, rna_f = rna_eq1_female)

saveRDS(rna_eq1_sex, "rna_eq1_sex.rds")

