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
library(variancePartition)

#Load data
base_meta<- read.table("/home/ckelsey4/rna_data/base_meta.txt")
rna_counts<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_counts_9Jan26.rds")
rna_kin<- readRDS("/home/ckelsey4/rna_data/rna_kin_matrix.rds")

base_meta<- base_meta %>%
  arrange(Sample_ID) %>%
  mutate(y = 1)

#Normalize RNA Count Data-------------------------------------------------------
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

rna_norm<- rna_norm[rownames(rna_norm) %in% base_meta$Sample_ID,]

#Run EMMA for EQ3---------------------------------------------------------------
run_emma<- function(meta, model){
  
  if (model == "eq1") {
    
    # Create model matrix
    mat<- model.matrix(~ trapped_age + sex + p_gene_counts + p_uniq_mapped, data = meta)
    re_eq<- "y ~ trapped_age + sex + p_gene_counts + p_uniq_mapped + (1|animal_ID)"
    vp_form<- ~ trapped_age + (1|sex) + p_gene_counts + p_uniq_mapped
    
  } else if (model == "eq2") {
    
    # Create model matrix
    mat<- model.matrix(~ within_age + mean_age + sex + p_gene_counts + p_uniq_mapped, data = meta)
    re_eq<- "y ~ within_age + mean_age + sex + p_gene_counts + p_uniq_mapped + (1|animal_ID)"
    vp_form<- ~ within_age + mean_age + (1|sex) + p_gene_counts + p_uniq_mapped
    
  } else if (model == "eq3") {
    
    # Create model matrix
    mat<- model.matrix(~ trapped_age + mean_age + sex + p_gene_counts + p_uniq_mapped, data = meta)
    re_eq<- "y ~ trapped_age + mean_age + sex + p_gene_counts + p_uniq_mapped + (1|animal_ID)"
    vp_form<- ~ trapped_age + mean_age + (1|sex) + p_gene_counts + p_uniq_mapped
    
  }
  
  #Generates random effects matrix
  re_mat <- lFormula(eval(re_eq), meta)
  re_matZ <- as.matrix(t(re_mat$reTrms$Zt))
  
  rna_kin<- rna_kin[colnames(re_matZ), colnames(re_matZ)]
  
  #Generate empty df to input emma output
  df<- data.frame(matrix(nrow=0, ncol=4*(ncol(mat))))
  
  #Run EMMA
  for (i in 1:ncol(rna_norm)) {
    
    em<- emmreml(rna_norm[, i], mat, re_matZ, rna_kin, varbetahat=T,varuhat=T, PEVuhat=T, test=T)
    
    vars<- rownames(em$betahat)
    beta_vals<- as.numeric(em$betahat)
    chi_sq_vals<- as.numeric(em$Xsqtestbeta)
    pvals<- as.numeric(em$pvalbeta[, "fdr"])
    se<- as.numeric(em$varbetahat)
    
    # return a one-row data.frame
    em_row<- data.frame(outcome = colnames(rna_norm)[i],
                        t(setNames(beta_vals, paste0("beta_", vars))),
                        t(setNames(chi_sq_vals, paste0("chi_square_", vars))),
                        t(setNames(pvals, paste0("pvalue_", vars))),
                        t(setNames(se, paste0("se_", vars))),
                        Vu = em$Vu, Ve = em$Ve, loglik = em$loglik,
                        stringsAsFactors = FALSE)
    
    df<- rbind(df, em_row)
    
    print(paste("gene", i, "out of", ncol(rna_norm), "done", sep = " "))
    
  }
  
  #Run variance partition
  vp<- fitExtractVarPartModel(t(rna_norm), vp_form, meta)
  vp<- sortCols(vp)
  
  return(list(df=df, vp=vp))
}

params<- c("eq1", "eq2", "eq3")

model_name<- params[SAMP]

result<- run_emma(meta = base_meta, model = model_name)

saveRDS(result, paste("/home/ckelsey4/rna_data/rna", model_name, sep = "_"))
