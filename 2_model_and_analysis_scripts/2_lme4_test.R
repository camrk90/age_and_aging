#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mem=20G
#SBATCH --array=1-2

SAMP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SAMP <- as.integer(SAMP)

library(tidyverse)
library(lme4)
library(limma)
library(edgeR)
library(EMMREML)

#Load data----------------------------------------------------------------------
base_meta<- read.table("/home/ckelsey4/rna_data/base_meta.txt")
rna_counts<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_counts_9Jan26.rds")
rna_kin<- readRDS("/home/ckelsey4/rna_data/rna_kin_matrix.rds")

base_meta<- base_meta %>%
  arrange(Sample_ID) %>%
  mutate(y = 1)

#Generate model functions-------------------------------------------------------
#emmaEQ3
run_emma<- function(df, model){
  
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
  
  if (model == "eq1") {
    
    # Create model matrix
    mat <- model.matrix(~ trapped_age + sex + p_gene_counts, data = df)
    re_eq <- "y ~ trapped_age + sex + p_gene_counts + (1|animal_ID)"
    
  } else if (model == "eq3") {
    
    # Create model matrix
    mat <- model.matrix(~ trapped_age + mean_age + sex + p_gene_counts, data = df)
    re_eq <- "y ~ trapped_age + mean_age + sex + p_gene_counts + (1|animal_ID)"
    
  }
  
  #Generates random effects matrix
  re_mat <- lFormula(eval(re_eq), df)
  re_matZ <- as.matrix(t(re_mat$reTrms$Zt))
  
  rna_kin<- rna_kin[colnames(re_matZ), colnames(re_matZ)]
  
  df<- data.frame(matrix(nrow=0, ncol=4*(ncol(mat))))
  
  for (i in 1:ncol(rna_norm)) {
    
    em<- emmreml(rna_norm[, i], mat, re_matZ, rna_kin, varbetahat=T,varuhat=T, PEVuhat=T, test=T)
    
    vars<- rownames(em$betahat)[1:4]
    beta_vals<- as.numeric(em$betahat[1:4])
    chi_sq_vals<- as.numeric(em$Xsqtestbeta[1:4])
    pvals<- as.numeric(em$pvalbeta[1:4, "fdr"])
    se<- as.numeric(em$varbetahat[1:4])
    
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
  return(df)
}

#lme4
run_lme4<- function(df, model){
  
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
  
  rna_norm<- as.data.frame(rna_norm[rownames(rna_norm) %in% df$Sample_ID,])
  
  results_list <- vector("list", ncol(rna_norm))
  
  if (model == "eq1") {
    
    for (i in 1:ncol(rna_norm)) {
      
      df2<- cbind(df, gene = rna_norm[, i])
      
      res<- tryCatch({
        
        rfx<- lmerTest::lmer(gene ~ trapped_age + sex + p_gene_counts + (1|animal_ID), 
                             data = df2)
        
        rfx_sum<- summary(rfx)
        cfs<- as.data.frame(rfx_sum[["coefficients"]])
        cfs$term<- rownames(cfs)
        colnames(cfs)<- c("estimate", "se", "df", "tval", "pval", "term")
        
        rfx_row<- pivot_wider(cfs,names_from = term,
                              values_from = c(estimate, se, df, tval, pval, term))
        
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
        
        #message(w$message)
        
        rfx_sum<- summary(rfx)
        cfs<- as.data.frame(rfx_sum[["coefficients"]])
        cfs$term<- rownames(cfs)
        colnames(cfs)<- c("estimate", "se", "df", "tval", "pval", "term")
        
        rfx_row<- pivot_wider(cfs,names_from = term,
                              values_from = c(estimate, se, df, tval, pval, term))
        
        rfx_row$note<- w$message
        rfx_row<- as.data.frame(rfx_row)
        rfx_row
        
      })
      
      results_list[[i]]<- res
      
    }
    
  } else if (model == "eq3") {
    
    for (i in 1:ncol(rna_norm)) {
      
      df2<- cbind(df, gene = rna_norm[, i])
      
      res<- tryCatch({
        
        rfx<- lmerTest::lmer(gene ~ trapped_age + mean_age + sex + p_gene_counts + (1|animal_ID), 
                             data = df2)
        
        rfx_sum<- summary(rfx)
        cfs<- as.data.frame(rfx_sum[["coefficients"]])
        cfs$term<- rownames(cfs)
        colnames(cfs)<- c("estimate", "se", "df", "tval", "pval", "term")
        
        rfx_row<- pivot_wider(cfs,names_from = term,
                              values_from = c(estimate, se, df, tval, pval, term))
        
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
        
        #message(w$message)
        
        rfx_sum<- summary(rfx)
        cfs<- as.data.frame(rfx_sum[["coefficients"]])
        cfs$term<- rownames(cfs)
        colnames(cfs)<- c("estimate", "se", "df", "tval", "pval", "term")
        
        rfx_row<- pivot_wider(cfs,names_from = term,
                              values_from = c(estimate, se, df, tval, pval, term))
        
        rfx_row$note<- w$message
        rfx_row<- as.data.frame(rfx_row)
        rfx_row
        
      })
      
      results_list[[i]]<- res
      
    }
    
  }
  
  # bind once at the end
  dd<- bind_rows(results_list)
  
  dd$region<- colnames(rna_norm)
  
  dd<- dd %>%
    dplyr::relocate(region, .before = "estimate_(Intercept)")
  
  return(dd)
}

#Run models---------------------------------------------------------------------
params<- c("eq1", "eq3")

model_name<- params[SAMP]

result_emma<- run_emma(df = base_meta, model = model_name)
saveRDS(result_emma, paste("/home/ckelsey4/rna_data/rna_emma", model_name, sep = "_"))

result_lme4<- run_lme4(df = base_meta, model = model_name)
saveRDS(result_lme4, paste("/home/ckelsey4/rna_data/rna_lme4", model_name, sep = "_"))




