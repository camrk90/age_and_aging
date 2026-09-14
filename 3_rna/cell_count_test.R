#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mem=20G

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

cell_counts<- cell_counts %>% 
  dplyr::rename(animal_ID = monkey_id)
cell_counts<- cell_counts %>% 
  dplyr::rename(trapping_ID = trapping_id)

base_meta<- inner_join(base_meta, cell_counts, by = c("animal_ID", "trapping_ID"))

base_meta<- base_meta %>%
  distinct(trapping_ID, .keep_all = T)

base_meta<- base_meta %>%
  group_by(animal_ID) %>%
  mutate(mean_age = mean(trapped_age)) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  arrange(trapping_ID)

library(ggplot2)

base_meta %>%
  distinct(animal_ID, .keep_all = T) %>%
  ggplot(aes(n)) + 
  geom_bar()

hist(base_meta$trapped_age)

#Run EMMA for EQ3---------------------------------------------------------------
run_emma<- function(df, cell_counts){

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
  
  if (cell_counts == T) {
    
    # Create model matrix
    mat <- model.matrix(~ trapped_age + sex + Seq_batch + p_gene_counts + p_uniq_mapped +
                          cd3_cd8_proportion + cd3_cd16_proportion + cd20_proportion, data = df)
    re_eq <- "y ~ trapped_age + mean_age + sex + Seq_batch + p_gene_counts + p_uniq_mapped + cd3_cd4_proportion + 
                          cd3_cd8_proportion + cd3_cd16_proportion + cd20_proportion + (1|animal_ID)"
    
    print(colnames(mat))
    
  } else {
    
    # Create model matrix
    mat <- model.matrix(~ trapped_age + sex + Seq_batch + p_gene_counts + p_uniq_mapped, data = df)
    re_eq <- "y ~ trapped_age + mean_age + sex + Seq_batch + p_gene_counts + p_uniq_mapped + (1|animal_ID)"
    
    print(colnames(mat))
    
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

rna_cell_counts<- run_emma(base_meta, cell_counts = T)
saveRDS(rna_cell_counts, "rna_cell_counts.rds")

rna_no_cell_counts<- run_emma(base_meta, cell_counts = F)
saveRDS(rna_no_cell_counts, "rna_no_cell_counts.rds")

#Rune LME4----------------------------------------------------------------------
run_lme4<- function(df, cell_counts){

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
  
  if (cell_counts == T) {
    
    for (i in 1:ncol(rna_norm)) {
      
      df2<- cbind(df, gene = rna_norm[, i])
      
      res<- tryCatch({
        
        rfx<- lmerTest::lmer(gene ~ trapped_age + mean_age + sex + p_gene_counts + cd3_cd4_proportion + 
                             cd3_cd8_proportion + cd3_cd16_proportion + cd20_proportion + (1|animal_ID), 
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
      
      #print(paste(i, "of", ncol(rna_norm), "done"))
      
    }
    
  } else {
    
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
      
      #print(paste(i, "of", ncol(rna_norm), "done"))
    
    }
    
  }
  
  # bind once at the end
  dd<- bind_rows(results_list)
  
  dd$region<- colnames(rna_norm)
  
  dd<- dd %>%
    dplyr::relocate(region, .before = "estimate_(Intercept)")
  
  return(dd)
}

lm_cells<- suppressMessages(run_lme4(base_meta, cell_counts = T))
saveRDS(lm_cells, "/home/ckelsey4/age_and_aging/models_out/rna_lm_cell_counts.rds")

lm_no_cells<- suppressMessages(run_lme4(base_meta, cell_counts = F))
saveRDS(lm_no_cells, "/home/ckelsey4/age_and_aging/models_out/rna_lm_no_cell_counts.rds")



