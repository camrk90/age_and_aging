library(tidyverse)
library(EMMREML)
library(lme4)
library(ggplot2)
library(variancePartition)
library(limma)
library(edgeR)
#load("/home/ckelsey4/Cayo_meth/rna_compare.RData")

rna_meta_full<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_metadata_wELA_18Jul25.rds")
rna_counts<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_counts_18Jul25.rds")
rna_kin<- readRDS("/scratch/ckelsey4/Cayo_meth/rna_kin_matrix.rds")

#Subset to control and truncate df to relevant variables
rna_meta<- rna_meta_full %>%
  filter(Stimulation == "H2O") %>%
  select(c(animal_ID, sex, trapped_age, trapping_ID, Sample_ID, Seq_batch))

#Generate within and between age
rna_meta<- rna_meta %>%
  group_by(animal_ID) %>%
  mutate(mean_age = mean(trapped_age)) %>%
  mutate(within_age = trapped_age - mean_age) %>%
  arrange(animal_ID, trapped_age)

#Make batch factor
rna_meta$Seq_batch<- gsub("batch", "", rna_meta$Seq_batch)
rna_meta$Seq_batch<- factor(rna_meta$Seq_batch, levels = c("1", "2", "3"))

#Subset counts to RIDs in metadata
rna_counts<- rna_counts[, colnames(rna_counts) %in% rna_meta$Sample_ID]

#Plot rna smaples 
rna_meta<- rna_meta %>%
  group_by(animal_ID) %>%
  mutate(min_age = min(trapped_age))
rna_meta$trapped_age<- round(rna_meta$trapped_age, 0)

rna_meta %>%
  ggplot(aes(x=trapped_age, y=reorder(animal_ID, min_age), colour=sex)) +
  geom_path(linewidth = 1.2, alpha = 0.8) +
  geom_point(colour="black") +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  #scale_colour_manual(values = c("purple", "purple4"), name = "Sex") +
  ylab("Individual") +
  xlab("Age") +
  theme_classic(base_size = 24) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#Variance Partition-------------------------------------------------------------
vp_model<- ~ within_age + mean_age + (1|sex) + (1|Seq_batch)

vp<- fitExtractVarPartModel(t(rna_counts), vp_model, rna_meta)

plotVarPart(vp)
plotPercentBars(vp[1:10, ])

#PCA----------------------------------------------------------------------------
rna_pca<- prcomp(cor(t(rna_counts), use="pairwise.complete.obs"))
summary(rna_pca)

pcs<- as.data.frame(rna_pca$x) 

#Check which pca's explain the most variance
summary(rna_pca)$importance[2, ]

pcs<- cbind(pcs[1:5], 
            dplyr::select(blood_metadata, 
                          c("age_at_sampling", "sex", "pid", "prep_date", "prep_year", 
                            "percent_mapped", "log_unique", "conversion_rate", "university", "global_meth")))

pcs %>% ggplot(aes(x = pid, y = age_at_sampling)) +
  geom_boxplot()

pcs %>% ggplot(aes(x = prep_year, y = age_at_sampling)) +
  geom_boxplot()

blood_metadata$university<- factor(blood_metadata$university, levels = c("uw", "asu"))

pc.matrix<- model.matrix(~ PC1 + PC2 + PC3 + age_at_sampling + sex + university + global_meth, data = pcs)
pc.matrix[, 2:8] %>% 
  cor(use="pairwise.complete.obs") %>%
  ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2)

#Normalize RNA Count Data-------------------------------------------------------
m_matrix <- model.matrix(~ within_age + mean_age + sex + Seq_batch, rna_meta)
rna_norm<- voom(calcNormFactors(DGEList(counts=rna_counts)), m_matrix, plot=FALSE)
rna_norm<- rna_norm[["E"]]
colnames(rna_norm)<- colnames(rna_counts)
rownames(rna_norm)<- rownames(rna_counts)
rna_norm<- t(rna_norm)

#Run EMMA for within_age--------------------------------------------------------
#Generate model matrix
within_matrix<- model.matrix(~ within_age + mean_age + sex + Seq_batch, data = rna_meta)

re_eq = "y ~ within_age + mean_age + sex + Seq_batch + (1|animal_ID)"

rna_meta$y<- 1

#This generates a random intercept matrix with dims 309 x 114
re_mat <- lFormula(eval(re_eq), rna_meta)
re_matZ <- as.matrix(t(re_mat$reTrms$Zt))

rna_kin<- rna_kin[colnames(re_matZ), colnames(re_matZ)]

df<- data.frame(matrix(nrow=0, ncol=16))

for (i in 1:ncol(rna_norm)) {
  
  em<- emmreml(rna_norm[, i], within_matrix, re_matZ, rna_kin, varbetahat=T,varuhat=T, PEVuhat=T, test=T)
  
  em_row<- as.data.frame(cbind(em[["betahat"]][1:4,], em[["Xsqtestbeta"]][1:4,],
                           em[["pvalbeta"]][1:4, "none"]))
  
  colnames(em_row)<- c("beta", "chi_square", "pvalue")
  
  em_row$vars<- rownames(em_row)
  
  em_row<- em_row %>%
    pivot_wider(names_from = vars, values_from = c(beta, chi_square, pvalue))
  
  em_row$Vu<- em[["Vu"]]
  
  em_row$Ve<- em[["Ve"]]
  
  em_row$loglik<- em[["loglik"]]
  
  em_row$outcome<- colnames(rna_norm)[i]
  
  em_row<- em_row %>% 
    relocate(outcome, .before=`beta_(Intercept)`)
  
  df<- rbind(df, em_row)
  
  print(paste("Gene", i, "out of", ncol(rna_norm), "done!", sep=" "))
  
}

#Run EMMA for Eq. 3-------------------------------------------------------------
#Generate model matrix
eq3_matrix<- model.matrix(~ trapped_age + mean_age + sex, data = rna_meta)

re_eq = "y ~ trapped_age + mean_age + sex + (1|animal_ID)"

rna_meta$y<- 1

#This generates a random intercept matrix with dims 309 x 114
re_mat <- lFormula(eval(re_eq), rna_meta)
re_matZ <- as.matrix(t(re_mat$reTrms$Zt))

df_eq3<- data.frame(matrix(nrow=0, ncol=16))

for (i in 1:ncol(rna_norm)) {
  
  em<- emmreml(rna_norm[, i], eq3_matrix, re_matZ, rna_kin, varbetahat=T,varuhat=T, PEVuhat=T, test=T)
  
  em_row<- as.data.frame(cbind(em[["betahat"]][1:4,], em[["Xsqtestbeta"]][1:4,],
                               em[["pvalbeta"]][1:4,"none"]))
  
  colnames(em_row)<- c("beta", "chi_square", "pvalue")
  
  em_row$vars<- rownames(em_row)
  
  em_row<- em_row %>%
    pivot_wider(names_from = vars, values_from = c(beta, chi_square, pvalue))
  
  em_row$Vu<- em[["Vu"]]
  
  em_row$Ve<- em[["Ve"]]
  
  em_row$loglik<- em[["loglik"]]
  
  em_row$outcome<- colnames(rna_norm)[i]
  
  em_row<- em_row %>% 
    relocate(outcome, .before=`beta_(Intercept)`)
  
  df_eq3<- rbind(df_eq3, em_row)
  
  print(paste("Gene", i, "out of", ncol(rna_norm), "done!", sep=" "))
  
}

#Run EMMA for chron_age-------------------------------------------------------------
#Generate model matrix
chron_matrix<- model.matrix(~ trapped_age + sex + Seq_batch, data = rna_meta)

re_eq = "y ~ trapped_age + sex + Seq_batch + (1|animal_ID)"

rna_meta$y<- 1

#This generates a random intercept matrix with dims 309 x 114
re_mat <- lFormula(eval(re_eq), rna_meta)
re_matZ <- as.matrix(t(re_mat$reTrms$Zt))

df_chron<- data.frame(matrix(nrow=0, ncol=16))

for (i in 1:ncol(rna_norm)) {
  
  em<- emmreml(rna_norm[, i], chron_matrix, re_matZ, rna_kin, varbetahat=T,varuhat=T, PEVuhat=T, test=T)
  
  em_row<- as.data.frame(cbind(em[["betahat"]][1:4,], em[["Xsqtestbeta"]][1:4,],
                               em[["pvalbeta"]][1:4,"none"]))
  
  colnames(em_row)<- c("beta", "chi_square", "pvalue")
  
  em_row$vars<- rownames(em_row)
  
  em_row<- em_row %>%
    pivot_wider(names_from = vars, values_from = c(beta, chi_square, pvalue))
  
  em_row$Vu<- em[["Vu"]]
  
  em_row$Ve<- em[["Ve"]]
  
  em_row$loglik<- em[["loglik"]]
  
  em_row$outcome<- colnames(rna_norm)[i]
  
  em_row<- em_row %>% 
    relocate(outcome, .before=`beta_(Intercept)`)
  
  df_chron<- rbind(df_chron, em_row)
  
  print(paste("Gene", i, "out of", ncol(rna_norm), "done!", sep=" "))
  
}

df %>%
  select(c(pvalue_within_age, pvalue_mean_age)) %>%
  pivot_longer(cols=c(pvalue_within_age, pvalue_mean_age), 
               values_to = "pval",
               names_to = "term") %>%
  ggplot(aes(pval, fill = term)) +
  geom_histogram(bins=30, colour="black", position = position_dodge()) +
  theme_classic(base_size=12)

df %>%
  ggplot(aes(pvalue_sexM)) +
  geom_histogram(bins=30, colour='black') +
  theme_classic(base_size=12)

df_chron %>%
  ggplot(aes(pvalue_trapped_age)) +
  geom_histogram(bins=30, colour='black') +
  theme_classic(base_size=12)

df_eq3 %>%
  select(c(fdr_trapped_age, fdr_mean_age)) %>%
  pivot_longer(cols=c(fdr_trapped_age, fdr_mean_age), 
               values_to = "pval",
               names_to = "term") %>%
  ggplot(aes(pval, fill = term)) +
  geom_histogram(bins=30, colour="black", position = position_dodge()) +
  theme_classic(base_size=12)

save.image("/home/ckelsey4/Cayo_meth/rna_compare.RData")


em<- emmreml(rna_counts[, 1], predictor_matrix, re_matZ, rna_kin, varbetahat=T,varuhat=T, PEVuhat=T, test=T)

# find number of available cores
ncores <- as.numeric(system("nproc", intern = TRUE))
# create cluster
clus <- makeCluster(ncores)
registerDoParallel(cores = ncores)
# initiate local cluster
clusterExport(clus,
              varlist = c("rna_norm", "rna_meta", "re_matZ", "rna_kin"))


emma_ge_nested_add <- data.frame(parApply(clus, rna_norm, 1, function(y){
  .libPaths(c("~/age_and_aging/",.libPaths()))
  emma <- EMMREML::emmreml(y = y, 
                           X = model.matrix(~ rna_meta$within_age + rna_meta$sex + rna_meta$mean_age + rna_meta$Seq_batch),
                           Z = re_matZ, 
                           K = rna_kin, 
                           varbetahat = T, 
                           varuhat = T, 
                           PEVuhat = T, 
                           test = T)
  return(c(emma$betahat, emma$pvalbeta[, "none"], emma$varbetahat))  
}))



