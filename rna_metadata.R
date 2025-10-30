library(tidyverse)
library(EMMREML)
library(lme4)

rna_meta<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_metadata_wELA_18Jul25.rds")
rna_counts<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_counts_18Jul25.rds")
rna_kin<- readRDS("/scratch/ckelsey4/Cayo_meth/rna_kin_matrix.rds")

#Subset to control and truncate df to relevant variables
rna_meta<- rna_meta %>%
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
rna_counts<- t(rna_counts[, colnames(rna_counts) %in% rna_meta$Sample_ID])

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

#Run PQLseq for within_age-------------------------------------------------------------
#Generate model matrix
predictor_matrix<- model.matrix(~ within_age + mean_age + sex + Seq_batch, data = rna_meta)

re_eq = "y ~ within_age + mean_age + sex + Seq_batch + (1|animal_ID)"

rna_meta$y<- 1

#This generates a random intercept matrix with dims 309 x 114
re_mat <- lFormula(eval(re_eq), rna_meta)
re_matZ <- as.matrix(t(re_mat$reTrms$Zt))

df<- data.frame(matrix(nrow=0, ncol=20))

for (i in 1:ncol(rna_counts)) {
  
  em<- emmreml(rna_counts[, i], predictor_matrix, re_matZ, rna_kin, varbetahat=T,varuhat=T, PEVuhat=T, test=T)
  
  em_row<- as.data.frame(cbind(em[["betahat"]][1:4,], em[["Xsqtestbeta"]][1:4,],
                           em[["pvalbeta"]][1:4,8], em[["pvalbeta"]][1:4,7]))
  
  colnames(em_row)<- c("beta", "chi_square", "pvalue", "fdr")
  
  em_row$vars<- rownames(em_row)
  
  em_row<- em_row %>%
    pivot_wider(names_from = vars, values_from = c(beta, chi_square, pvalue, fdr))
  
  em_row$Vu<- em[["Vu"]]
  
  em_row$Ve<- em[["Ve"]]
  
  em_row$loglik<- em[["loglik"]]
  
  em_row$outcome<- colnames(rna_counts)[i]
  
  em_row<- em_row %>% 
    relocate(outcome, .before=`beta_(Intercept)`)
  
  df<- rbind(df, em_row)
  
}



