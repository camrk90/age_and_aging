library(tidyverse)
library(EMMREML)
library(lme4)

rna_meta<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_metadata_wELA_18Jul25.rds")
rna_counts<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_counts_18Jul25.rds")

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
rna_counts<- rna_counts[, colnames(rna_counts) %in% rna_meta$Sample_ID]


#Run PQLseq for within_age-------------------------------------------------------------
#Generate model matrix
predictor_matrix<- model.matrix(~ within_age + mean_age + sex + Seq_batch, data = rna_meta)

re_eq = "y ~ within_age + mean_age + sex + Seq_batch + (1|animal_ID)"

rna_meta$y<- 1

#This generates a random intercept matrix with dims 309 x 114
re_mat <- lFormula(eval(re_eq), rna_meta)
re_matZ <- t(as.matrix(re_mat$reTrms$Zt))

test<- emmreml(rna_counts, predictor_matrix, re_matZ, )











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


