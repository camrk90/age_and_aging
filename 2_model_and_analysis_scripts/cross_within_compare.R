library(tidyverse)
library(ggplot2)
library(ggcorrplot)
library(ggvenn)
library(ggeffects)
#library(variancePartition)
library(lme4)
library(fgsea)
library(UpSetR)
library(biomaRt)
library(GenomicFeatures)
library(GenomicRanges)
library(msigdbr)

load("/scratch/ckelsey4/Cayo_meth/cross_within_compare.RData")

#Define import function
import_pqlseq<- function(x, y){
  
  #Generate list of file names
  file_list<- list.files(pattern = x)
  file_order<- str_split_i(file_list, "_", y)
  
  #Import glm models as list
  model_list<- lapply(file_list, readRDS)
  
  #Rename list elements
  names(model_list)<- file_order
  model_list<- model_list[1:21]
  
  #Bind model list to df and add rownames
  model<- do.call(rbind, model_list)
  model$outcome<- str_split_i(model$outcome, "\\.", 3)
  model$outcome2<- model$outcome
  
  #Separate region coordinates into start and end, delete the chr col, and move region col to front
  model<- model %>% 
    separate_wider_delim(outcome2, names=c("chr", "chromStart", "chromEnd"), delim = "_") %>%
    relocate(c(chr, chromStart, chromEnd), .after = outcome)
  
  #Add length col and filter by length
  model<- model %>%
    mutate(length = 1+(as.numeric(chromEnd) - as.numeric(chromStart))) %>%
    relocate(length, .after=outcome)
}

######################################
###        Import Metadata         ###
######################################
#Metadata
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

######################################
###          Import data           ###
######################################
#Import longitudinal pqlseq files
setwd('/scratch/ckelsey4/Cayo_meth/glmer_model_compare')

#Chronological Age
chron_age_files<- 'wb_pqlseq2_agechron'
chron_age_pqlseq<- import_pqlseq(chron_age_files, y = 4)

#Age Within
eq2_age_w_files<- 'wb_pqlseq2_within_age'
eq2_age_w_pqlseq<- import_pqlseq(eq2_age_w_files, y = 5)
eq2_age_w_pqlseq<- eq2_age_w_pqlseq[,c(1, 7:15)]

#Mean Age
eq2_age_m_files<- 'wb_pqlseq2_mean_age'
eq2_age_m_pqlseq<- import_pqlseq(eq2_age_m_files, y = 5)
eq2_age_m_pqlseq<- eq2_age_m_pqlseq[,c(1, 7:15)]

#Eq. 3
eq3_age_files<- 'wb_pqlseq2_eq3'
eq3_age_pqlseq<- import_pqlseq(eq3_age_files, y = 4)
eq3_age_pqlseq<- eq3_age_pqlseq[,c(1, 7:15)]

#Eq. 3 Mean Age
eq3_age_m_files<- 'wb_pqlseq2_eq3_m'
eq3_age_m_pqlseq<- import_pqlseq(eq3_age_m_files, y = 5)
eq3_age_m_pqlseq<- eq3_age_m_pqlseq[,c(1, 7:15)]

#Join model dataframes
age_full<- inner_join(chron_age_pqlseq, eq2_age_w_pqlseq, by = "outcome")
age_full<- inner_join(age_full, eq2_age_m_pqlseq, by = "outcome")
age_full<- inner_join(age_full, eq3_age_pqlseq, by = "outcome")
age_full<- inner_join(age_full, eq3_age_m_pqlseq, by = "outcome")

#Cross Sectional Models
setwd('/scratch/ckelsey4/Cayo_meth/cross_models')

long_cross_files<- 'cs_pqlseq2_age_'
long_cross_pqlseq<- import_pqlseq(long_cross_files, y = 4)

long_outcome<- long_cross_pqlseq %>%
  dplyr::select(outcome)

long_cross_pqlseq<- long_cross_pqlseq %>%
  dplyr::select(-c(outcome, length, chr, chromStart, chromEnd, n, )) %>%
  mutate_if(is.character, as.numeric)

long_cross_pqlseq<- cbind(long_outcome, long_cross_pqlseq)

rm(long_outcome)

age_full<- inner_join(age_full, long_cross_pqlseq, by = "outcome")

#Sort chromosome factors
sorted_labels<- str_sort(unique(age_full$chr), numeric=T)

age_full<- age_full %>% 
  mutate(chr = factor(chr, levels = sorted_labels)) %>%
  arrange(chr) %>%
  mutate(region_range = paste(chromStart, "-", chromEnd, sep = " ")) %>%
  relocate(region_range, .after = chromEnd)

age_full$within_cross<- "Both Insignificant"
age_full$within_cross[age_full$fdr_eq2_w_age < 0.05 & age_full$fdr_cross < 0.05]<- "Both Significant"
age_full$within_cross[age_full$fdr_eq2_w_age < 0.05 & age_full$fdr_cross > 0.05]<- "Within Age Significant"
age_full$within_cross[age_full$fdr_eq2_w_age > 0.05 & age_full$fdr_cross < 0.05]<- "Cross Age Significant"

age_full$within_chron<- "Both Insignificant"
age_full$within_chron[age_full$fdr_eq2_w_age < 0.05 & age_full$fdr_chron_age < 0.05]<- "Both Significant"
age_full$within_chron[age_full$fdr_eq2_w_age < 0.05 & age_full$fdr_chron_age > 0.05]<- "Within Age Significant"
age_full$within_chron[age_full$fdr_eq2_w_age > 0.05 & age_full$fdr_chron_age < 0.05]<- "Chron Age Significant"

age_full<- age_full %>%
  mutate(within_chron_diff = abs(beta_chron_age) - abs(beta_eq2_w_age),
         eq3_chron_diff = abs(beta_chron_age) - abs(beta_eq3_age),
         within_chron_ratio = abs(beta_eq2_w_age)/abs(beta_chron_age),
         eq3_chron_ratio = abs(beta_eq3_age)/abs(beta_chron_age))

#Generate age df subset
age_trunc<- age_full %>%
  dplyr::select(c(outcome, region_range, chr, 
           beta_cross, fdr_cross, #cross sectional age
           beta_chron_age, fdr_chron_age, #chron_age
           beta_eq2_w_age, fdr_eq2_w_age, #eq2 within age
           beta_eq2_m_age, fdr_eq2_m_age, #eq2 mean age
           beta_eq3_age, fdr_eq3_age, #eq3 age
           beta_eq3_age_m, fdr_eq3_age_m, #eq3 mean age 
           within_chron_diff, eq3_chron_diff, #diffs
           within_chron_ratio, eq3_chron_ratio)) #ratios

rm(eq2_age_w_pqlseq);rm(eq2_age_m_pqlseq);rm(chron_age_pqlseq)
rm(long_cross_pqlseq);rm(eq3_age_pqlseq);rm(eq3_age_m_pqlseq)

######################################
###       Descriptive Stats        ###
######################################
#Cross-age age distribution
blood_metadata %>%
  ggplot(aes(age_at_sampling, fill=individual_sex)) +
  geom_histogram(position = "dodge", colour = "black") +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  scale_y_continuous(breaks = seq(0, 30, 5)) +
  coord_cartesian(ylim = c(0, 30)) +
  scale_fill_manual(values = c("goldenrod1", "goldenrod4"), name = "Sex") +
  theme_classic(base_size=24) +
  theme(legend.key.height= unit(2, 'cm')) +
  theme(panel.background = element_rect(colour = "black", linewidth=2)) +
  ylab("Count") +
  xlab("Age")

blood_metadata %>%
  ggplot(aes(individual_sex, fill=individual_sex)) +
  geom_bar(position = "dodge", colour = "black") +
  scale_fill_manual(values = c("goldenrod1", "goldenrod4"), name = "Sex") +
  theme_classic(base_size=24) +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(colour = "black", linewidth=2)) +
  ylab("Count") +
  xlab("Sex")

ggsave("/home/ckelsey4/Cayo_meth/age_aging/cross_age.svg", plot = last_plot())

#Longitudinal data distribution
long_data<- long_data %>%
  group_by(monkey_id) %>%
  mutate(min_age = min(age_at_sampling))
long_data$age_at_sampling<- round(long_data$age_at_sampling, 0)

samples_dist<- long_data %>%
  ggplot(aes(x=age_at_sampling, y=reorder(monkey_id, min_age), colour=individual_sex)) +
  geom_path(linewidth = 0.5) +
  geom_point(colour="black", size = 0.25) +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  coord_cartesian(xlim = c(0, 30)) +
  scale_colour_manual(values = c("red3", "pink2"), name = "Sex") +
  ylab("Individual") +
  xlab("Age") +
  theme_classic() +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        legend.position = "none")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/samples_dist.svg", 
       samples_dist, 
       height = 80, width = 60, units = "mm")

#Age distribution
long_data %>%
  ggplot(aes(x=round(age_at_sampling), fill=as.factor(individual_sex))) +
  geom_bar(colour='black', position = 'dodge') +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  coord_cartesian(xlim = c(0, 30)) +
  scale_fill_manual(values = c("red3", "pink2"), name = "Sex") +
  ylab("Count") +
  xlab("Age") +
  theme_classic(base_size = 24) +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(colour = "black", linewidth=3))

#N Samples
n_samples<- long_data %>%
  ggplot(aes(x=n, fill=as.factor(individual_sex))) +
  geom_bar(position = 'dodge') +
  scale_fill_manual(values = c("red3", "pink2"), name = "Sex") +
  #ylab("Count") +
  #xlab("N Samples") +
  theme_classic() +
  theme(panel.background = element_rect(colour = "black", linewidth=0.5),
        axis.line = element_line(colour = "black", linewidth = 0.25),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        legend.position = "none")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/n_samples.svg", 
       n_samples, 
       height = 16, width = 24, units = "mm")


age_full %>%
  select(c(beta_within_age, beta_chron_age, beta_long_cross, beta_mean_age)) %>%
  cor(use="pairwise.complete.obs") %>%
  ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=5, sig.level = 0.05, insig = "blank")

#Variance Partition-------------------------------------------------------------
#Import m/cov rds
regions_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_filtered.rds")
regions_m<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_m_filtered.rds")

regions_m<- do.call(rbind, regions_m)
regions_cov<- do.call(rbind, regions_cov)

p_meth<- regions_m/regions_cov

ratio_matrix<- as.matrix(p_meth)
ratio_matrix[!is.finite(ratio_matrix)]<- 0
ratio_matrix[is.na(ratio_matrix)]<- 0
ratio_matrix[is.nan(ratio_matrix)]<- 0

pca_blood<- prcomp(cor(p_meth, use="pairwise.complete.obs"))
pcs<- as.data.frame(pca_blood$x) 
summary(pca_blood)$importance[2, ]

pcs<- pcs[rownames(pcs) %in% blood_metadata$lid_pid,]
df<- blood_metadata[blood_metadata$lid_pid %in% rownames(pcs),]



meta<- long_data[long_data$lid_pid %in% colnames(p_meth),]

ratio_matrix<- ratio_matrix[,meta$lid_pid]

vp_model<- ~ within.age + mean.age + (1|individual_sex) + perc_unique

vp<- fitExtractVarPartModel(ratio_matrix, vp_model, meta)

plotVarPart(vp)


#P-Value Distributions----------------------------------------------------------
age_full %>%
  ggplot(aes(pvalue_cross)) +
  geom_histogram(bins=30, colour='black', fill = 'goldenrod2') +
  theme_classic(base_size=24) +
  theme(legend.position = "none") +
  xlab("P-Value") +
  ylab("Count")

age_full %>%
  ggplot(aes(pvalue_chron_age)) +
  geom_histogram(bins=30, colour='black', fill = 'steelblue2') +
  theme_classic(base_size=24) +
  theme(legend.position = "none") +
  xlab("P-Value") +
  ylab("Count")

age_full %>%
  ggplot(aes(pvalue_eq2_w_age)) +
  geom_histogram(bins=30, colour='black', fill = 'green4') +
  theme_classic(base_size=24) +
  theme(legend.position = "none") +
  xlab("P-Value") +
  ylab("Count")

age_full %>%
  ggplot(aes(pvalue_eq3_age)) +
  geom_histogram(bins=30, colour='black', fill = 'purple') +
  theme_classic(base_size=24) +
  theme(legend.position = "none") +
  xlab("P-Value") +
  ylab("Count")

age_full %>%
  select(c(pvalue_eq2_w_age, pvalue_eq2_mean_age)) %>%
  pivot_longer(cols=c(pvalue_eq2_w_age, pvalue_eq2_mean_age), 
               values_to = "pval",
               names_to = "term") %>%
  ggplot(aes(pval, fill = term)) +
  geom_histogram(bins=30, colour="black", position = position_dodge()) +
  scale_fill_manual(values = c("black", "green4")) +
  theme_classic(base_size=24) +
  theme(legend.position = "none") +
  xlab("P-Value") +
  ylab("Count")

age_full %>%
  select(c(pvalue_eq3_age, pvalue_eq3_age_m)) %>%
  pivot_longer(cols=c(pvalue_eq3_age, pvalue_eq3_age_m), 
               values_to = "pval",
               names_to = "term") %>%
  ggplot(aes(pval, fill = term)) +
  geom_histogram(bins=30, colour="black", position = position_dodge()) +
  scale_fill_manual(values = c("purple", "black")) +
  theme_classic(base_size=24) 
  #theme(legend.position = "none")


#FDR Thresholds
thresholds <- c(0.01, 0.05, 0.10, 0.15, 0.20)

results <- sapply(age_trunc[, c("fdr_cross", "fdr_chron_age", "fdr_eq2_w_age", "fdr_eq3_age")], function(pvals) {
  sapply(thresholds, function(t) sum(pvals < t, na.rm = TRUE))
})

results_df <- as.data.frame(results)
rownames(results_df) <- paste0("FDR_", thresholds)
results_df$fdr_threshold<- rownames(results_df)

results_df<- results_df %>%
  pivot_longer(cols = c(fdr_cross, fdr_chron_age, fdr_eq2_w_age, fdr_eq3_age),
               names_to = "model",
               values_to = "count")
results_df$fdr_threshold<- as.numeric(str_split_i(results_df$fdr_threshold, "_", 2))

results_df %>%
  ggplot(aes(fdr_threshold, count, colour = model)) +
  geom_point() +
  geom_path() +
  scale_colour_manual(values = c("darkgoldenrod2","steelblue1", "green4", "purple"), name = "") +
  theme_classic() +
  theme(panel.background = element_rect(colour = "black", linewidth=0.5),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  theme(panel.grid.major = element_line(color = "grey90", size = 0.5),
        panel.grid.minor = element_line(color = "grey98", size = 0.5)) 

#### Coefficient comparisons ---------------------------------------------------
compare_plot<- function(df, fdr1, fdr2, var1, var2, c1, c2, plot_type) {
  
  df<- df %>%
    filter({{fdr1}} < .05 & {{fdr2}} < .05) %>%
    mutate(diff = abs({{var2}}) - abs({{var1}}),
           ratio = abs({{var2}})/abs({{var1}}))
  
  print(paste("The median effect size difference =", median(df$ratio), sep = " "))
  correlation<- cor(df %>% pull({{var1}}), df %>% pull({{var2}}))
  
  print(paste("The correlation between", 
              deparse(substitute(var1)), "and", 
              deparse(substitute(var2)), "=", 
              correlation,
              sep = " "))
  print(nrow(df))
  print(nrow(df[df$ratio < 1,]))
  
  if (plot_type == "scatter"){
    
    df %>%
      ggplot(aes({{var1}}, {{var2}}, colour = diff)) +
      geom_point(size = 0.25, alpha = 0.8) +
      geom_abline() +
      geom_smooth(method = "lm") +
      geom_vline(xintercept=0, linetype="dashed") +
      geom_hline(yintercept=0, linetype="dashed") +
      scale_color_gradient2(low = c2, mid = "grey70", high = c1, midpoint = 0, name = "") +
      theme_classic() +
      theme(legend.key.width = unit(1, 'cm'), 
            legend.position = "none") +
      theme(panel.background = element_rect(colour = "black", linewidth=1),
            axis.line = element_line(colour = "black", linewidth = 0.5),
            axis.title = element_blank(),
            axis.text = element_blank(),
            plot.margin = margin(0, 0, 0, 0, "pt"),
            aspect.ratio = 1) +
      xlim(-0.20, 0.20) +
      ylim(-0.20, 0.20) 
    
  } else if (plot_type == "hist") {
    
    df %>%
      ggplot(aes(ratio, fill = after_stat(x))) +
      geom_histogram(bins = 50) +
      geom_vline(xintercept=1, linetype="dashed") +
      geom_vline(xintercept=median(df$ratio), linetype="dashed", colour = 'red') +
      scale_fill_gradient2(low = c2, mid = "grey70", high = c1, midpoint = 1, name = "") +
      theme_classic() +
      theme(legend.position = "none",
            panel.background = element_rect(colour = "black", linewidth=1),
            axis.line = element_line(colour = "black", linewidth = 0.5),
            axis.title = element_blank(),
            axis.text = element_blank(),
            plot.margin = margin(0, 0, 0, 0, "pt"),
            aspect.ratio = 1) +
      scale_x_continuous(breaks = seq(0, 7, by=1)) +
      coord_cartesian(xlim = c(0,7))
    
  }
}

compare_plot(age_trunc, fdr_chron_age, fdr_cross, 
             beta_chron_age, beta_cross, 
             "darkgoldenrod2", "steelblue2", "scatter")

compare_plot(age_full, fdr_chron_age, fdr_cross,
             beta_chron_age, beta_cross,
             "darkgoldenrod2", "steelblue2", "hist")

within_chron_plot<- compare_plot(age_trunc, fdr_chron_age, fdr_eq2_w_age,
                                 beta_chron_age, beta_eq2_w_age,
                                 "green4", "steelblue2", "scatter")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/within_chron_scatter.svg", 
       within_chron_plot, 
       height = 45, width = 45, units = "mm")

within_chron_hist<- compare_plot(age_trunc, fdr_chron_age, fdr_eq2_w_age,
                                beta_chron_age, beta_eq2_w_age,
                                "green4", "steelblue2", plot_type = "hist")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/within_chron_hist2.svg", 
       within_chron_hist, 
       height = 45, width = 45, units = "mm")


eq3_chron_plot<- compare_plot(age_trunc, fdr_chron_age, fdr_eq3_age,
                              beta_chron_age, beta_eq3_age,
                              "purple", "steelblue2", "scatter")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq3_chron_scatter.svg", 
       eq3_chron_plot, 
       height = 45, width = 45, units = "mm")

eq3_chron_hist<- compare_plot(age_trunc, fdr_chron_age, fdr_eq3_age,
                              beta_chron_age, beta_eq3_age,
                              "purple", "steelblue2", "hist")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq3_chron_hist.svg", 
       eq3_chron_hist, 
       height = 45, width = 45, units = "mm")

compare_plot(age_trunc, fdr_eq2_w_age, fdr_eq3_age,
             beta_eq2_w_age, beta_eq3_age,
             "purple", "green4", "scatter")

compare_plot(age_full, fdr_eq2_w_age, fdr_eq3_age,
             beta_eq2_w_age, beta_eq3_age,
             "purple", "green4", "hist")

age.cross.count<- nrow(age_full[age_full$fdr_cross < 0.05,])
age.chron.count<- nrow(age_full[age_full$fdr_chron_age < 0.05,])
age.w.count<- nrow(age_full[age_full$fdr_eq2_w_age < 0.05,])
age.eq3.count<- nrow(age_full[age_full$fdr_eq3_age < 0.05,])
counts<- data.frame(predictor = as.factor(c('Eq2. Within Age', 'Eq3 Age', 'Chron Age', 'Cross Age')),
                    count = c(age.w.count, age.eq3.count, age.chron.count, age.cross.count))

counts<- counts %>%
  mutate(predictor = as.factor(predictor)) %>%
  mutate(predictor=fct_reorder(predictor, count, .desc=T),
         perc_signif = count/nrow(age_full))

df<- age_trunc %>%
  pivot_longer(cols = c(beta_eq2_w_age, beta_eq3_age, beta_chron_age, beta_cross),
               values_to = 'beta',
               names_to = 'var') %>%
  group_by(var) %>%
  summarize(mean = mean(beta),
            sd = sd(beta),
            var = var(beta))

rm(age.w.count);rm(age.cross.count);rm(age.chron.count);rm(age.eq3.count)

signif_counts<- counts %>%
  filter(!predictor == "Cross Age") %>%
  ggplot(aes(predictor, count, fill = predictor)) +
  geom_bar(stat = 'identity') +
  #geom_text(label=counts$count, vjust=-1, size=5) +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  xlab("Predictor") +
  ylab("Count") +
  scale_fill_manual(values = c("purple", 'green4', 'steelblue2'))
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/signif_counts.svg", 
       signif_counts, 
       height = 20, width = 30, units = "mm")

#Bar plot for all 4 models including cross-sectional
counts %>%
  ggplot(aes(predictor, count, fill = predictor)) +
  geom_bar(stat = 'identity', colour="black") +
  geom_text(label=counts$count, vjust=-1, size=5) +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        aspect.ratio = 1) +
  xlab("Predictor") +
  ylab("Count") +
  scale_fill_manual(values = c("purple", 'green4', 'steelblue2', 'darkgoldenrod2'))

#Distribution of effect sizes for each variable
beta_dist<- age_full %>%
  dplyr::select(c(beta_eq2_w_age, beta_eq3_age, beta_chron_age)) %>%
  pivot_longer(cols = c(beta_eq2_w_age, beta_eq3_age, beta_chron_age),
               values_to = 'beta',
               names_to = 'var') %>%
  ggplot(aes(beta, fill=var)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_fill_manual(values = c('steelblue2', 'green4', "purple"),
                    labels = c("Chron Age", "Eq2. Within Age", "Eq3. Age")) +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  xlim(-0.4, 0.2) +
  ylab("Density") +
  xlab("Beta")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/dnam_beta_dist.svg", 
       beta_dist, 
       height = 45, width = 60, units = "mm")

t.test(age_trunc$beta_eq3_age, age_trunc$beta_chron_age)
t.test(age_trunc$beta_eq2_w_age, age_trunc$beta_chron_age)

## Significant regions Venn diagram
cross.age<- age_trunc$outcome[age_trunc$fdr_cross< 0.05]
chron.age<- age_trunc$outcome[age_trunc$fdr_chron_age< 0.05]
age.w<- age_trunc$outcome[age_trunc$fdr_eq2_w_age < 0.05]
eq2.btwn<- age_trunc$outcome[age_trunc$fdr_eq2_mean_age < 0.05]
eq3<- age_trunc$outcome[age_trunc$fdr_eq3_age < 0.05]

venn_list2<- list(cross.age, chron.age, age.w)
names(venn_list2)<- c("cross.age", "chron.age", "age.w")

ggvenn(venn_list2,
       text_size = 8,
       show_percentage = F)

venn_all<- list(cross.age, chron.age, age.w, eq3)
names(venn_all)<- c("cross.age", "chron.age", "age.w", "eq3")

upset(fromList(venn_all), order.by = "freq", 
      text.scale = c(2, 2, 2, 1, 2, 1.5), 
      line.size = 1, point.size = 2)
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/n_samples.svg", 
       n_samples, 
       height = 16, width = 24, units = "mm")

######################################
###      JOIN INTERSECT FILES      ###   
######################################
##Import annotation files-------------------------------------------------------
re_anno<- read_csv("/scratch/ckelsey4/Cayo_meth/re_annotations.csv")
re_anno<- re_anno %>%
  filter(chr != "Y")
chmm_intersect<- read_csv("/scratch/ckelsey4/Cayo_meth/chmm_annotations.csv")
chmm_intersect<- chmm_intersect %>%
  filter(chr != "Y")
promoters<- read_csv("/scratch/ckelsey4/Cayo_meth/promoters.csv")
promoters<- promoters %>%
  filter(chr != "Y")

#Promoters----------------------------------------------------------------------
pqlseq_prom<- left_join(promoters, age_trunc, by = c("region_range", "chr"))
pqlseq_prom<- pqlseq_prom %>%
  drop_na() %>%
  distinct(anno, .keep_all = T)
pqlseq_prom$anno_class<- "Promoter"

#Create new promoter df where 'anno' = "Promoter" not the gene name for factor issues later
pqlseq_prom2<- pqlseq_prom %>%
  mutate(anno = "Promoter")

#CHMM---------------------------------------------------------------------------
#Join pqlseq model and chmm 
pqlseq_chmm<- left_join(chmm_intersect, age_trunc, by = c("region_range", "chr"))

#Filter out regions with models that didn't converge resulting in NAs in the annotation join
pqlseq_chmm<- pqlseq_chmm %>%
  drop_na()

#Set annotations as factor and reorder
chmm_ordered<- as.factor(str_sort(unique(pqlseq_chmm$anno), numeric = TRUE))
pqlseq_chmm$anno<- factor(pqlseq_chmm$anno, levels = rev(chmm_ordered))

#Create column of broad categories for annotations
pqlseq_chmm$anno_class<- "TSSs"
pqlseq_chmm$anno_class[pqlseq_chmm$anno %in% chmm_ordered[3:5]]<- "Active Tr."
pqlseq_chmm$anno_class[pqlseq_chmm$anno %in% chmm_ordered[6:8]]<- "Enhancers"
pqlseq_chmm$anno_class[pqlseq_chmm$anno %in% chmm_ordered[9:15]]<- "Quiescent"

#Set classes as factors
class_factors<- c("TSSs", "Active Tr.", "Enhancers", "Quiescent")
pqlseq_chmm$anno_class<- factor(pqlseq_chmm$anno_class, levels = rev(class_factors))

#REPEAT ELEMENTS----------------------------------------------------------------
#Join repeats annotations and glm_models df
pqlseq_re<- left_join(re_anno, age_trunc, by = c("region_range", "chr"))

pqlseq_re<- pqlseq_re %>%
  drop_na()

#Remove non-sensical annotations (NA, Unknown etc)
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "Unknown",]
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "DNA?",]
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "LTR?",]
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "RC?",]
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass == "Unspecified",]

pqlseq_re$anno_class<- "Simple Repeats"
pqlseq_re$anno_class[pqlseq_re$repClass %in% c("SINE", "LINE", "LTR", "Retroposon")]<- "TE Class I"
pqlseq_re$anno_class[pqlseq_re$repClass %in% "DNA"]<- "TE Class II"
pqlseq_re<- pqlseq_re[!pqlseq_re$repClass %in% c("RC", "rRNA", "snRNA", "tRNA", "srpRNA", "scRNA", "Low_complexity"),]

re_ordered<- c("Simple_repeat", "Satellite", "SINE", "LINE", "LTR", "Retroposon", "DNA")

pqlseq_re<- pqlseq_re %>%
  dplyr::rename(anno = repClass) %>%
  dplyr::select(-repName) %>%
  dplyr::relocate(anno, .after = anno_end) %>%
  dplyr::select(-c(range))

pqlseq_re$anno<- factor(pqlseq_re$anno, levels = rev(re_ordered))

#Bind annotation dfs together---------------------------------------------------
pqlseq_anno<- rbind(pqlseq_chmm, pqlseq_re, pqlseq_prom2)

annos<- c("Promoter", "TSSs", "Active Tr.", "Enhancers", "Quiescent", "Simple Repeats",
          "TE Class I", "TE Class II")

pqlseq_anno$anno_source<- "Repeat Elements"
pqlseq_anno$anno_source[pqlseq_anno$anno_class == "TSSs" | pqlseq_anno$anno_class == "Active Tr." |
                          pqlseq_anno$anno_class == "Enhancers" | pqlseq_anno$anno_class == "Quiescent" | 
                          pqlseq_anno$anno_class == "Promoter"]<- "Transcription"

pqlseq_anno<- pqlseq_anno %>%
  arrange(anno_source, anno_class)

#Rearrange factors to sort by type then log_or
pqlseq_anno$anno_class<- factor(pqlseq_anno$anno_class, levels = rev(annos))

pqlseq_anno$unique_cpg<- paste(pqlseq_anno$chr, pqlseq_anno$cpg_loc, sep="_")

pqlseq_anno<- pqlseq_anno %>%
  dplyr::relocate(c(anno_class, anno_source, unique_cpg), .after=anno)

rm(chmm_intersect);rm(re_anno);rm(promoters)

#Annotation proportions
pqlseq_anno$cross_signif<- "Non-Significant"
pqlseq_anno$cross_signif[pqlseq_anno$fdr_cross < 0.05 & pqlseq_anno$beta_cross < 0]<- "Age-Hypomethylated"
pqlseq_anno$cross_signif[pqlseq_anno$fdr_cross < 0.05 & pqlseq_anno$beta_cross > 0]<- "Age-Hypermethylated"

pqlseq_anno$cross_signif<- factor(pqlseq_anno$cross_signif, 
                        levels = c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))

pqlseq_anno$chron_signif<- "Non-Significant"
pqlseq_anno$chron_signif[pqlseq_anno$fdr_chron_age < 0.05 & pqlseq_anno$beta_chron_age < 0]<- "Age-Hypomethylated"
pqlseq_anno$chron_signif[pqlseq_anno$fdr_chron_age < 0.05 & pqlseq_anno$beta_chron_age > 0]<- "Age-Hypermethylated"

pqlseq_anno$chron_signif<- factor(pqlseq_anno$chron_signif, 
                                  levels = c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))

pqlseq_anno$eq2_w_signif<- "Non-Significant"
pqlseq_anno$eq2_w_signif[pqlseq_anno$fdr_eq2_w_age < 0.05 & pqlseq_anno$beta_eq2_w_age < 0]<- "Age-Hypomethylated"
pqlseq_anno$eq2_w_signif[pqlseq_anno$fdr_eq2_w_age < 0.05 & pqlseq_anno$beta_eq2_w_age > 0]<- "Age-Hypermethylated"

pqlseq_anno$eq2_w_signif<- factor(pqlseq_anno$eq2_w_signif, 
                                  levels = c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))

pqlseq_anno$eq2_m_signif<- "Non-Significant"
pqlseq_anno$eq2_m_signif[pqlseq_anno$fdr_eq2_mean_age < 0.05 & pqlseq_anno$beta_eq2_mean_age < 0]<- "Age-Hypomethylated"
pqlseq_anno$eq2_m_signif[pqlseq_anno$fdr_eq2_mean_age < 0.05 & pqlseq_anno$beta_eq2_mean_age > 0]<- "Age-Hypermethylated"

pqlseq_anno$eq2_m_signif<- factor(pqlseq_anno$eq2_m_signif, 
                                    levels = c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))

pqlseq_anno$eq3_age_signif<- "Non-Significant"
pqlseq_anno$eq3_age_signif[pqlseq_anno$fdr_eq3_age < 0.05 & pqlseq_anno$beta_eq3_age < 0]<- "Age-Hypomethylated"
pqlseq_anno$eq3_age_signif[pqlseq_anno$fdr_eq3_age < 0.05 & pqlseq_anno$beta_eq3_age > 0]<- "Age-Hypermethylated"

pqlseq_anno$eq3_age_signif<- factor(pqlseq_anno$eq3_age_signif, 
                                   levels = c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))

#save pqlseq_anno df to .rds for gsea script
saveRDS(pqlseq_anno, "/scratch/ckelsey4/Cayo_meth/pqlseq_anno.rds")

#Plot annotation proportions----------------------------------------------------
generate_proportion<- function(df, x, c1, c2, c3){
  
  d1<- df %>% 
    distinct(unique_cpg, .keep_all = T) %>%
    group_by(anno, {{x}}) %>% 
    summarise(count = n()) %>% 
    mutate(perc = count/sum(count))
  
  print(d1)
  
  d2<- df %>%
    distinct(unique_cpg, .keep_all = T) %>%
    group_by({{x}}) %>%
    summarise(count = n()) %>%
    mutate(perc = count/sum(count))
  
  d2$anno<- "All"
  
  d3<- rbind(d1, d2)
  annos2<- unique(d3$anno)
  d3$anno<- factor(d3$anno, levels = annos2)
  
  d3 %>%
    filter(is.na(anno)) %>%
    dplyr::select(count)
  
  col1 <- eval(substitute(x), d2)
  
  d3$percent<- d3$perc*100
  
  d3 %>%
    arrange(anno) %>%
    ggplot(aes(x = percent, y=anno, fill = factor({{x}}))) +
    geom_bar(stat="identity", width = 0.7, colour="black") +
    theme_classic(base_size=24) +
    geom_vline(xintercept = (1-d2$perc[col1 == "Age-Hypermethylated"])*100, linetype = 'dashed') +
    geom_vline(xintercept = d2$perc[col1 == "Age-Hypomethylated"]*100, linetype = 'dashed') +
    theme(legend.position = "top") +
    scale_fill_manual(values = c(c1, c2, c3), name = "") +
    ylab("Annotation") +
    xlab("Percentage")
  
}

generate_proportion(pqlseq_anno, cross_signif, "darkgoldenrod1", "gray90", "darkgoldenrod4")

generate_proportion(pqlseq_anno, chron_signif, "steelblue1", "gray90", "steelblue4")

generate_proportion(pqlseq_anno, eq2_w_signif, "green1", "gray90", "green4")

generate_proportion(pqlseq_anno, eq3_age_signif, "purple1", "gray90", "purple4")


#Promoter Gene Enrichment Analyses----------------------------------------------
#Generate hallmark gene set
hallmark.msigdb = msigdbr(species = "Macaca mulatta", category = "H")
hallmark_list = split(x = hallmark.msigdb$ensembl_gene, f = hallmark.msigdb$gs_name)

proms_chron<- pqlseq_prom %>%
  dplyr::select(anno, beta_chron_age) %>%
  arrange(desc(beta_chron_age))

proms_chron2<- proms_chron$beta_chron_age
names(proms_chron2) = proms_chron$anno

proms_chron_gsea<- fgsea(pathways = hallmark_list, 
                   stats = proms_chron2)

proms_eq3<- pqlseq_prom %>%
  dplyr::select(anno, beta_eq3_age) %>%
  arrange(desc(beta_eq3_age))

proms_eq3_2<- proms_eq3$beta_eq3_age
names(proms_eq3_2) = proms_eq3$anno

proms_chron_gsea<- fgsea(pathways = hallmark_list, 
                         stats = proms_eq3_2,
                         minSize = 1,
                         maxSize = 1000,
                         eps = 0.0)

#### Rank Enrichment ####
gsea<- readRDS('/scratch/ckelsey4/Cayo_meth/within_gsea.rds')

gsea$type<- "Repeat Elements"
gsea$type[gsea$pathway %in% chmm_ordered]<- "Chromatin States"
gsea$type[gsea$pathway == "Promoter"]<- "Chromatin States"

gsea_chmm<- gsea %>%
  filter(type == "Chromatin States") %>%
  mutate(pathway = factor(pathway, levels = pathway))

gsea_chmm %>%
  ggplot(aes(x=pathway, y=NES, fill = NES < 0)) +
  #geom_point(aes(alpha=padj<0.05)) +
  #geom_line(aes(group = type)) +
  geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  scale_fill_manual(values = c("green4", 'steelblue2')) +
  theme_classic(base_size =20) +
  #theme(legend.position = "none") +
  #ylim(c(-1, 7)) +
  ylab("NES") +
  xlab("Annotation") +
  coord_flip()

#### Fishers Exact ####
enrichment<- function(model_df, model_type1, model_type2, var_type1, var_type2){
  
  df_list<- list()
  
  for(i in unique(model_df$anno)) {
    
    df2<- model_df[model_df$anno == i,]
    df2<- df2 %>% distinct(unique_cpg, .keep_all = T)
    
    a<- nrow(df2[df2[[model_type1]] == var_type1 & df2[[model_type2]] == var_type2,])
    c<- nrow(df2[!c(df2[[model_type1]] == var_type1 & df2[[model_type2]] == var_type2),])
    
    df3<- model_df[!model_df$unique_cpg %in% df2$unique_cpg,]
    df3<- df3 %>% distinct(unique_cpg, .keep_all = T)
    
    b<- nrow(df3[df3[[model_type1]] == var_type1 & df3[[model_type2]] == var_type2,])
    d<- nrow(df3[!c(df3[[model_type1]] == var_type1 & df3[[model_type2]] == var_type2),])
    
    #Generate contingency table
    c_table<- data.frame("x" = c(a, b),
                         "y" = c(c, d),
                         row.names = c(paste(i, "Y", sep=""), paste(i, "N", sep="")))
    
    colnames(c_table) = c("Is EQ2/Hypo", "Is NOT EQ2/hypo")
    
    if (all.equal(sum(c_table), length(unique(model_df$unique_cpg))) == T){
      print(paste("Contingency table sum for", i, "matches unique cpg_loc length"))
      
      df_list[[length(df_list)+1]] = c_table
      
      print(c_table)
    } else {
      print(paste("Contingency table sum for", i, "DOES NOT MATCH unique cpg_loc length"))
    }
  }
  #name table list
  names(df_list)<- unique(model_df$anno)
  
  #Fisher test for each table and tidy with broom
  ft<- lapply(df_list, fisher.test)
  ft<- lapply(ft, broom::tidy)
  
  ft<- do.call(rbind, ft)
  ft<- ft %>%
    mutate(annotation = rownames(ft))
  
  #FDR p-val adjustment
  ft<- ft %>%
    mutate(padj = p.adjust(p.value)) %>%
    mutate_at(vars(annotation), as.factor)
  
  #Log estimates and CIs
  ft<- ft %>%
    mutate(log_or = log(estimate),
           log_ci.lo = log(conf.low),
           log_ci.hi = log(conf.high))
  
  annos_order<- str_sort(ft$annotation, numeric = TRUE)
  
  ft$type<- var_type2
  
  ft$source<- "RE's"
  ft$source[ft$annotation %in% chmm_ordered]<- "CHMM"
  
  return(ft)
}

#Generate col for Eq.2 vs Eq.1 signif
pqlseq_anno$within_chron<- "Both Insignificant"
pqlseq_anno$within_chron[pqlseq_anno$fdr_eq2_w_age < 0.05 & pqlseq_anno$fdr_chron_age < 0.05]<- "Both Significant"
pqlseq_anno$within_chron[pqlseq_anno$fdr_eq2_w_age < 0.05 & pqlseq_anno$fdr_chron_age > 0.05]<- "Within Age Significant"
pqlseq_anno$within_chron[pqlseq_anno$fdr_eq2_w_age > 0.05 & pqlseq_anno$fdr_chron_age < 0.05]<- "Chron Age Significant"

#Generate col for Eq.3 vs Eq.1 signif
pqlseq_anno$eq3_chron<- "Both Insignificant"
pqlseq_anno$eq3_chron[pqlseq_anno$fdr_eq3_age < 0.05 & pqlseq_anno$fdr_chron_age < 0.05]<- "Both Significant"
pqlseq_anno$eq3_chron[pqlseq_anno$fdr_eq3_age < 0.05 & pqlseq_anno$fdr_chron_age > 0.05]<- "Eq3 Age Significant"
pqlseq_anno$eq3_chron[pqlseq_anno$fdr_eq3_age > 0.05 & pqlseq_anno$fdr_chron_age < 0.05]<- "Chron Age Significant"

#Hypo Enrichment
within_chron_hypo<- enrichment(pqlseq_anno, 
                                "within_chron", "eq2_w_signif",
                                "Within Age Significant", "Age-Hypomethylated")
 
within_chron_hypo_chmm<- within_chron_hypo %>% filter(source == "CHMM")  %>% mutate(model = "EQ2")
within_chron_hypo_chmm$annotation<- factor(within_chron_hypo_chmm$annotation, 
                                            levels = rev(chmm_ordered)) 

eq3_chron_hypo<- enrichment(pqlseq_anno, 
                               "eq3_chron", "eq3_age_signif",
                               "Eq3 Age Significant", "Age-Hypomethylated")

eq3_chron_hypo_chmm<- eq3_chron_hypo %>% filter(source == "CHMM") %>% mutate(model = "EQ3")
eq3_chron_hypo_chmm$annotation<- factor(eq3_chron_hypo_chmm$annotation, 
                                          levels = rev(chmm_ordered))

chron_hypo<- enrichment(pqlseq_anno, "chron_signif", "Age-Hypomethylated")
chron_hypo$anno_source<- gsub("Chromatin States", "CHMM", chron_hypo$anno_source)
chron_hypo_chmm<- chron_hypo %>% filter(anno_source == "CHMM")
chron_hypo_chmm<- chron_hypo_chmm %>% 
  dplyr::rename(source = anno_source) %>%
  mutate(model = "CHRON")
chron_hypo_chmm$annotation<- factor(chron_hypo_chmm$annotation, 
                                        levels = rev(chmm_ordered))

hypo_enrich<- rbind(within_chron_hypo_chmm, eq3_chron_hypo_chmm, chron_hypo_chmm)

hypo_enrich_plot<- hypo_enrich %>%
  ggplot(aes(x=annotation, y=estimate, colour = model, alpha=padj<0.05)) +
  geom_point(aes(shape = model),size = 1) +
  geom_line(aes(group = model)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  #scale_colour_manual(values = c("purple1", "purple4"), name = "") +
  geom_errorbar(ymin = hypo_enrich$conf.low, ymax = hypo_enrich$conf.high, width = 0.3) +
  scale_colour_manual(values = c("steelblue2", "green4", "purple"), name = "") +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=0.5),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  ylab("Odds Ratio") +
  xlab("Annotation") +
  scale_y_continuous(breaks = seq(0,10,2), limits = c(0, 10))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/hypo_enrich.svg", 
       hypo_enrich_plot, 
       height = 40, width = 105, units = "mm")

#Hyper Enrichment
within_chron_hyper<- enrichment(pqlseq_anno, 
                               "within_chron", "eq2_w_signif",
                               "Within Age Significant", "Age-Hypermethylated")

within_chron_hyper_chmm<- within_chron_hyper %>% filter(source == "CHMM")  %>% mutate(model = "EQ2")
within_chron_hyper_chmm$annotation<- factor(within_chron_hyper_chmm$annotation, 
                                           levels = rev(chmm_ordered)) 

eq3_chron_hyper<- enrichment(pqlseq_anno, 
                            "eq3_chron", "eq3_age_signif",
                            "Eq3 Age Significant", "Age-Hypermethylated")

eq3_chron_hyper_chmm<- eq3_chron_hyper %>% filter(source == "CHMM") %>% mutate(model = "EQ3")
eq3_chron_hyper_chmm$annotation<- factor(eq3_chron_hyper_chmm$annotation, 
                                        levels = rev(chmm_ordered))

chron_hyper<- enrichment(pqlseq_anno, "chron_signif", "Age-hypermethylated")
chron_hyper$anno_source<- gsub("Chromatin States", "CHMM", chron_hyper$anno_source)
chron_hyper_chmm<- chron_hyper %>% filter(anno_source == "CHMM")
chron_hyper_chmm<- chron_hyper_chmm %>% 
  dplyr::rename(source = anno_source) %>%
  mutate(model = "CHRON")
chron_hyper_chmm$annotation<- factor(chron_hyper_chmm$annotation, 
                                    levels = rev(chmm_ordered))

hyper_enrich<- rbind(within_chron_hyper_chmm, eq3_chron_hyper_chmm, chron_hyper_chmm)

hyper_enrich_plot<- hyper_enrich %>%
  ggplot(aes(x=annotation, y=estimate, colour = model, alpha=padj<0.05)) +
  geom_point(aes(shape = model),size = 1) +
  geom_line(aes(group = model)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  #scale_colour_manual(values = c("purple1", "purple4"), name = "") +
  geom_errorbar(ymin = hyper_enrich$conf.low, ymax = hyper_enrich$conf.high, width = 0.3) +
  scale_colour_manual(values = c("steelblue2", "green4", "purple"), name = "") +
  theme_classic() +
  theme(legend.position = "none",
    panel.background = element_rect(colour = "black", linewidth=0.5),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "pt")) +
  ylab("Odds Ratio") +
  xlab("Annotation") +
  scale_y_continuous(breaks = seq(0,11,2), limits = c(0, 11))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/hyper_enrich.svg", 
       hyper_enrich_plot, 
       height = 40, width = 105, units = "mm")



#Genes
mm_mart<- useEnsembl(biomart="genes", dataset="mmulatta_gene_ensembl")
mm_genes<- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
             mart = mm_mart)
colnames(mm_genes)<- c("anno", "gene_name")

pqlseq_prom3<- inner_join(pqlseq_prom, mm_genes, by = "anno")
pqlseq_prom3<- pqlseq_prom3 %>%
  filter(!gene_name == "") %>%
  filter(!grepl("RNA", gene_name)) %>%
  filter(!grepl("Metazoa", gene_name))
pqlseq_prom3$unique_cpg<- paste(pqlseq_prom3$chr, pqlseq_prom3$cpg_loc, sep="_")

#regions were imported into the environment using the run pqlseq script
#one DNAm region overlaps two genes (MAP2K1 and TIPIN) so need to figure that out
#Maybe run pqlseq using the proms as defined as 2kb upstream instead 
regions_cov<- do.call(rbind, regions_cov)
regions_m<- do.call(rbind, regions_m)
perc_meth<- regions_m/regions_cov
rownames(perc_meth)<- str_split_i(rownames(perc_meth), "\\.", 4)
regs<- c("7_42975452_42975501", "3_91395305_91395354", "7_42975452_42975501")
perc_meth<- perc_meth[rownames(perc_meth) %in% regs,]

regions_cov[rownames(regions_cov) == ""]

#Enrichment for individual models
enrichment<- function(model_df, model_type, var_type){
  
  df_list<- list()
  
  for(i in unique(model_df$anno)) {
    
    df2<- model_df %>%
      distinct(unique_cpg, .keep_all = T)
    
    df3<- model_df %>%
      distinct(unique_cpg, .keep_all = T) %>%
      filter(!unique_cpg %in% df2$unique_cpg)
    
      a<- nrow(df2[df2[[model_type]] == var_type & df2$anno == i,])
      b<- nrow(df2[df2[[model_type]] == var_type & df2$anno != i,])
      
      c<- nrow(df2[df2[[model_type]] != var_type & df2$anno == i,])
      d<- nrow(df2[df2[[model_type]] != var_type & df2$anno != i,])
      
      #Generate contingency table
      c_table<- data.frame("x" = c(a, b),
                           "y" = c(c, d),
                           row.names = c(paste(i, "Y", sep=""), paste(i, "N", sep="")))
      
      colnames(c_table) = c(paste("Is", var_type), paste("Is", "NOT", var_type))
      
      if (all.equal(sum(c_table), length(unique(model_df$unique_cpg)))){
        print(paste("Contingency table sum for", i, "matches unique cpg_loc length"))
     
      df_list[[length(df_list)+1]] = c_table
      
      print(c_table)
    }
  }
  #name table list
  names(df_list)<- unique(model_df$anno)
  
  #Fisher test for each table and tidy with broom
  ft<- lapply(df_list, fisher.test)
  ft<- lapply(ft, broom::tidy)
  
  ft<- do.call(rbind, ft)
  ft<- ft %>%
    mutate(annotation = rownames(ft))
  
  #FDR p-val adjustment
  ft<- ft %>%
    mutate(padj = p.adjust(p.value)) %>%
    mutate_at(vars(annotation), as.factor)
  
  #Log estimates and CIs
  ft<- ft %>%
    mutate(log_or = log(estimate),
           log_ci.lo = log(conf.low),
           log_ci.hi = log(conf.high))
  
  annos_order<- str_sort(ft$annotation, numeric = TRUE)
  
  ft$anno_source<- "Repeat Elements"
  ft$anno_source[ft$annotation %in% chmm_ordered]<- "Chromatin States"

  ft<- ft %>%
    arrange(anno_source, annotation)
  
  #Rearrange factors to sort by type then log_or
  ft$annotation<- factor(ft$annotation, levels = rev(annos_order))
  
  ft$type<- var_type
  
  return(list(ft=ft, c_tables=df_list))
}

cross_hypo<- enrichment(pqlseq_anno, "cross_signif",  "Age-Hypomethylated")
cross_hyper<- enrichment(pqlseq_anno, "cross_signif",  "Age-Hypermethylated")

cross_enrich<- rbind(cross_hypo$ft, cross_hyper$ft)

#Chronological Age
chron_hypo<- enrichment(pqlseq_anno, "chron_signif", "Age-Hypomethylated")
chron_hyper<- enrichment(pqlseq_anno, "chron_signif", "Age-Hypermethylated")

chron_enrich<- rbind(chron_hypo, chron_hyper)


#Eq3 Fisher Enrichment
eq3_hypo<- enrichment(pqlseq_anno, "eq3_age_signif", "Age-Hypomethylated")
eq3_hyper<- enrichment(pqlseq_anno, "eq3_age_signif",  "Age-Hypermethylated")

eq3_enrich<- rbind(eq3_hypo, eq3_hyper)

#Within age fisher enrichment
within_hypo<- enrichment(pqlseq_anno, "eq2_w_signif", "Age-Hypomethylated")
within_hyper<- enrichment(pqlseq_anno, "eq2_w_signif", "Age-Hypermethylated")

within_enrich<- rbind(within_hypo, within_hyper)

#Between age fisher enrichment
between_hypo<- enrichment(pqlseq_anno, "eq2_m_signif", "Age-Hypomethylated")
between_hyper<- enrichment(pqlseq_anno, "eq2_m_signif", "Age-Hypermethylated")

between_enrich<- rbind(between_hypo, between_hyper)

#Combine all enrichment dfs
within_enrich$model<- "within"
chron_enrich$model<- "chron"
cross_enrich$model<- "cross"
eq3_enrich$model<- "eq3"
between_enrich$model<- "between"

all_enrich<- rbind(cross_enrich, chron_enrich, within_enrich, eq3_enrich)

chmm_enrich<- all_enrich %>%
  filter(anno_source == "Chromatin States")

chmm_enrich$annotation<- factor(chmm_enrich$annotation, levels = rev(chmm_ordered))

chmm_enrich %>%
  ggplot(aes(x=annotation, y=estimate, colour = model, alpha=padj<0.05)) +
  geom_point(aes(shape = type), size = 2) +
  geom_line(aes(group = model)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_colour_manual(values = c("darkgoldenrod1", "steelblue1", "green4", "purple"), name = "") +
  geom_errorbar(ymin = chmm_enrich$conf.low, ymax = chmm_enrich$conf.high, width = 0.3) +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  #ylim(c(-6, 4)) +
  ylab("Odds Ratio") +
  xlab("Annotation") +
  facet_wrap(vars(type), nrow=2)

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/enrich_hypo.svg", 
       enrich_hypo, 
       height = 50, width = 105, units = "mm")

re_enrich<- all_enrich %>%
  filter(anno_source == "Repeat Elements")

re_enrich$annotation<- factor(re_enrich$annotation, levels = re_ordered)

re_enrich %>%
  ggplot(aes(x=annotation, y=estimate, colour = model)) +
  geom_point(aes(alpha=padj<0.05)) +
  #geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  geom_line(aes(group = model)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_colour_manual(values = c("steelblue1", "darkgoldenrod1", "purple","green4"), name = "") +
  geom_errorbar(ymin = re_enrich$conf.low, ymax = re_enrich$conf.high, width = 0.3) +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  #theme(legend.position = "top") +
  #ylim(c(-6, 4)) +
  ylab("Odds Ratio") +
  xlab("Annotation") +
  facet_wrap(vars(type), ncol = 1, scales = "free_y")

within_enrich_chmm %>%
  ggplot(aes(x=annotation, y=estimate, fill = type)) +
  #geom_point(aes(alpha=padj<0.05)) +
  geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  #geom_line(aes(group = type)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("green1", "green4"), name = "") +
  #geom_errorbar(ymin = within_enrich$log_ci.lo, ymax = within_enrich$log_ci.hi, width = 0.3, position = position_dodge(0.5)) +
  theme_classic(base_size = 24) +
  theme(legend.position = "top") +
  #ylim(c(-6, 4)) +
  ylab("Odds Ratio") +
  xlab("Annotation") +
  coord_flip()

within<- pqlseq_anno %>%
  select(unique_cpg, beta_eq2_w_age) %>%
  arrange(desc(beta_eq2_w_age))

within2<- within$beta_eq2_w_age
names(within2) = within$unique_cpg

#Enrichment for Hallmark set
within_gsea<- fgseaSimple(pathways = func_region_list, 
                           stats = within2,
                           minSize = 5,
                           maxSize = length(within2) - 1,
                           nperm = 100)

within_gsea<- within_gsea %>%
  arrange(pathway, NES)

annos_sorted<- str_sort(within_gsea$pathway, numeric = T)
annos_sorted<- annos_sorted[-20]
annos_sorted<- c("Promoter", annos_sorted)
annos_sorted<- factor(annos_sorted, levels = annos_sorted)

within_gsea$pathway<- factor(within_gsea$pathway, levels = rev(annos_sorted))

within_gsea %>%
  ggplot(aes(x=pathway, y=NES, colour = NES < 0)) +
  geom_point(aes(alpha=padj<0.05)) +
  #geom_line(aes(group = type)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #geom_errorbar(ymin = within_gsea$log_ci.lo, ymax = within_gsea$log_ci.hi, width = 0.3, position = position_dodge(0.5)) +
  theme_classic(base_size =20) +
  #theme(legend.position = "none") +
  #ylim(c(-1, 7)) +
  ylab("NES") +
  xlab("Annotation") +
  coord_flip()


#### Enrichment for regions only significant in eq2/eq3 ####
enrichment<- function(model_df, model_type, var_type){
  
  df_list<- list()
  
  for(i in unique(model_df$anno)) {
    
    df2<- model_df %>%
      distinct(unique_cpg, .keep_all = T)
    
    df3<- model_df %>%
      distinct(unique_cpg, .keep_all = T) %>%
      filter(!unique_cpg %in% df2$unique_cpg)
    
    a<- nrow(df2[df2[[model_type]] == var_type & df2$anno == i,])
    b<- nrow(df2[df2[[model_type]] == var_type & df2$anno != i,])
    
    c<- nrow(df2[df2[[model_type]] != var_type & df2$anno == i,])
    d<- nrow(df2[df2[[model_type]] != var_type & df2$anno != i,])
    
    #Generate contingency table
    c_table<- data.frame("x" = c(a, b),
                         "y" = c(c, d),
                         row.names = c(paste(i, "Y", sep=""), paste(i, "N", sep="")))
    
    colnames(c_table) = c(paste("Is", var_type), paste("Is", "NOT", var_type))
    
    if (all.equal(sum(c_table), length(unique(model_df$unique_cpg)))){
      print(paste("Contingency table sum for", i, "matches unique cpg_loc length"))
      
      df_list[[length(df_list)+1]] = c_table
      
      print(c_table)
    }
  }
  #name table list
  names(df_list)<- unique(model_df$anno)
  
  #Fisher test for each table and tidy with broom
  ft<- lapply(df_list, fisher.test)
  ft<- lapply(ft, broom::tidy)
  
  ft<- do.call(rbind, ft)
  ft<- ft %>%
    mutate(annotation = rownames(ft))
  
  #FDR p-val adjustment
  ft<- ft %>%
    mutate(padj = p.adjust(p.value)) %>%
    mutate_at(vars(annotation), as.factor)
  
  #Log estimates and CIs
  ft<- ft %>%
    mutate(log_or = log(estimate),
           log_ci.lo = log(conf.low),
           log_ci.hi = log(conf.high))
  
  annos_order<- str_sort(ft$annotation, numeric = TRUE)
  
  ft$anno_source<- "Repeat Elements"
  ft$anno_source[ft$annotation %in% chmm_ordered]<- "Chromatin States"
  
  ft<- ft %>%
    arrange(anno_source, annotation)
  
  #Rearrange factors to sort by type then log_or
  ft$annotation<- factor(ft$annotation, levels = rev(annos_order))
  
  ft$type<- var_type
  
  return(ft)
}

cross_hypo<- enrichment(pqlseq_anno, "cross_signif",  "Age-Hypomethylated")
cross_hyper<- enrichment(pqlseq_anno, "cross_signif",  "Age-Hypermethylated")

#Save workspace image
save.image("/scratch/ckelsey4/Cayo_meth/cross_within_compare.RData")
