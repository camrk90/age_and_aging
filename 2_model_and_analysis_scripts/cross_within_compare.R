library(tidyverse)
library(ggplot2)
library(ggcorrplot)
library(ggpubr)
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
library(broom)

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
  #model$outcome<- str_split_i(model$outcome, "\\.", 3)
  #model$outcome2<- model$outcome
  
  #Separate region coordinates into start and end, delete the chr col, and move region col to front
  #model<- model %>% 
    #separate_wider_delim(outcome2, names=c("chr", "chromStart", "chromEnd"), delim = "_") %>%
    #relocate(c(chr, chromStart, chromEnd), .after = outcome)
  
  #Add length col and filter by length
 # model<- model %>%
   # mutate(length = 1+(as.numeric(chromEnd) - as.numeric(chromStart))) %>%
   # relocate(length, .after=outcome)
}

######################################
###        Import Metadata         ###
######################################
#Metadata
long_data<- read.table("/scratch/ckelsey4/Cayo_meth/long_data_adjusted.txt")

blood_metadata<- read.table("/scratch/ckelsey4/Cayo_meth/blood_metadata_full.txt")

long_ids<- unique(long_data$monkey_id)

overlap_lids<- blood_metadata[blood_metadata$monkey_id %in% long_ids,]

overlap_lids<- overlap_lids %>%
  group_by(monkey_id) %>%
  sample_n(1)

lids_to_remove<- long_data[!long_data$lid_pid %in% overlap_lids$lid_pid,]
blood_metadata<- blood_metadata[!blood_metadata$lid_pid %in% lids_to_remove$lid_pid,]

blood_metadata<- blood_metadata %>%
  filter(age_at_sampling > 1)

#Plot samples collected during each month and each trapping period
long_data<- long_data %>% 
  mutate(year = year(processing_timestamp), month = month(processing_timestamp, label = T))
long_data %>% ggplot(aes(month)) + geom_bar() + facet_wrap(vars(year))

#Import Model data--------------------------------------------------------------
#Import longitudinal pqlseq files
setwd('/scratch/ckelsey4/Cayo_meth/glmer_model_compare')

#Chronological Age
chron_age_files<- 'wb_pqlseq2_agechron'
chron_age_rrbs<- import_pqlseq(chron_age_files, y = 4)
chron_age_rrbs<- chron_age_rrbs[,c(1, 3:8)]

#Age Within
eq2_age_w_files<- 'wb_pqlseq2_within_age'
eq2_age_w_rrbs<- import_pqlseq(eq2_age_w_files, y = 5)
eq2_age_w_rrbs<- eq2_age_w_rrbs[,c(1, 3:8)]

#Mean Age
eq2_age_m_files<- 'wb_pqlseq2_mean_age'
eq2_age_m_rrbs<- import_pqlseq(eq2_age_m_files, y = 5)
eq2_age_m_rrbs<- eq2_age_m_rrbs[,c(1, 3:8)]

#Eq. 3
eq3_age_files<- 'wb_pqlseq2_eq3'
eq3_age_rrbs<- import_pqlseq(eq3_age_files, y = 4)
eq3_age_rrbs<- eq3_age_rrbs[,c(1, 3:8)]

#Eq. 3 Mean Age
eq3_age_m_files<- 'wb_pqlseq2_eq3_m'
eq3_age_m_rrbs<- import_pqlseq(eq3_age_m_files, y = 5)
eq3_age_m_rrbs<- eq3_age_m_rrbs[,c(1, 3:8)]

#Join model dataframes
age_full<- inner_join(chron_age_rrbs, eq2_age_w_rrbs, by = "outcome")
age_full<- inner_join(age_full, eq2_age_m_rrbs, by = "outcome")
age_full<- inner_join(age_full, eq3_age_rrbs, by = "outcome")
age_full<- inner_join(age_full, eq3_age_m_rrbs, by = "outcome")

#Cross Sectional Models
setwd('/scratch/ckelsey4/Cayo_meth/cross_models')

long_cross_files<- 'cs_pqlseq2_age_'
long_cross_pqlseq<- import_pqlseq(long_cross_files, y = 4)

long_cross_pqlseq<- long_cross_pqlseq %>% 
  dplyr::select(-c(h2_cross, sigma2_cross)) 

age_full<- inner_join(age_full, long_cross_pqlseq, by = "outcome")

age_full$outcome<- str_split_i(age_full$outcome, "\\.", 3)

age_full<- age_full %>%
  mutate(chr = str_split_i(outcome, "_", 1),
         chromStart = str_split_i(outcome, "_", 2),
         chromEnd = str_split_i(outcome, "_", 3),
           mutate(across(2:38, as.numeric))) %>%
  relocate(chr, chromStart, chromEnd, .before = outcome)

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

#Twist--------------------------------------------------------------------------
#Chronological Age
chron_age_files<- 'twist_pqlseq2_agechron'
chron_age_twist<- import_pqlseq(chron_age_files, y = 4)

#Age Within
eq2_age_w_files<- 'twist_pqlseq2_within_age'
eq2_age_w_twist<- import_pqlseq(eq2_age_w_files, y = 5)
eq2_age_w_twist<- eq2_age_w_twist[,c(1, 3:10)]

#Mean Age
eq2_age_m_files<- 'twist_pqlseq2_mean_age'
eq2_age_m_twist<- import_pqlseq(eq2_age_m_files, y = 5)
eq2_age_m_twist<- eq2_age_m_twist[,c(1, 3:10)]

#Eq. 3
eq3_age_files<- 'twist_pqlseq2_eq3'
eq3_age_twist<- import_pqlseq(eq3_age_files, y = 4)
eq3_age_twist<- eq3_age_twist[,c(1, 3:10)]

#Eq. 3 Mean Age
eq3_age_m_files<- 'twist_pqlseq2_eq3_m'
eq3_age_m_twist<- import_pqlseq(eq3_age_m_files, y = 5)
eq3_age_m_twist<- eq3_age_m_twist[,c(1, 3:10)]

#Join model dataframes
twist_full<- inner_join(chron_age_twist, eq2_age_w_twist, by = "outcome")
twist_full<- inner_join(twist_full, eq2_age_m_twist, by = "outcome")
twist_full<- inner_join(twist_full, eq3_age_twist, by = "outcome")
twist_full<- inner_join(twist_full, eq3_age_m_twist, by = "outcome")

#Plot Age and Sample Distributions----------------------------------------------
#Cross-age age distribution
cs_samples<- blood_metadata %>%
  ggplot(aes(age_at_sampling, fill=individual_sex)) +
  geom_histogram(position = "dodge", colour = "black") +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  scale_y_continuous(breaks = seq(0, 30, 5)) +
  coord_cartesian(ylim = c(0, 30)) +
  scale_fill_manual(values = c("red3", "pink2"), name = "Sex") +
  theme_classic(base_size=6) +
  theme(legend.position = "bottom",
        legend.key.height = unit(1, "mm"),
        legend.title = element_blank(),
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt")) +
  ylab("Count") +
  xlab("Age")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/cs_samples_dnam.svg", 
       cs_samples, 
       height = 45, width = 65, units = "mm")

#Longitudinal data distribution
long_data<- long_data %>%
  group_by(monkey_id) %>%
  mutate(min_age = min(age_at_sampling))
long_data$age_at_sampling<- round(long_data$age_at_sampling, 0)

samples_dist_dnam<- long_data %>%
  ggplot(aes(x=age_at_sampling, y=reorder(monkey_id, min_age), colour=individual_sex)) +
  geom_path(linewidth = 0.5) +
  geom_point(colour="black", size = 0.25) +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  coord_cartesian(xlim = c(0, 30)) +
  scale_colour_manual(values = c("red3", "pink2"), name = "Sex") +
  theme_classic(base_size = 6) +
  theme(legend.position = "none", 
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin = margin(1, 1, 1, 1, "pt")) +
  xlab("Age") +
  ylab("Individual")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/samples_dist_dnam.svg", 
       samples_dist_dnam, 
       height = 85, width = 55, units = "mm")

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
  ylab("Count") +
  xlab("N Samples") +
  theme_classic(base_size=6) +
  theme(panel.background = element_rect(colour = "black", linewidth=0.5),
        axis.line = element_line(colour = "black", linewidth = 0.25),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        legend.position = "none")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/n_samples_dnam.svg", 
       n_samples, 
       height = 20, width = 30, units = "mm")


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
#P-values
twist_full %>%
  mutate(eq2_signif = ifelse(pvalue_eq2_m_age < .05, "Y", "N")) %>%
  ggplot(aes(beta_eq2_m_age, fill = eq2_signif)) +
  geom_histogram(bins=100, position = "identity", alpha = 0.8) +
  theme_classic(base_size = 6) +
  theme(legend.key.width = unit(3, 'mm'), 
        legend.key.height = unit(1, 'mm'),
        legend.position = "top",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  xlab(expression(beta["Eq.2"] ~ "Between-Age")) +
  ylab("Count")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq2_m_hist_dnam.svg", 
       eq2_m_hist, 
       height = 50, width = 50, units = "mm")

twist_full %>%
  mutate(eq3_signif = ifelse(pvalue_eq3_age_m < .05, "Y", "N")) %>%
  ggplot(aes(beta_eq3_age_m, fill = eq3_signif)) +
  geom_histogram(bins=100, alpha = 0.8) +
  theme_classic(base_size = 6) +
  theme(legend.key.width = unit(3, 'mm'), 
        legend.key.height = unit(1, 'mm'),
        legend.position = "top",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  xlab(expression(beta["Eq.3"] ~ "Between-Age")) +
  ylab("Count")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq3_m_hist_dnam.svg", 
       eq3_m_hist, 
       height = 50, width = 50, units = "mm")

twist_full %>%
  ggplot(aes(beta_eq2_w_age, -log10(fdr_eq2_m_age), colour = -log10(fdr_eq2_m_age) > -log10(0.05))) +
  geom_point(alpha = 0.5, size = 0.05) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  #scale_colour_manual(values = c('steelblue4', 'steelblue1')) +
  theme_classic(base_size = 7) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1) +
  xlab("Beta Eq.2 Between-Age") +
  ylab("-log10(FDR)")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq2_m_volcano_dnam.svg", 
       eq2_m_volcano, 
       height = 50, width = 50, units = "mm")

eq3_m_volcano<- age_trunc %>%
  ggplot(aes(beta_eq3_age_m, -log10(fdr_eq3_age_m), colour = -log10(fdr_eq3_age_m) < -log10(0.05))) +
  geom_point(alpha = 0.5, size = 0.05) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic(base_size = 7) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1) +
  xlab(paste(expression(beta["Eq.3"]), "Between-Age", sep = " ")) +
  ylab("-log10(FDR)")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq3_m_volcano_dnam.svg", 
       eq3_m_volcano, 
       height = 50, width = 50, units = "mm")

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
  filter(model != "fdr_eq3_age") %>%
  ggplot(aes(fdr_threshold, count, colour = model)) +
  geom_point() +
  geom_path() +
  scale_colour_manual(values = c("darkgoldenrod2","steelblue1", "green4"), name = "") +
  theme_classic(base_size = 6) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5))

beta_thresholds<- c(0.25, 0.20, 0.15, 0.10, 0.05, 0.025, 0.01, 0.001)

res <- map_dfr(thresholds, function(x) {
  
  df <- age_trunc %>%
    filter(fdr_eq3_age_m < x)
  
  test <- tidy(cor.test(df$beta_chron_age, df$beta_eq2_w_age))
  
  test$threshold <- x
  
  test
})

res %>%
  ggplot(aes(threshold*-1, estimate)) +
  geom_point() +
  geom_path() +
  geom_text(label = res$parameter, hjust = 1.5, vjust = 1, size =3) +
  geom_errorbar(ymin = res$conf.low, ymax = res$conf.high, width = 0.01) +
  theme_classic(base_size = 12) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  scale_x_continuous(breaks = seq(0, -0.30, -0.05)) +
  ylab("Pearson Correlation") +
  xlab(expression(beta["Eq.3"] ~ "FDR Threshold"))

#### Coefficient comparisons ---------------------------------------------------
#Generate plot function
compare_plot<- function(df, fdr1, fdr2, var1, var2, plot_type) {
  
  v1<-deparse(substitute(var1))
  v2<-deparse(substitute(var2))
  
  df<- df %>%
    filter({{fdr1}} < .05 | {{fdr2}} < .05) %>%
    mutate(diff = abs({{var2}}) - abs({{var1}}),
           ratio = abs({{var2}})/abs({{var1}}))
  
  print(paste("The median", v2, "/", v1, "=", median(df$ratio), sep = " "))
  print(paste("The median", v2, "-", v1, "=", median(df$diff), sep = " "))
  correlation<- cor(df %>% pull({{var1}}), df %>% pull({{var2}}))
  
  print(paste("The correlation between", 
              v1, "and", 
              v2, "=", 
              correlation,
              sep = " "))
  
  print(paste("Total # of shared significant genes = ", nrow(df)))
  print(paste("Total # of genes where", v1, ">", v2, "=", nrow(df[df$ratio < 1,]), sep = " "))
  
  if (plot_type == "scatter"){
    
    df %>%
      ggplot(aes({{var1}}, {{var2}}, colour = diff)) +
      geom_point(size = 1, alpha = 0.8) +
      geom_abline() +
      geom_smooth(method = "lm", linewidth = 0.5) +
      geom_vline(xintercept=0, linetype="dashed") +
      geom_hline(yintercept=0, linetype="dashed") +
      theme_classic(base_size = 6) +
      theme(legend.key.width = unit(5, 'mm'), 
            legend.key.height = unit(1, 'mm'),
            legend.position = "top") +
      theme(panel.background = element_rect(colour = "black", linewidth=1),
            axis.line = element_line(colour = "black", linewidth = 0.5),
            plot.margin = margin(1, 1, 1, 1, "pt"),
            aspect.ratio = 1,
            panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
            panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
      stat_cor()
    
  } else if (plot_type == "hist") {
    
    df %>%
      ggplot(aes(diff, fill = after_stat(x))) +
      geom_histogram(bins = 50) +
      geom_vline(xintercept=0, linetype="dashed") +
      geom_vline(xintercept=median(df$diff), linetype="dashed", colour = 'red') +
      theme_classic(base_size = 6) +
      theme(legend.key.width = unit(5, 'mm'), 
            legend.key.height = unit(1, 'mm'),
            legend.position = "top") +
      theme(panel.background = element_rect(colour = "black", linewidth=1),
            axis.line = element_line(colour = "black", linewidth = 0.5),
            plot.margin = margin(1, 1, 1, 1, "pt"),
            aspect.ratio = 1,
            panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
            panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) 
    
  }
}

#CS vs Eq.1
cor.test(age_trunc$beta_cross, age_trunc$beta_chron_age)
cs_chron_scatter<- compare_plot(age_trunc, fdr_chron_age, fdr_cross, 
                                beta_chron_age, beta_cross, "scatter")

cs_chron_scatter<- cs_chron_scatter +
  scale_color_gradient2(low = "steelblue2", mid = "grey70", high = "darkgoldenrod2", midpoint = 0, name = "") +
  xlim(-0.1, 0.1) +
  ylim(-0.1, 0.1) +
  xlab(expression(beta["Eq.1"])) +
  ylab(expression(beta["CS"]))
cs_chron_scatter

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/cs_chron_scatter_dnam.svg", 
       cs_chron_scatter, 
       height = 50, width = 50, units = "mm")

cs_chron_hist<- compare_plot(age_trunc, fdr_chron_age, fdr_cross, 
                             beta_chron_age, beta_cross, "hist")

cs_chron_hist<- cs_chron_hist +
  scale_fill_gradient2(low = "steelblue2", mid = "grey70", high = "darkgoldenrod2", midpoint = 0, name = "") +
  xlab(expression(abs(beta["Eq.1"]) - abs(beta["CS"]))) +
  ylab("Count") +
  scale_x_continuous(breaks = seq(-0.05, 0.15, 0.05), limits = c(-0.06, 0.11))
cs_chron_hist

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/cs_chron_hist_dnam.svg", 
       cs_chron_hist, 
       height = 50, width = 50, units = "mm")

#Eq.1 vs Eq.2
cor.test(age_trunc$beta_eq2_w_age, age_trunc$beta_chron_age)
within_chron_plot<- compare_plot(age_trunc, fdr_chron_age, fdr_eq2_w_age,
                                 beta_chron_age, beta_eq2_w_age, "scatter")
within_chron_plot<- within_chron_plot + 
  scale_color_gradient2(low = "steelblue2", mid = "grey70", high = "green4", midpoint = 0, name = "") +
  xlim(-0.2, 0.2) +
  ylim(-0.2, 0.2) +
  xlab("Beta Eq.1") +
  ylab("Beta Eq.2")
within_chron_plot

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/within_chron_scatter.svg", 
       within_chron_plot, 
       height = 50, width = 50, units = "mm")

within_chron_hist<- compare_plot(age_trunc, fdr_chron_age, fdr_eq2_w_age,
                                beta_chron_age, beta_eq2_w_age, "hist")
within_chron_hist<- within_chron_hist + 
  scale_fill_gradient2(low = "steelblue2", mid = "grey70", high = "green4", midpoint = 0, name = "") +
  xlab("Abs(Eq.2) - Abs(Eq.1)") +
  ylab("Count") +
  xlim(-0.1, 0.2)
within_chron_hist

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/within_chron_hist.svg", 
       within_chron_hist, 
       height = 50, width = 50, units = "mm")

#Eq.1 vs Eq.3
eq3_chron_plot<- compare_plot(age_trunc, fdr_chron_age, fdr_eq3_age,
                              beta_chron_age, beta_eq3_age,"scatter")
eq3_chron_plot<- within_chron_plot + 
  scale_color_gradient2(low = "steelblue2", mid = "grey70", high = "purple", midpoint = 0, name = "") +
  xlim(-0.2, 0.2) +
  ylim(-0.2, 0.2) +
  xlab(expression(beta["Eq.1"])) +
  ylab(expression(beta["Eq.3"]))
eq3_chron_plot

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq3_chron_scatter.svg", 
       eq3_chron_plot, 
       height = 45, width = 45, units = "mm")

eq3_chron_hist<- compare_plot(age_trunc, fdr_chron_age, fdr_eq3_age,
                              beta_chron_age, beta_eq3_age,
                              "purple", "steelblue2", "hist")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq3_chron_hist.svg", 
       eq3_chron_hist, 
       height = 45, width = 45, units = "mm")

#Eq.2 vs Eq.3
eq2_eq3<- compare_plot(age_trunc, fdr_eq2_w_age, fdr_eq3_age,
             beta_eq2_w_age, beta_eq3_age, "scatter")

eq2_eq3<- eq2_eq3 +
  scale_color_gradient2(low = "purple", mid = "grey70", high = "green4", midpoint = 0, name = "") +
  theme(legend.position = "none") +
  xlim(-0.2, 0.2) +
  ylim(-0.2, 0.2) +
  xlab(expression(beta["Eq.2"])) +
  ylab(expression(beta["Eq.3"]))
eq2_eq3

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq2_eq3_scatterplot_dnam.svg", 
       eq2_eq3, 
       height = 45, width = 45, units = "mm")

#Counts and Beta Distributions
count_signif_regions<- function(x) {
  
  df<- x %>%
    dplyr::select(starts_with("fdr_"))
  vars<- gsub("fdr_", "", colnames(df))
  
  counts <- colSums(df < 0.05, na.rm = TRUE)
  
  counts<- data.frame(predictor = vars,
                      count = counts)
  
  counts<- counts %>%
    arrange(counts) %>%
    mutate(predictor = factor(predictor, levels = predictor))
  
  counts<- counts %>%
    mutate(perc_signif = count/nrow(df))
  
  counts_plot<- counts %>%
    ggplot(aes(predictor, count, fill = predictor)) +
    geom_bar(stat = 'identity') +
    geom_text(label=counts$count, vjust=-0.25, size = 1.5) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          panel.background = element_rect(colour = "black", linewidth=1),
          axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.title.x = element_blank(),
          plot.margin = margin(1, 1, 1, 1, "pt")) +
    xlab("Predictor") +
    ylab("Count")
  
  return(list(plot = counts_plot, df = counts))
  
}

counts<- count_signif_regions(age_trunc)


signif_counts<-counts %>%
  ggplot(aes(predictor, count, fill = predictor)) +
  geom_bar(stat = 'identity') +
  geom_text(label=counts[,2], vjust=-0.25, size = 1.5) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.title.x = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "pt")) +
  scale_y_continuous(breaks = seq(0, 25000, 5000)) +
  xlab("Predictor") +
  ylab("Count") +
  scale_fill_manual(values = c("purple", 'green4', "steelblue2", "darkgoldenrod2", "grey30"))
signif_counts
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/signif_counts.svg", 
       signif_counts, 
       height = 45, width = 65, units = "mm")

#Distribution of effect sizes for each variable
beta_dist<- age_full %>%
  dplyr::select(c(beta_eq3_age, beta_eq2_w_age, beta_eq2_m_age, beta_chron_age, beta_cross)) %>%
  pivot_longer(cols = c(beta_eq3_age, beta_eq2_w_age, beta_eq2_m_age, beta_chron_age, beta_cross),
               values_to = 'beta',
               names_to = 'var') %>%
  mutate(var = factor(var, levels = rev(c("beta_eq2_m_age", "beta_chron_age", "beta_cross",
                                      "beta_eq2_w_age", "beta_eq3_age")))) %>%
  ggplot(aes(var, beta, fill=var)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.25) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_fill_manual(values = c("purple","green4",'darkgoldenrod2','steelblue2','grey30')) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        axis.title.x = element_blank()) +
  scale_x_discrete(labels = c("Eq.3 W.","Eq.2 W.","CS","Eq.1","Eq.2 B.")) +
  ylab("Beta")
beta_dist
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/beta_dist_dnam.svg", 
       beta_dist, 
       height = 45, width = 65, units = "mm")

## Significant regions Venn diagram
cross.age<- age_trunc$outcome[age_trunc$fdr_cross< 0.05]
chron.age<- age_trunc$outcome[age_trunc$fdr_chron_age< 0.05]
age.w<- age_trunc$outcome[age_trunc$fdr_eq2_w_age < 0.05]
eq2.btwn<- age_trunc$outcome[age_trunc$fdr_eq2_m_age < 0.05]
eq3<- age_trunc$outcome[age_trunc$fdr_eq3_age < 0.05]

venn_all<- list(cross.age, chron.age, age.w, eq3, eq2.btwn)
names(venn_all)<- c("cross.age", "chron.age", "age.w", "eq3", "eq2.btwn")

ggvenn(venn_list2,
       text_size = 8,
       show_percentage = F)

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

annos<- c("TSSs", "Active Tr.", "Enhancers", "Quiescent", "Simple Repeats",
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

rm(chmm_intersect);rm(re_anno);rm(pqlseq_prom2)

#Generate a chmm factor vector that contains promoters 
tx_ordered<- c("0_Promoter", levels(chmm_ordered))
tx_ordered<- str_sort(tx_ordered, numeric = T)
tx_ordered[1]<- str_split_i(tx_ordered[1], "_", 2)
tx_ordered<- factor(tx_ordered, levels = tx_ordered)


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
pqlseq_anno$eq2_m_signif[pqlseq_anno$fdr_eq2_m_age < 0.05 & pqlseq_anno$beta_eq2_m_age < 0]<- "Age-Hypomethylated"
pqlseq_anno$eq2_m_signif[pqlseq_anno$fdr_eq2_m_age < 0.05 & pqlseq_anno$beta_eq2_m_age > 0]<- "Age-Hypermethylated"

pqlseq_anno$eq2_m_signif<- factor(pqlseq_anno$eq2_m_signif, 
                                    levels = c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))

pqlseq_anno$eq3_age_signif<- "Non-Significant"
pqlseq_anno$eq3_age_signif[pqlseq_anno$fdr_eq3_age < 0.05 & pqlseq_anno$beta_eq3_age < 0]<- "Age-Hypomethylated"
pqlseq_anno$eq3_age_signif[pqlseq_anno$fdr_eq3_age < 0.05 & pqlseq_anno$beta_eq3_age > 0]<- "Age-Hypermethylated"

pqlseq_anno$eq3_age_signif<- factor(pqlseq_anno$eq3_age_signif, 
                                   levels = c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))

#save pqlseq_anno df to .rds for gsea script
saveRDS(pqlseq_anno, "/scratch/ckelsey4/Cayo_meth/pqlseq_anno.rds")

#Enrichment Analyses------------------------------------------------------------
##Plot annotation proportions---------------------------------------------------
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

generate_proportion(pqlseq_anno, eq2_m_signif, "green1", "gray90", "green4")

generate_proportion(pqlseq_anno, eq3_age_signif, "purple1", "gray90", "purple4")

##Rank Enrichment----------------------------------------------------------------
gsea<- readRDS('/scratch/ckelsey4/Cayo_meth/within_gsea.rds')

gsea$type<- "Repeat Elements"
gsea$type[gsea$pathway %in% tx_ordered]<- "Chromatin States"

gsea<- gsea %>%
  mutate(pathway = factor(pathway, levels = pathway))

gsea %>%
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

##Fishers Enrichment-------------------------------------------------------------
#Enrichment function for single model
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
  
  ft$source<- "RE"
  ft$source[ft$annotation %in% chmm_ordered | ft$annotation == "Promoter"]<- "CHMM"
  
  ft<- ft %>%
    arrange(source, annotation)
  
  #Rearrange factors to sort by type then log_or
  ft$annotation<- factor(ft$annotation, levels = rev(annos_order))
  
  ft$direction<- var_type
  
  ft<- ft %>% mutate(model = model_type)
  ft_chmm<- ft %>% filter(source == "CHMM") 
  ft_re<- ft %>% filter(source == "RE")
  
  return(list(ft_chmm=ft_chmm, ft_re=ft_re, c_tables=df_list))
}

#Enrichment function for regions only significant in one model
enrichment_uniq<- function(model_df, model_type1, model_type2, var_type1, var_type2){
  
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
  
  ft$direction<- var_type2
  
  ft$source<- "RE"
  ft$source[ft$annotation %in% chmm_ordered | ft$annotation %in% "Promoter"]<- "CHMM"
  
  ft<- ft %>% mutate(model = "Eq.2")
  ft_chmm<- ft %>% filter(source == "CHMM") 
  ft_re<- ft %>% filter(source == "RE")
  
  return(list(ft_chmm=ft_chmm, ft_re=ft_re, c_tables=df_list))
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
within_chron_hypo<- enrichment_uniq(pqlseq_anno, 
                                    "within_chron", "eq2_w_signif",
                                    "Within Age Significant", "Age-Hypomethylated")

chron_hypo<- enrichment(pqlseq_anno, "chron_signif", "Age-Hypomethylated")

cross_hypo<- enrichment(pqlseq_anno, "cross_signif", "Age-Hypomethylated")

#Hyper Enrichment
within_chron_hyper<- enrichment_uniq(pqlseq_anno, 
                                    "within_chron", "eq2_w_signif",
                                    "Within Age Significant", "Age-Hypermethylated")

chron_hyper<- enrichment(pqlseq_anno, "chron_signif", "Age-Hypermethylated")

cross_hyper<- enrichment(pqlseq_anno, "cross_signif", "Age-Hypermethylated")

#Chromatin States
full_enrich_chmm<- rbind(within_chron_hyper[['ft_chmm']], chron_hyper[['ft_chmm']], cross_hyper[['ft_chmm']],
                          within_chron_hypo[['ft_chmm']], chron_hypo[['ft_chmm']], cross_hypo[['ft_chmm']])

full_enrich_chmm$annotation<- factor(full_enrich_chmm$annotation, levels = chmm_ordered)

chmm_enrich_plot<- full_enric_chmm %>%
  ggplot(aes(x=annotation, y=estimate, colour = model, alpha=padj<0.05)) +
  geom_point(size = 1) +
  geom_line(aes(group = model)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_errorbar(ymin = full_enrich$conf.low, ymax = full_enrich$conf.high, width = 0.3) +
  scale_colour_manual(values = c("darkgoldenrod2", "steelblue2", "green4"), name = "") +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=0.5),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("Odds Ratio") +
  xlab("Chromatin State") +
  scale_y_continuous(breaks = seq(0,12,2), limits = c(0, 11)) +
  facet_wrap(vars(direction), ncol = 1)

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/full_enrich_dnam.svg", 
       full_enrich_plot, 
       height = 90, width = 105, units = "mm")

#Repeat Elements
full_enrich_re<- rbind(within_chron_hyper[['ft_re']], chron_hyper[['ft_re']], cross_hyper[['ft_re']],
                         within_chron_hypo[['ft_re']], chron_hypo[['ft_re']], cross_hypo[['ft_re']])

full_enrich_re$annotation<- factor(full_enrich_re$annotation, levels = re_ordered)

re_enrich_plot<- full_enrich_re %>%
  ggplot(aes(x=annotation, y=estimate, colour = model, alpha=padj<0.05)) +
  geom_point(size = 1) +
  geom_line(aes(group = model)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_errorbar(ymin = full_enrich_re$conf.low, ymax = full_enrich_re$conf.high, width = 0.3) +
  scale_colour_manual(values = c("darkgoldenrod2", "steelblue2", "green4"), name = "") +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=0.5),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("Odds Ratio") +
  xlab("Repeat Element") +
  #scale_y_continuous(breaks = seq(0,12,2), limits = c(0, 12)) +
  facet_wrap(vars(direction), ncol = 1, scales = "free_y")
re_enrich_plot

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/re_enrich_dnam.svg", 
       re_enrich_plot, 
       height = 90, width = 55, units = "mm")

#Eq2 Between
eq2_m_hypo<- enrichment(pqlseq_anno, "eq2_m_signif", "Age-Hypomethylated")
eq2_m_hyper<- enrichment(pqlseq_anno, "eq2_m_signif", "Age-Hypermethylated")

eq2_m_full<- rbind(eq2_m_hyper[['ft_chmm']], eq2_m_hypo[['ft_chmm']])
eq2_m_full<- eq2_m_full %>%
  filter(annotation != "Promoter") %>%
  mutate(model = "Eq.2 Btwn")
eq2_m_full<- rbind(eq2_m_full, chron_hypo[['ft_chmm']], chron_hyper[['ft_chmm']])

eq2_m_full$annotation<- factor(eq2_m_full$annotation, levels = chmm_ordered)

eq2_m_enrich<- eq2_m_full %>%
  ggplot(aes(x=annotation, y=estimate, colour = model, alpha=padj<0.05)) +
  geom_point(size = 1) +
  geom_line(aes(group = model)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_errorbar(ymin = eq2_m_full$conf.low, ymax = eq2_m_full$conf.high, width = 0.3) +
  theme_classic(base_size = 6) +
  theme(
        panel.background = element_rect(colour = "black", linewidth=0.5),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("Odds Ratio") +
  xlab("CHMM") +
  #scale_y_continuous(breaks = seq(0,12,2), limits = c(0, 12)) +
  facet_wrap(vars(direction), ncol = 1, scales = "free_y")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq2m_enrich_dnam.svg", 
       eq2_m_enrich, 
       height = 90, width = 105, units = "mm")

#Save workspace image
save.image("/scratch/ckelsey4/Cayo_meth/cross_within_compare.RData")
