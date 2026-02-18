library(tidyverse)
library(ggplot2)
library(ggcorrplot)
library(ggvenn)
library(ggeffects)
#library(variancePartition)
library(lme4)
library(fgsea)
#library(BiocParallel)
library(UpSetR)

#register(MulticoreParam())

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
  
  #Filter for true convergences
  model<- model %>%
    filter(converged == "TRUE")
  
  #Generate df of adjusted pvalues
  model_fdr<- p.adjust(model$pvalue, method = "fdr")
  
  #Bind padj cols to model df and relocate
  model<- cbind(model, model_fdr)
  model<- model %>%
    dplyr::rename(fdr = model_fdr) %>%
    relocate(fdr, .after = pvalue) %>%
    dplyr::select(-c(elapsed_time, converged, h2, sigma2))
}

######################################
###        Import Metadata         ###
######################################
#Metadata
blood_metadata<- read.table("/scratch/ckelsey4/Cayo_meth/blood_metadata_full.txt")
long_data<- read.table("/scratch/ckelsey4/Cayo_meth/long_data_adjusted.txt")

long_ids<- unique(long_data$monkey_id)

overlap_lids<- blood_metadata[blood_metadata$monkey_id %in% long_ids,]

overlap_lids<- overlap_lids %>%
  group_by(monkey_id) %>%
  sample_n(1)

lids_to_remove<- long_data[!long_data$lid_pid %in% overlap_lids$lid_pid,]
blood_metadata<- blood_metadata[!blood_metadata$lid_pid %in% lids_to_remove$lid_pid,]

blood_metadata<- blood_metadata %>%
  filter(age_at_sampling > 1)

long_data<- long_data %>%
  group_by(monkey_id) %>%
  mutate(n = n()) %>%
  ungroup()

long_data<- long_data %>%
  filter(age_at_sampling > 1) %>%
  #filter(n > 1) %>%
  dplyr::rename(perc_unique = unique) %>%
  drop_na() %>%
  arrange(lid_pid)

rm(lids_to_remove);rm(overlap_lids);rm(long_ids)

######################################
###          Import data           ###
######################################
#Import longitudinal pqlseq files
setwd('/scratch/ckelsey4/Cayo_meth/glmer_model_compare')

#Chronological Age
chron_age_files<- 'wb_pqlseq2_agechron'
chron_age_pqlseq<- import_pqlseq(chron_age_files, y = 4)
colnames(chron_age_pqlseq)<- c("outcome", "length", "chr", "chromStart", "chromEnd", "n",
                               paste(names(chron_age_pqlseq[,7:12]), "chron_age", sep = "_"))

#Age Within
eq2_age_w_files<- 'wb_pqlseq2_within_age'
eq2_age_w_pqlseq<- import_pqlseq(eq2_age_w_files, y = 5)
eq2_age_w_pqlseq<- eq2_age_w_pqlseq[,c(1, 7:12)]
colnames(eq2_age_w_pqlseq)<- c("outcome", paste(names(eq2_age_w_pqlseq[,2:7]), "eq2_w_age", sep = "_"))

#Mean Age
eq2_age_m_files<- 'wb_pqlseq2_mean_age'
eq2_age_m_pqlseq<- import_pqlseq(eq2_age_m_files, y = 5)
eq2_age_m_pqlseq<- eq2_age_m_pqlseq[,c(1, 7:12)]
colnames(eq2_age_m_pqlseq)<- c("outcome", paste(names(eq2_age_m_pqlseq[,2:7]), "eq2_mean_age", sep = "_"))

#Eq. 3
eq3_age_files<- 'wb_pqlseq2_eq3'
eq3_age_pqlseq<- import_pqlseq(eq3_age_files, y = 4)
eq3_age_pqlseq<- eq3_age_pqlseq[,c(1, 7:12)]
colnames(eq3_age_pqlseq)<- c("outcome", paste(names(eq3_age_pqlseq[,2:7]), "eq3_age", sep = "_"))

#Eq. 3 Mean Age
eq3_age_m_files<- 'wb_pqlseq2_eq3_m'
eq3_age_m_pqlseq<- import_pqlseq(eq3_age_m_files, y = 5)
eq3_age_m_pqlseq<- eq3_age_m_pqlseq[,c(1, 7:12)]
colnames(eq3_age_m_pqlseq)<- c("outcome", paste(names(eq3_age_m_pqlseq[,2:7]), "eq3_age_m", sep = "_"))

test_files<- "within_age_no_singles"
test_pqlseq<- import_pqlseq(eq3_age_m_files, y = 5)

test<- as.data.frame(cbind(test_pqlseq$outcome, test_pqlseq$beta, test_pqlseq$fdr, eq2_age_w_pqlseq$beta_eq2_w_age, eq2_age_w_pqlseq$beta_eq2_w_age))
colnames(test)<- c("outcome", "beta_test", "fdr_test", "beta_eq2", "fdr_eq2")

test<- test %>%
  mutate(across(.cols = starts_with("beta_"), .fns = as.numeric)) %>%
  mutate(across(.cols = starts_with("fdr_"), .fns = as.numeric)) %>%
  mutate(diff = abs(beta_eq2) - abs(beta_test))

test %>%
  ggplot(aes(beta_test, beta_eq2, colour = diff)) +
  geom_point(alpha = 0.3)

#Join model dataframes
age_full<- inner_join(chron_age_pqlseq, eq2_age_w_pqlseq, by = "outcome")
age_full<- inner_join(age_full, eq2_age_m_pqlseq, by = "outcome")
age_full<- inner_join(age_full, eq3_age_pqlseq, by = "outcome")
age_full<- inner_join(age_full, eq3_age_m_pqlseq, by = "outcome")

#Cross Sectional Models
setwd('/scratch/ckelsey4/Cayo_meth/cross_models')

long_cross_files<- '_long'
long_cross_pqlseq<- import_pqlseq(long_cross_files, y = 3)

long_outcome<- long_cross_pqlseq %>%
  select(outcome)

long_cross_pqlseq<- long_cross_pqlseq %>%
  select(-c(outcome, length, chr, chromStart, chromEnd, n)) %>%
  mutate_if(is.character, as.numeric)

long_cross_pqlseq<- cbind(long_outcome, long_cross_pqlseq)

colnames(long_cross_pqlseq)<- c("outcome", paste(colnames(long_cross_pqlseq[,2:7]), "cross", sep = "_"))

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
age_full$within_chron<- factor(age_full$within_chron, levels = age_full$within_chron)

#Generate age df subset
age_trunc<- age_full %>%
  select(c(outcome, region_range, chr, 
           beta_cross, fdr_cross, #cross sectional age
           beta_chron_age, fdr_chron_age, #chron_age
           beta_eq2_w_age, fdr_eq2_w_age, #eq2 within age
           beta_eq2_mean_age, fdr_eq2_mean_age, #eq2 mean age
           beta_eq3_age, fdr_eq3_age, #eq3 age
           beta_eq3_age_m, fdr_eq3_age_m)) #eq3 mean age 

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

long_data %>%
  ggplot(aes(x=age_at_sampling, y=reorder(monkey_id, min_age), colour=individual_sex)) +
  geom_path(linewidth = 1.2, alpha = 0.8) +
  geom_point(colour="black", cex = 1) +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  coord_cartesian(xlim = c(0, 30)) +
  scale_colour_manual(values = c("red3", "pink2"), name = "Sex") +
  ylab("Individual") +
  xlab("Age") +
  theme_classic(base_size = 24) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  #theme(legend.key.height= unit(2, 'cm')) +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(colour = "black", linewidth=2))

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
long_data %>%
  ggplot(aes(x=n, fill=as.factor(individual_sex))) +
  geom_bar(colour='black', position = 'dodge') +
  scale_fill_manual(values = c("red3", "pink2"), name = "Sex") +
  ylab("Count") +
  xlab("N Samples") +
  theme_classic(base_size = 24) +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(colour = "black", linewidth=3))

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


#Plot Distributions-------------------------------------------------------------
#P-Values
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

compare_plot<- function(df, fdr1, fdr2, var1, var2, c1, c2, lab1, lab2, plot_type) {
  
  df<- df %>%
    filter({{fdr1}} < .05 & {{fdr2}} < .05) %>%
    mutate(diff = abs({{var2}}) - abs({{var1}}))
  
  #eval(substitute(df_lm<- lm(var1 ~ var2, data=df)))
  
  #print(summary(df_lm))
  
  if (plot_type == "scatter"){
    
    df %>%
      ggplot(aes({{var1}}, {{var2}}, colour = diff)) +
      geom_point(cex = 1,
                 shape = 1,
                 alpha = 0.8) +
      geom_abline() +
      #geom_abline(slope = df_lm[["coefficients"]][[2]], 
      #intercept = df_lm[["coefficients"]][[1]],
      #colour = "red") +
      geom_smooth(method = "lm") +
      geom_vline(xintercept=0, linetype="dashed") +
      geom_hline(yintercept=0, linetype="dashed") +
      scale_color_gradient2(low = c2, mid = "grey70", high = c1, midpoint = 0, name = "") +
      theme_classic(base_size=24) +
      theme(legend.key.width = unit(2, 'cm'), legend.position = "top") +
      theme(panel.background = element_rect(colour = "black", linewidth=3)) +
      theme(aspect.ratio = 1) +
      xlim(-0.20, 0.20) +
      ylim(-0.20, 0.20) +
      xlab(lab1) +
      ylab(lab2)
    
  } else if (plot_type == "hist") {
    
    df %>%
      ggplot(aes(diff, fill = after_stat(x))) +
      geom_histogram(bins = 50, colour="black") +
      geom_vline(xintercept=0, linetype="dashed") +
      scale_fill_gradient2(low = c2, mid = "grey70", high = c1, midpoint = 0, name = "") +
      theme_classic(base_size=24) +
      theme(legend.position = "none") +
      theme(panel.background = element_rect(colour = "black", linewidth=3)) +
      theme(aspect.ratio = 1) +
      xlab(paste(lab2, "-", lab1, sep=" "))
    
  }
}

#### Coefficient comparisons ---------------------------------------------------
compare_plot(age_trunc, fdr_chron_age, fdr_cross, 
             beta_chron_age, beta_cross, 
             "darkgoldenrod2", "steelblue2", 
             "Chron Age", "Cross Age","scatter")

compare_plot(age_full, fdr_chron_age, fdr_cross,
             beta_chron_age, beta_cross,
             "darkgoldenrod2", "steelblue2",
             "Chron Age", "Cross Age", "hist")

compare_plot(age_full, fdr_chron_age, fdr_eq2_w_age,
             beta_chron_age, beta_eq2_w_age,
             "green4", "steelblue2",
             "Chron Age", "Eq2. Within Age", "scatter")

compare_plot(age_full, fdr_cross, fdr_eq2_w_age,
             beta_cross, beta_eq2_w_age,
             "green4", "steelblue2",
             "Chron Age", "Eq2. Within Age", "hist")

compare_plot(age_full, fdr_chron_age, fdr_eq3_age,
             beta_chron_age, beta_eq3_age,
             "purple", "steelblue2", 
             "Chron Age", "Eq3. Age", "scatter")

compare_plot(age_full, fdr_chron_age, fdr_eq3_age,
             beta_chron_age, beta_eq3_age,
             "purple", "steelblue2", 
             "Chron Age", "Eq3. Age", "hist")

compare_plot(age_full, fdr_eq2_w_age, fdr_eq3_age,
             beta_eq2_w_age, beta_eq3_age,
             "purple", "green4",
             "Eq2. Within Age", "Eq3. Age", "scatter")

compare_plot(age_full, fdr_eq2_w_age, fdr_eq3_age,
             beta_eq2_w_age, beta_eq3_age,
             "purple", "green4",
             "Eq2. Within Age", "Eq3. Age", "hist")

age.cross.count<- nrow(age_full[age_full$fdr_cross < 0.05,])
age.chron.count<- nrow(age_full[age_full$fdr_chron_age < 0.05,])
age.w.count<- nrow(age_full[age_full$fdr_eq2_w_age < 0.05,])
age.eq3.count<- nrow(age_full[age_full$fdr_eq3_age < 0.05,])
counts<- data.frame(count = c(age.w.count, age.eq3.count, age.chron.count, age.cross.count),
                    predictor = as.factor(c('Eq2. Within Age', 'Eq3 Age', 'Chron Age', 'Cross Age')))

counts<- counts %>%
  mutate(predictor = as.factor(predictor)) %>%
  mutate(predictor=fct_reorder(predictor, count, .desc=T))

rm(age.w.count);rm(age.cross.count);rm(age.chron.count);rm(age.eq3.count)

counts %>%
  ggplot(aes(predictor, count, fill = predictor)) +
  geom_bar(stat = 'identity', colour="black") +
  geom_text(label=counts$count, vjust=-1, size=5) +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_text(angle = 15, hjust=0.9),
        legend.position = "none") +
  xlab("Predictor") +
  ylab("Count") +
  scale_fill_manual(values = c("purple", 'green4', 'steelblue2', 'darkgoldenrod2'))

#Distribution of effect sizes for each variable
age_full %>%
  dplyr::select(c(beta_eq2_w_age, beta_eq3_age, beta_chron_age, beta_cross)) %>%
  pivot_longer(cols = c(beta_eq2_w_age, beta_eq3_age, beta_chron_age, beta_cross),
               values_to = 'beta',
               names_to = 'var') %>%
  ggplot(aes(beta, fill=var)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_fill_manual(values = c('steelblue2', 'darkgoldenrod2', 'green4', "purple"),
                    labels = c("Chron Age", "Cross Age", "Eq2. Within Age", "Eq3. Age")) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none") +
  xlim(-0.25, 0.25) +
  ylab("Density") +
  xlab("Beta")

## Significant regions Venn diagram
cross.age<- age_trunc$outcome[age_trunc$fdr_cross< 0.05]
chron.age<- age_trunc$outcome[age_trunc$fdr_chron_age< 0.05]
age.w<- age_trunc$outcome[age_trunc$fdr_eq2_w_age < 0.05]
eq3<- age_trunc$outcome[age_trunc$fdr_eq3_age < 0.05]

venn_list1<- list(cross.age, chron.age, age.w)
names(venn_list1)<- c("cross.age", "chron.age", "age.w")

ggvenn(venn_list1,
       text_size = 8,
       show_percentage = F)

venn_list2<- list(chron.age, age.w, eq3)
names(venn_list2)<- c("chron.age", "age.w", "eq3")

ggvenn(venn_list2,
       text_size = 8,
       show_percentage = F)

venn_all<- list(cross.age, chron.age, age.w, eq3)
names(venn_all)<- c("cross.age", "chron.age", "age.w", "eq3")

upset(fromList(venn_list2), order.by = "freq", text.scale = c(2, 2, 2, 1, 2, 1.5), line.size = 1, point.size = 2)

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
    select(count)
  
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


#Enrichment Analyses------------------------------------------------------------
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
pqlseq_anno<- pqlseq_anno %>%
  mutate('chron-eq2w' = abs(beta_chron_age) - abs(beta_eq2_w_age),
         'chron-eq3' = abs(beta_chron_age) - abs(beta_eq3_age))

func_region<- pqlseq_anno %>%
  select(unique_cpg, anno)

func_region_list<- split(x=func_region$unique_cpg, f=func_region$anno)
func_region_list<- func_region_list[1:29]

run_gsea<- function(){
  
  df<- pqlseq_anno %>%
    select(unique_cpg, `chron-eq2w`) %>%
    arrange(desc(`chron-eq2w`))
  
  chron_within2<- chron_within$`chron-eq2w`
  names(chron_within2) = chron_within$unique_cpg
  
  #Enrichment for Hallmark set
  chron_gsea<- fgsea(pathways = func_region_list, 
                     stats = chron_within2,
                     minSize = 1,
                     maxSize = 500000,
                     eps = 0.0,
                     scoreType = "neg")
  
  chron_gsea<- chron_gsea %>%
    arrange(pathway, NES)
  
}

chron_within<- pqlseq_anno %>%
  select(unique_cpg, `chron-eq2w`) %>%
  arrange(desc(`chron-eq2w`))

chron_within2<- chron_within$`chron-eq2w`
names(chron_within2) = chron_within$unique_cpg

#Enrichment for Hallmark set
chron_gsea<- fgsea(pathways = func_region_list, 
                   stats = chron_within2,
                   minSize = 1,
                   maxSize = 500000,
                   eps = 0.0)

chron_gsea<- chron_gsea %>%
  arrange(pathway, NES)

chron_gsea %>%
  ggplot(aes(x=reorder(pathway, NES), y=NES)) +
  #geom_point(aes(alpha=padj<0.05)) +
  #geom_line(aes(group = type)) +
  geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  #geom_errorbar(yin = test_full$log_ci.lo, ymax = test_full$log_ci.hi, width = 0.3, position = position_dodge(0.5)) +
  #scale_colour_manual(values = c("darkgoldenrod2", 'steelblue2', "green4", 'purple')) +
  theme_classic(base_size =20) +
  #theme(legend.position = "none") +
  #ylim(c(-1, 7)) +
  ylab("NES") +
  xlab("Annotation") +
  coord_flip()
  
chron_eq3<- pqlseq_anno %>%
  select(unique_cpg, `chron-eq3`) %>%
  arrange(desc(`chron-eq3`))

chron_eq3_2<- chron_eq3$`chron-eq3`
names(chron_eq3_2) = chron_eq3$unique_cpg

#Enrichment for Hallmark set
eq3_gsea<- fgsea(pathways = func_region_list, 
                   stats = chron_eq3_2,
                   minSize = 1,
                   maxSize = 500000,
                   eps = 0.0)

eq3_gsea<- eq3_gsea %>%
  arrange(pathway, NES)

chron_gsea %>%
  ggplot(aes(x=reorder(pathway, NES), y=NES)) +
  #geom_point(aes(alpha=padj<0.05)) +
  #geom_line(aes(group = type)) +
  geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  #geom_errorbar(yin = test_full$log_ci.lo, ymax = test_full$log_ci.hi, width = 0.3, position = position_dodge(0.5)) +
  #scale_colour_manual(values = c("darkgoldenrod2", 'steelblue2', "green4", 'purple')) +
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
  
  #ft$anno_source<- "Repeat Elements"
  #ft$anno_source[ft$annotation %in% chmm_ordered]<- "Chromatin States"
  
  #ft<- ft %>%
  #arrange(anno_source, annotation)
  
  #Rearrange factors to sort by type then log_or
  #ft$annotation<- factor(ft$annotation, levels = rev(annos_order))
  
  #ft$type<- var_type
  
  return(ft)
}

pqlseq_anno$within_chron<- "Both Insignificant"
pqlseq_anno$within_chron[pqlseq_anno$fdr_eq2_w_age < 0.05 & pqlseq_anno$fdr_chron_age < 0.05]<- "Both Significant"
pqlseq_anno$within_chron[pqlseq_anno$fdr_eq2_w_age < 0.05 & pqlseq_anno$fdr_chron_age > 0.05]<- "Within Age Significant"
pqlseq_anno$within_chron[pqlseq_anno$fdr_eq2_w_age > 0.05 & pqlseq_anno$fdr_chron_age < 0.05]<- "Chron Age Significant"

#DMRs Significant for EQ2 only
within_chron_hypo<- enrichment(pqlseq_anno, 
                               "within_chron", "eq2_w_signif",
                               "Within Age Significant", "Age-Hypomethylated")
within_chron_hypo$type<- "Age Hypomethylated"

within_chron_hypo$source<- "RE's"
within_chron_hypo$source[within_chron_hypo$annotation %in% chmm_ordered]<- "CHMM"

within_chron_hypo$annotation<- factor(within_chron_hypo$annotation, levels = rev(within_chron_hypo$annotation))

within_chron_hyper<- enrichment(pqlseq_anno, 
                                "within_chron", "eq2_w_signif",
                                "Within Age Significant", "Age-Hypermethylated")
within_chron_hyper$type<- "Age Hypermethylated"

within_chron_hyper$source<- "RE's"
within_chron_hyper$source[within_chron_hyper$annotation %in% chmm_ordered]<- "CHMM"

within_chron_hyper$annotation<- factor(within_chron_hyper$annotation, levels = rev(within_chron_hyper$annotation))

within_chron_enrich<- rbind(within_chron_hypo, within_chron_hyper)

within_chron_enrich_chmm<- within_chron_enrich %>%
  filter(source == "CHMM")

within_chron_enrich_chmm %>%
  ggplot(aes(x=annotation, y=estimate, colour = type)) +
  geom_point(aes(alpha=padj<0.05, shape = type), size = 2) +
  #geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  geom_line(aes(group = type)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_colour_manual(values = c("green1", "green4"), name = "") +
  geom_errorbar(ymin = within_chron_enrich_chmm$conf.low, ymax = within_chron_enrich_chmm$conf.high, width = 0.3) +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  #theme(legend.position = "none") +
  #ylim(c(-6, 4)) +
  ylab("Odds Ratio") +
  xlab("Annotation")

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
  
  return(ft)
}

cross_hypo<- enrichment(pqlseq_anno, "cross_signif",  "Age-Hypomethylated")
cross_hyper<- enrichment(pqlseq_anno, "cross_signif",  "Age-Hypermethylated")

cross_enrich<- rbind(cross_hypo, cross_hyper)

cross_enrich$anno_source[cross_enrich$annotation == "Promoter"]<- "Chromatin States"

cross_enrich$class<- "Simple Repeats"
cross_enrich$class[cross_enrich$annotation %in% c("SINE", "LINE", "LTR", "Retroposon")]<- "TE Class I"
cross_enrich$class[cross_enrich$annotation %in% "DNA"]<- "TE Class II"

cross_enrich_re<- cross_enrich %>% filter(anno_source == "Repeat Elements")
cross_enrich_re$annotation<- factor(cross_enrich_re$annotation, levels = rev(re_ordered))

cross_enrich_chmm<- cross_enrich %>% filter(anno_source != "Repeat Elements")

cross_enrich_chmm$class[cross_enrich_chmm$annotation %in% chmm_ordered[1:2]]<- "TSSs"
cross_enrich_chmm$class[cross_enrich_chmm$annotation %in% chmm_ordered[3:5]]<- "Active Tr."
cross_enrich_chmm$class[cross_enrich_chmm$annotation %in% chmm_ordered[6:8]]<- "Enhancers"
cross_enrich_chmm$class[cross_enrich_chmm$annotation %in% chmm_ordered[9:15]]<- "Quiescent"

chmm_ordered<- as.factor(c("Promoter", as.character(chmm_ordered)))

cross_enrich_chmm$annotation<- factor(cross_enrich_chmm$annotation, levels = rev(chmm_ordered))

cross_enrich_chmm %>%
  ggplot(aes(x=annotation, y=estimate, fill = type)) +
  #geom_point(aes(alpha=padj<0.05)) +
  geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  #geom_line(aes(group = type)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("darkgoldenrod1", "darkgoldenrod4"), name = "") +
  #geom_errorbar(ymin = within_enrich$log_ci.lo, ymax = within_enrich$log_ci.hi, width = 0.3, position = position_dodge(0.5)) +
  theme_classic(base_size = 24) +
  #theme(legend.position = "none") +
  #ylim(c(-6, 4)) +
  ylab("Odds Ratio") +
  xlab("Annotation") +
  coord_flip()

cross_enrich_re %>%
  ggplot(aes(x=annotation, y=estimate, fill = type)) +
  #geom_point(aes(alpha=padj<0.05)) +
  geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  #geom_line(aes(group = type)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("darkgoldenrod1", "darkgoldenrod4"), name = "") +
  #geom_errorbar(ymin = within_enrich$log_ci.lo, ymax = within_enrich$log_ci.hi, width = 0.3, position = position_dodge(0.5)) +
  theme_classic(base_size = 24) +
  #theme(legend.position = "none") +
  #ylim(c(-6, 4)) +
  ylab("Odds Ratio") +
  xlab("Annotation") +
  coord_flip()

#Chronological Age
chron_hypo<- enrichment(pqlseq_anno, "chron_signif", "Age-Hypomethylated")
chron_hyper<- enrichment(pqlseq_anno, "chron_signif", "Age-Hypermethylated")

chron_enrich<- rbind(chron_hypo, chron_hyper)

chron_enrich$anno_source[chron_enrich$annotation == "Promoter"]<- "Chromatin States"

chron_enrich$class<- "Simple Repeats"
chron_enrich$class[chron_enrich$annotation %in% c("SINE", "LINE", "LTR", "Retroposon")]<- "TE Class I"
chron_enrich$class[chron_enrich$annotation %in% "DNA"]<- "TE Class II"

chron_enrich_re<- chron_enrich %>% filter(anno_source == "Repeat Elements")
chron_enrich_re$annotation<- factor(chron_enrich_re$annotation, levels = rev(re_ordered))

chron_enrich_chmm<- chron_enrich %>% filter(anno_source != "Repeat Elements")

chron_enrich_chmm$class[chron_enrich_chmm$annotation %in% chmm_ordered[1:2]]<- "TSSs"
chron_enrich_chmm$class[chron_enrich_chmm$annotation %in% chmm_ordered[3:5]]<- "Active Tr."
chron_enrich_chmm$class[chron_enrich_chmm$annotation %in% chmm_ordered[6:8]]<- "Enhancers"
chron_enrich_chmm$class[chron_enrich_chmm$annotation %in% chmm_ordered[9:15]]<- "Quiescent"

chron_enrich_chmm$annotation<- factor(chron_enrich_chmm$annotation, levels = rev(chmm_ordered))

chron_enrich_chmm %>%
  ggplot(aes(x=annotation, y=estimate, fill = type)) +
  #geom_point(aes(alpha=padj<0.05)) +
  geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  #geom_line(aes(group = type)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("steelblue1", "steelblue4"), name = "") +
  #geom_errorbar(ymin = within_enrich$log_ci.lo, ymax = within_enrich$log_ci.hi, width = 0.3, position = position_dodge(0.5)) +
  theme_classic(base_size = 24) +
  #theme(legend.position = "none") +
  #ylim(c(-6, 4)) +
  ylab("Odds Ratio") +
  xlab("Annotation") +
  coord_flip()

chron_enrich_re %>%
  ggplot(aes(x=annotation, y=estimate, fill = type)) +
  #geom_point(aes(alpha=padj<0.05)) +
  geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  #geom_line(aes(group = type)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("steelblue1", "steelblue4"), name = "") +
  #geom_errorbar(ymin = within_enrich$log_ci.lo, ymax = within_enrich$log_ci.hi, width = 0.3, position = position_dodge(0.5)) +
  theme_classic(base_size = 24) +
  #theme(legend.position = "none") +
  #ylim(c(-6, 4)) +
  ylab("Odds Ratio") +
  xlab("Annotation") +
  coord_flip()

#Eq3 Fisher Enrichment
eq3_hypo<- enrichment(pqlseq_anno, "eq3_age_signif", "Age-Hypomethylated")
eq3_hyper<- enrichment(pqlseq_anno, "eq3_age_signif",  "Age-Hypermethylated")

eq3_enrich<- rbind(eq3_hypo, eq3_hyper)

eq3_enrich$anno_source[eq3_enrich$annotation == "Promoter"]<- "Chromatin States"

eq3_enrich$class<- "Simple Repeats"
eq3_enrich$class[eq3_enrich$annotation %in% c("SINE", "LINE", "LTR", "Retroposon")]<- "TE Class I"
eq3_enrich$class[eq3_enrich$annotation %in% "DNA"]<- "TE Class II"

eq3_enrich_re<- eq3_enrich %>% filter(anno_source == "Repeat Elements")
eq3_enrich_re$annotation<- factor(eq3_enrich_re$annotation, levels = rev(re_ordered))

eq3_enrich_chmm<- eq3_enrich %>% filter(anno_source != "Repeat Elements")

eq3_enrich_chmm$class[eq3_enrich_chmm$annotation %in% chmm_ordered[1:2]]<- "TSSs"
eq3_enrich_chmm$class[eq3_enrich_chmm$annotation %in% chmm_ordered[3:5]]<- "Active Tr."
eq3_enrich_chmm$class[eq3_enrich_chmm$annotation %in% chmm_ordered[6:8]]<- "Enhancers"
eq3_enrich_chmm$class[eq3_enrich_chmm$annotation %in% chmm_ordered[9:15]]<- "Quiescent"

eq3_enrich_chmm$annotation<- factor(eq3_enrich_chmm$annotation, levels = rev(chmm_ordered))

#Within age fisher enrichment
within_hypo<- enrichment(pqlseq_anno, "eq2_w_signif", "Age-Hypomethylated")
within_hyper<- enrichment(pqlseq_anno, "eq2_w_signif", "Age-Hypermethylated")

within_enrich<- rbind(within_hypo, within_hyper)

within_enrich$anno_source[within_enrich$annotation == "Promoter"]<- "Chromatin States"

within_enrich$class<- "Simple Repeats"
within_enrich$class[within_enrich$annotation %in% c("SINE", "LINE", "LTR", "Retroposon")]<- "TE Class I"
within_enrich$class[within_enrich$annotation %in% "DNA"]<- "TE Class II"

within_enrich_re<- within_enrich %>% filter(anno_source == "Repeat Elements")
within_enrich_re$annotation<- factor(within_enrich_re$annotation, levels = rev(re_ordered))

within_enrich_chmm<- within_enrich %>% filter(anno_source != "Repeat Elements")

within_enrich_chmm$class[within_enrich_chmm$annotation %in% chmm_ordered[1:2]]<- "TSSs"
within_enrich_chmm$class[within_enrich_chmm$annotation %in% chmm_ordered[3:5]]<- "Active Tr."
within_enrich_chmm$class[within_enrich_chmm$annotation %in% chmm_ordered[6:8]]<- "Enhancers"
within_enrich_chmm$class[within_enrich_chmm$annotation %in% chmm_ordered[9:15]]<- "Quiescent"

within_enrich_chmm$annotation<- factor(within_enrich_chmm$annotation, levels = rev(chmm_ordered))

#Between age fisher enrichment
between_hypo<- enrichment(pqlseq_anno, "eq2_m_signif", "Age-Hypomethylated")
between_hyper<- enrichment(pqlseq_anno, "eq2_m_signif", "Age-Hypermethylated")

between_enrich<- rbind(between_hypo, between_hyper)

between_enrich$anno_source[between_enrich$annotation == "Promoter"]<- "Chromatin States"

between_enrich$class<- "Simple Repeats"
between_enrich$class[between_enrich$annotation %in% c("SINE", "LINE", "LTR", "Retroposon")]<- "TE Class I"
between_enrich$class[between_enrich$annotation %in% "DNA"]<- "TE Class II"

between_enrich_re<- between_enrich %>% filter(anno_source == "Repeat Elements")
between_enrich_re$annotation<- factor(between_enrich_re$annotation, levels = rev(re_ordered))

between_enrich_chmm<- between_enrich %>% filter(anno_source != "Repeat Elements")

between_enrich_chmm$class[between_enrich_chmm$annotation %in% chmm_ordered[1:2]]<- "TSSs"
between_enrich_chmm$class[between_enrich_chmm$annotation %in% chmm_ordered[3:5]]<- "Active Tr."
between_enrich_chmm$class[between_enrich_chmm$annotation %in% chmm_ordered[6:8]]<- "Enhancers"
between_enrich_chmm$class[between_enrich_chmm$annotation %in% chmm_ordered[9:15]]<- "Quiescent"

between_enrich_chmm$annotation<- factor(between_enrich_chmm$annotation, levels = rev(chmm_ordered))

within_enrich$model<- "within"
chron_enrich$model<- "chron"
cross_enrich$model<- "cross"
eq3_enrich$model<- "eq3"
between_enrich$model<- "between"

all_enrich<- rbind(within_enrich, eq3_enrich, chron_enrich)

chmm_enrich<- all_enrich %>%
  filter(anno_source == "Chromatin States")

chmm_enrich$annotation<- factor(chmm_enrich$annotation, levels = chmm_ordered)

chmm_enrich %>%
  ggplot(aes(x=annotation, y=estimate, colour = model, alpha=padj<0.05)) +
  geom_point() +
  #geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  geom_line(aes(group = model)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_colour_manual(values = c("steelblue1", "purple","green4"), name = "") +
  geom_errorbar(ymin = chmm_enrich$conf.low, ymax = chmm_enrich$conf.high, width = 0.3) +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  #theme(legend.position = "top") +
  #ylim(c(-6, 4)) +
  ylab("Odds Ratio") +
  xlab("Annotation") +
  facet_wrap(vars(type), ncol = 1, scales = "free_y")

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




saveRDS(pqlseq_anno, "/scratch/ckelsey4/Cayo_meth/pqlseq_anno.rds")

#Save workspace image
save.image("/scratch/ckelsey4/Cayo_meth/cross_within_compare.RData")
