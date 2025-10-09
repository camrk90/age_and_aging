library(tidyverse)
library(ggplot2)
library(ggcorrplot)
library(ggvenn)
library(variancePartition)
library(lme4)
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
###          Import Data           ###
######################################
#Metadata
blood_metadata<- read.table("/scratch/ckelsey4/Cayo_meth/blood_metadata_full.txt")
long_data<- read.table("/scratch/ckelsey4/Cayo_meth/long_data_adjusted.txt")

blood_metadata<- blood_metadata[!blood_metadata$lid_pid %in% long_data$lid_pid,]

blood_metadata<- blood_metadata %>%
  drop_na(age_at_sampling) %>%
  filter(age_at_sampling >= 3)

long_data<- long_data %>%
  arrange(lid_pid) %>%
  filter(age_at_sampling > 1) %>%
  dplyr::rename(perc_unique = unique) %>%
  drop_na()

#Import longitudinal pqlseq files
setwd('/scratch/ckelsey4/Cayo_meth/glmer_model_compare')

#Chronological Age
chron_age_files<- 'wb_pqlseq2_agechron'
chron_age_pqlseq<- import_pqlseq(chron_age_files, y = 4)

#Age Within
age_w_files<- 'wb_pqlseq2_within_age'
age_w_pqlseq<- import_pqlseq(age_w_files, y = 5)

#Mean Age
age_m_files<- 'wb_pqlseq2_mean_age'
age_m_pqlseq<- import_pqlseq(age_m_files, y = 5)

#Eq. 3
age_eq3_files<- 'wb_pqlseq2_eq3'
age_eq3_pqlseq<- import_pqlseq(age_eq3_files, y = 4)

#Rename cols for each df to indicate variable
colnames(age_w_pqlseq)<- c("outcome", "length", "chr", "chromStart", "chromEnd", "n", 
                           paste(names(age_w_pqlseq[,7:12]), "within", "age", sep = "_"))
colnames(age_m_pqlseq)<- c(paste(names(age_m_pqlseq), "mean_age", sep = "_"))
colnames(chron_age_pqlseq)<- c(paste(names(chron_age_pqlseq), "chron_age", sep = "_"))

#Cbind cols for age and sex dfs
age_w_pqlseq<- age_w_pqlseq[age_w_pqlseq$outcome %in% chron_age_pqlseq$outcome_chron_age,]
age_m_pqlseq<- age_m_pqlseq[age_m_pqlseq$outcome %in% chron_age_pqlseq$outcome_chron_age,]
all.equal(age_w_pqlseq$outcome, chron_age_pqlseq$outcome_chron_age)
pqlseq_model<- cbind(age_w_pqlseq, age_m_pqlseq[,7:12], chron_age_pqlseq[,7:12])

#Cross Sectional Models
setwd('/scratch/ckelsey4/Cayo_meth/cross_models')

long_cross_files<- '_long'
long_cross_pqlseq<- import_pqlseq(long_cross_files, y = 3)

short_cross_files<- '_short'
short_cross_pqlseq<- import_pqlseq(short_cross_files, y = 3)

pqlseq_model<- pqlseq_model[pqlseq_model$outcome %in% long_cross_pqlseq$outcome,]
long_cross_pqlseq<- long_cross_pqlseq[long_cross_pqlseq$outcome %in% pqlseq_model$outcome,]
short_cross_pqlseq<- short_cross_pqlseq[short_cross_pqlseq$outcome %in% long_cross_pqlseq$outcome,]

all.equal(short_cross_pqlseq$outcome, long_cross_pqlseq$outcome)

long_cross_pqlseq<- long_cross_pqlseq %>%
  select(-c(outcome, length, chr, chromStart, chromEnd, n)) %>%
  mutate_if(is.character, as.numeric)
colnames(long_cross_pqlseq)<- paste(colnames(long_cross_pqlseq), "long", "cross", sep = "_")

short_cross_pqlseq<- short_cross_pqlseq %>%
  select(-c(outcome, length, chr, chromStart, chromEnd, n))
colnames(short_cross_pqlseq)<- paste(colnames(short_cross_pqlseq), "short", "cross", sep = "_")

age_full<- cbind(pqlseq_model, long_cross_pqlseq, short_cross_pqlseq)

#Sort chromosome factors
sorted_labels<- str_sort(unique(age_full$chr), numeric=T)

age_full<- age_full %>% 
  mutate(chr = factor(chr, levels = sorted_labels)) %>%
  arrange(chr) %>%
  mutate(region_range = paste(chromStart, "-", chromEnd, sep = " ")) %>%
  relocate(region_range, .after = chromEnd)

age_full$within_cross<- "Both Insignificant"
age_full$within_cross[age_full$fdr_within_age < 0.05 & age_full$fdr_long_cross < 0.05]<- "Both Significant"
age_full$within_cross[age_full$fdr_within_age < 0.05 & age_full$fdr_long_cross > 0.05]<- "Within Age Significant"
age_full$within_cross[age_full$fdr_within_age > 0.05 & age_full$fdr_long_cross < 0.05]<- "Cross Age Significant"

age_full$within_chron<- "Both Insignificant"
age_full$within_chron[age_full$fdr_within_age < 0.05 & age_full$fdr_chron_age < 0.05]<- "Both Significant"
age_full$within_chron[age_full$fdr_within_age < 0.05 & age_full$fdr_chron_age > 0.05]<- "Within Age Significant"
age_full$within_chron[age_full$fdr_within_age > 0.05 & age_full$fdr_chron_age < 0.05]<- "Chron Age Significant"

#Generate age df subset
age_trunc<- age_full %>%
  select(c(outcome, region_range, chr, beta_within_age, fdr_within_age, beta_mean_age, 
           fdr_mean_age, beta_chron_age, fdr_chron_age, beta_long_cross, fdr_long_cross,
           within_cross, within_chron))

rm(age_w_pqlseq);rm(age_m_pqlseq);rm(chron_age_pqlseq);rm(long_cross_pqlseq)

######################################
###       Descriptive Stats        ###
######################################
#Cross-age age distribution
blood_metadata %>%
  ggplot(aes(age_at_sampling, fill=individual_sex)) +
  geom_histogram(position = "dodge", colour = "black") +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  scale_fill_manual(values = c("steelblue1", "steelblue4"), name = "Sex") +
  theme_classic(base_size=24) +
  ylab("Count") +
  xlab("Age")

#Longitudinal data distribution
long_data %>%
  ggplot(aes(age_at_sampling, fill=individual_sex)) +
  geom_histogram(position = "dodge", colour = "black") +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  scale_fill_manual(values = c("purple", "purple4"), name = "Sex") +
  theme_classic(base_size=24) +
  ylab("Count") +
  xlab("Age")

long_data<- long_data %>%
  group_by(monkey_id) %>%
  mutate(min_age = min(age_at_sampling))
long_data$age_at_sampling<- round(long_data$age_at_sampling, 0)

long_data %>%
  ggplot(aes(x=age_at_sampling, y=reorder(monkey_id, min_age), colour=individual_sex)) +
  geom_path(linewidth = 1.2, alpha = 0.8) +
  geom_point(colour="black") +
  scale_x_continuous(breaks = seq(0, 30, by=5)) +
  scale_colour_manual(values = c("purple", "purple4"), name = "Sex") +
  ylab("Individual") +
  xlab("Age") +
  theme_classic(base_size = 24) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

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

######################################
###       Plot Distributions       ###
######################################
compare_plot<- function(df, fdr1, fdr2, var1, var2, c1, c2, lab1, lab2, plot_type) {
  
  df<- df %>%
    filter({{fdr1}} < .05 & {{fdr2}} < .05) %>%
    mutate(diff = abs({{var2}}) - abs({{var1}}))
  
  #eval(substitute(df_lm<- lm(var1 ~ var2, data=df)))
  
  #print(summary(df_lm))
  
  if (plot_type == "scatter"){
    
    df %>%
      ggplot(aes({{var1}}, {{var2}}, colour = diff)) +
      geom_point() +
      geom_abline() +
      #geom_abline(slope = df_lm[["coefficients"]][[2]], 
      #intercept = df_lm[["coefficients"]][[1]],
      #colour = "red") +
      geom_smooth(method = "lm") +
      geom_vline(xintercept=0, linetype="dashed") +
      geom_hline(yintercept=0, linetype="dashed") +
      scale_color_gradient2(low = c2, mid = "grey70", high = c1, midpoint = 0, name = "") +
      theme_classic(base_size=32) +
      theme(legend.key.height= unit(2, 'cm')) +
      theme(panel.background = element_rect(colour = "black", linewidth=3)) +
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
      theme_classic(base_size=32) +
      theme(legend.key.height= unit(2, 'cm')) +
      theme(panel.background = element_rect(colour = "black", linewidth=3)) +
      xlab(paste(lab2, "-", lab1, sep=" "))
    
  }
}

#Full CS set plots
compare_plot(age_full, fdr_long_cross, fdr_chron_age,
             beta_long_cross, beta_chron_age,
             "darkgoldenrod2", "steelblue2",
             "Cross Age", "Chronological Age", "scatter")

compare_plot(age_full, fdr_long_cross, fdr_chron_age,
             beta_long_cross, beta_chron_age,
             "darkgoldenrod2", "steelblue2",
             "Cross Age", "Chronological Age", "hist")

compare_plot(age_full, fdr_long_cross, fdr_mean_age,
             beta_long_cross, beta_mean_age,
             "chartreuse3", "steelblue2",
             "Full Cross Age", "Mean Age")

compare_plot(age_full, fdr_long_cross, fdr_within_age,
             beta_long_cross, beta_within_age,
             "purple", "steelblue2",
             "Cross Age", "Within Age", "scatter")

compare_plot(age_full, fdr_long_cross, fdr_within_age,
             beta_long_cross, beta_within_age,
             "purple", "steelblue2",
             "Cross Age", "Within Age", "hist")

compare_plot(age_full, fdr_chron_age, fdr_within_age,
             beta_chron_age, beta_within_age,
             "purple", "darkgoldenrod2",
             "Chronological Age", "Within Age", "scatter")

compare_plot(age_full, fdr_chron_age, fdr_within_age,
             beta_chron_age, beta_within_age,
             "purple", "darkgoldenrod2",
             "Chron Age", "Within Age", "hist")

age.w.count<- nrow(age_full[age_full$fdr_within_age < 0.05,])
age.mean.count<- nrow(age_full[age_full$fdr_mean_age < 0.05,])
age.chron.count<- nrow(age_full[age_full$fdr_chron_age < 0.05,])
age.long.count<- nrow(age_full[age_full$fdr_long_cross < 0.05,])
counts<- data.frame(count = c(age.w.count, age.mean.count, age.chron.count, age.long.count),
                    predictor = as.factor(c('Within Age', 'Mean Age', 'Chron Age', 'Cross Age')))

counts<- counts %>%
  mutate(predictor = as.factor(predictor)) %>%
  mutate(predictor=fct_reorder(predictor, count, .desc=T))
counts$predictor

rm(age.w.count);rm(age.mean.count);rm(age.chron.count);rm(age.long.count)

counts %>%
  ggplot(aes(predictor, count, fill = predictor)) +
  geom_bar(stat = 'identity', colour="black") +
  geom_text(label=counts$count, vjust=-1, size=5) +
  theme_classic(base_size = 32) +
  theme(axis.text.x = element_text(angle = 15, hjust=0.9)) +
  xlab("Predictor") +
  ylab("Count") +
  scale_fill_manual(values = c("purple", 'steelblue2', 'darkgoldenrod2', 'chartreuse3'))

#Distribution of effect sizes for each variable
age_full %>%
  dplyr::select(c(beta_within_age, beta_mean_age, beta_chron_age, beta_long_cross)) %>%
  pivot_longer(cols = c(beta_within_age, beta_mean_age, beta_chron_age, beta_long_cross),
               values_to = 'beta',
               names_to = 'var') %>%
  ggplot(aes(beta, fill=var)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_fill_manual(values = c('darkgoldenrod2', 'steelblue2', 'chartreuse3', "purple"),
                    labels = c("Chron Age", "Cross Age", "Mean Age", "Within Age")) +
  theme_classic(base_size = 24) +
  xlim(-0.25, 0.25)

######################################
###           CONCORDANCE          ###   
######################################

age_trunc %>% 
  filter(within_chron != "Both Insignificant") %>%
  ggplot(aes(beta_chron_age, beta_within_age, colour = within_chron)) +
  geom_point(alpha = 0.5) +
  geom_abline() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_classic(base_size=24) +
  theme(panel.background = element_rect(colour = "black", linewidth=3)) +
  xlim(-0.25, 0.25) +
  ylim(-0.25, 0.25) +
  xlab("Chronological Age") +
  ylab("Within Age")

age_trunc %>% 
  filter(within_cross != "Both Insignificant") %>%
  ggplot(aes(beta_long_cross, beta_within_age, colour = within_cross)) +
  geom_point(alpha = 0.5) +
  geom_abline() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_classic(base_size=24) +
  theme(panel.background = element_rect(colour = "black", linewidth=3)) +
  #theme(legend.position = "none") +
  xlim(-0.25, 0.25) +
  ylim(-0.25, 0.25) +
  xlab("Cross Age") +
  ylab("Within Age")

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

#CHMM---------------------------------------------------------------------------
#Join pqlseq model and chmm 
pqlseq_chmm<- left_join(chmm_intersect, age_trunc, by = c("region_range", "chr"))

#Filter out regions with models that didn't converge resulting in NAs in the annotation join
pqlseq_chmm<- pqlseq_chmm %>%
  drop_na()

#Set annotations as factor and reorder
annotations_ordered<- str_sort(unique(pqlseq_chmm$anno), numeric = TRUE)
pqlseq_chmm$anno<- factor(pqlseq_chmm$anno, levels = rev(annotations_ordered))

#Create column of broad categories for annotations
#pqlseq_chmm$anno_class<- "A"
pqlseq_chmm$anno_class<- "TSSs"
pqlseq_chmm$anno_class[pqlseq_chmm$anno %in% annotations_ordered[3:5]]<- "Active Tr."
pqlseq_chmm$anno_class[pqlseq_chmm$anno %in% annotations_ordered[6:8]]<- "Enhancers"
pqlseq_chmm$anno_class[pqlseq_chmm$anno %in% annotations_ordered[9:15]]<- "Quiescent"

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
pqlseq_re$anno_class[pqlseq_re$repClass %in% c("rRNA", "snRNA", "tRNA", "srpRNA", "scRNA")]<- "Str. RNAs"

pqlseq_re<- pqlseq_re %>%
  dplyr::rename(anno = repClass) %>%
  dplyr::select(-repName) %>%
  dplyr::relocate(anno, .after = anno_end) %>%
  dplyr::select(-c(range))

#Bind annotation dfs together---------------------------------------------------
pqlseq_anno<- rbind(pqlseq_chmm, pqlseq_re, pqlseq_prom)

annos<- c("Promoter", "TSSs", "Active Tr.", "Enhancers", "Quiescent", "Simple Repeats",
          "TE Class I", "TE Class II", "Str. RNAs")

pqlseq_anno$anno_source<- "Repeat Elements"
pqlseq_anno$anno_source[pqlseq_anno$anno_class == "Transcription Start Sites" | pqlseq_anno$anno_class == "Active Transcription" |
                          pqlseq_anno$anno_class == "Enhancer Regions" | pqlseq_anno$anno_class == "Quiescent States" | 
                          pqlseq_anno$anno_class == "Promoter"]<- "Transcription"

pqlseq_anno<- pqlseq_anno %>%
  arrange(anno_source, anno_class)

#Rearrange factors to sort by type then log_or
pqlseq_anno$anno_class<- factor(pqlseq_anno$anno_class, levels = rev(annos))

pqlseq_anno$unique_cpg<- paste(pqlseq_anno$chr, pqlseq_anno$cpg_loc, sep="_")

pqlseq_anno<- pqlseq_anno %>%
  dplyr::relocate(c(anno_class, anno_source, unique_cpg), .after=anno)

#Annotation proportions
pqlseq_anno$cross_signif<- "Non-Significant"
pqlseq_anno$cross_signif[pqlseq_anno$fdr_long_cross < 0.05 & pqlseq_anno$beta_long_cross < 0]<- "Age-Hypomethylated"
pqlseq_anno$cross_signif[pqlseq_anno$fdr_long_cross < 0.05 & pqlseq_anno$beta_long_cross > 0]<- "Age-Hypermethylated"

pqlseq_anno$cross_signif<- factor(pqlseq_anno$cross_signif, 
                        levels = c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))

pqlseq_anno$chron_signif<- "Non-Significant"
pqlseq_anno$chron_signif[pqlseq_anno$fdr_chron_age < 0.05 & pqlseq_anno$beta_chron_age < 0]<- "Age-Hypomethylated"
pqlseq_anno$chron_signif[pqlseq_anno$fdr_chron_age < 0.05 & pqlseq_anno$beta_chron_age > 0]<- "Age-Hypermethylated"

pqlseq_anno$chron_signif<- factor(pqlseq_anno$chron_signif, 
                                  levels = c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))

pqlseq_anno$within_signif<- "Non-Significant"
pqlseq_anno$within_signif[pqlseq_anno$fdr_within_age < 0.05 & pqlseq_anno$beta_within_age < 0]<- "Age-Hypomethylated"
pqlseq_anno$within_signif[pqlseq_anno$fdr_within_age < 0.05 & pqlseq_anno$beta_within_age > 0]<- "Age-Hypermethylated"

pqlseq_anno$within_signif<- factor(pqlseq_anno$within_signif, 
                                  levels = c("Age-Hypermethylated", "Non-Significant", "Age-Hypomethylated"))

pqlseq_anno$direction<- "Non-Significant"
pqlseq_anno$direction[(pqlseq_anno$beta_within_age < 0 & pqlseq_anno$beta_long_cross < 0) & 
                        (pqlseq_anno$fdr_within_age < .05 & pqlseq_anno$fdr_long_cross < .05)]<-"Both Negative"
pqlseq_anno$direction[(pqlseq_anno$beta_within_age > 0 & pqlseq_anno$beta_long_cross > 0) & 
                        (pqlseq_anno$fdr_within_age < .05 & pqlseq_anno$fdr_long_cross < .05)]<- "Both Positive"
pqlseq_anno$direction[(pqlseq_anno$beta_within_age < 0 & pqlseq_anno$beta_long_cross > 0) & 
                        (pqlseq_anno$fdr_within_age < .05 & pqlseq_anno$fdr_long_cross < .05)]<- "Within Negative, Cross Positive"
pqlseq_anno$direction[(pqlseq_anno$beta_within_age > 0 & pqlseq_anno$beta_long_cross < 0) & 
                        (pqlseq_anno$fdr_within_age < .05 & pqlseq_anno$fdr_long_cross < .05)]<- "Within Positive, Cross Negative"

#Plot annotation proportions----------------------------------------------------
generate_proportion<- function(df, x){
  
  d1<- df %>% 
    distinct(unique_cpg, .keep_all = T) %>%
    group_by(anno_class, {{x}}) %>% 
    summarise(count = n()) %>% 
    mutate(perc = count/sum(count))
  
  print(d1)
  
  d2<- df %>%
    distinct(unique_cpg, .keep_all = T) %>%
    group_by({{x}}) %>%
    summarise(count = n()) %>%
    mutate(perc = count/sum(count))
  
  d2$anno_class<- "All"
  
  d3<- rbind(d1, d2)
  annos2<- unique(d3$anno_class)
  d3$anno_class<- factor(d3$anno_class, levels = annos2)
  
  col1 <- eval(substitute(x), d2)
  
  d3$percent<- d3$perc*100
  
  d3 %>%
    ggplot(aes(x = percent, y=anno_class, fill = factor({{x}}))) +
    geom_bar(stat="identity", width = 0.7, colour="black") +
    theme_classic(base_size=32) +
    #geom_vline(xintercept = (1-d2$perc[col1 == "Age-Hypermethylated"])*100, linetype = 'dashed') +
    #geom_vline(xintercept = d2$perc[col1 == "Age-Hypomethylated"]*100, linetype = 'dashed') +
    #theme(legend.position = "none") +
    #scale_fill_manual(values = c(c1, c2, c3)) +
    ylab("Annotation") +
    xlab("Percentage")
  
}

generate_proportion(pqlseq_anno, cross_signif, "steelblue2", "gray90", "steelblue4")

generate_proportion(pqlseq_anno, within_signif, "purple4", "gray90", "purple1")

generate_proportion(pqlseq_anno, chron_signif, "darkgoldenrod4", "gray90", "darkgoldenrod1")

generate_proportion(pqlseq_anno, within_chron)


#Enrichment Analyses------------------------------------------------------------

enrichment<- function(model_df, cross_type, var_type){
  
  df_list<- list()
  
  for(i in unique(model_df$anno_class)) {
    
    df2<- model_df %>%
      distinct(unique_cpg, .keep_all = T)
    
    df3<- model_df %>%
      distinct(unique_cpg, .keep_all = T) %>%
      filter(!unique_cpg %in% df2$unique_cpg)
      
      #Counts for fdr < 0.05 & cpg is Within Positive, Cross Negative
      a<- nrow(df2[df2[[cross_type]] == var_type & df2$anno_class == i,])
      b<- nrow(df2[df2[[cross_type]] == var_type & df2$anno_class != i,])
      
      #Counts for NOT fdr < 0.05 & cpg is Within Positive, Cross Negative
      c<- nrow(df2[df2[[cross_type]] != var_type & df2$anno_class == i,])
      d<- nrow(df2[df2[[cross_type]] != var_type & df2$anno_class != i,])
      
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
  names(df_list)<- unique(model_df$anno_class)
  
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
  
  ft$anno_source<- "Repeat Elements"
  ft$anno_source[ft$annotation == "TSSs" | ft$annotation == "Active Tr." |
                   ft$annotation == "Enhancers" | ft$annotation == "Quiescent" | 
                   ft$annotation == "Promoter"]<- "Transcription"
  ft<- ft %>%
    arrange(anno_source, log_or)
  hyper_levels<- as.character(ft$annotation)
  
  #Rearrange factors to sort by type then log_or
  ft$annotation<- factor(ft$annotation, levels = rev(annos))
  
  ft$type<- var_type
  
  return(ft)
}

unique(pqlseq_anno$within_cross)

both_enrich<- enrichment(pqlseq_anno, "Both Significant")
within_enrich<- enrichment(pqlseq_anno, "Within Age Significant")
cross_enrich<- enrichment(pqlseq_anno, "Cross Age Significant")

both_enrich<- enrichment(pqlseq_anno,  "Both Significant")
within_enrich<- enrichment(pqlseq_anno, "within_chron", "Within Age Significant")
chron_enrich<- enrichment(pqlseq_anno, "within_chron", "Chron Age Significant")

test_full<- rbind(both_enrich, within_enrich, chron_enrich)

test_full %>%
  ggplot(aes(x=annotation, y=estimate, fill=type, alpha=padj<0.05)) +
  geom_col(position = position_dodge(0.5), colour="black") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  #geom_errorbar(ymin = test$conf.low, ymax = test$conf.high, width = 0.3, position = position_dodge(0.7)) +
  #scale_fill_manual(values = c("hotpink", "hotpink3")) +
  theme_classic(base_size = 32) +
  #theme(legend.position = "none") +
  #ylim(c(-1, 7)) +
  ylab("Odds Ratio") +
  xlab("Annotation") +
  coord_flip()

#Save workspace image
save.image("/scratch/ckelsey4/Cayo_meth/cross_within_compare.RData")

######################################
###          Hip Flexion           ###
######################################

hip_flexion<- blood_metadata<- read_csv("/scratch/ckelsey4/Cayo_meth/hip_flexion.csv", col_names = T)
hip_flexion<- hip_flexion %>%
  group_by(individual_code) %>%
  mutate(n = n()) %>%
  filter(n > 1) %>%
  mutate(between_age = mean(age),
         within_age = age - between_age) %>%
  relocate(within_age, .after = age) %>%
  relocate(between_age, .after = within_age)

hip_extension_within<- lmer(hip_extension_deg ~ within_age + between_age + individual_sex + (1|individual_code), 
                            data = hip_flexion)

hip_extension_chron<- lmer(hip_extension_deg ~ age + individual_sex + (1|individual_code), 
                           data = hip_flexion)

summary(hip_extension_within)[["coefficients"]]
summary(hip_extension_chron)[["coefficients"]]

hip_rotation_within<- lmer(hip_external_rotation_deg ~ within_age + between_age + individual_sex + (1|individual_code), 
                           data = hip_flexion)

hip_rotation_chron<- lmer(hip_external_rotation_deg ~ age + individual_sex + (1|individual_code), 
                          data = hip_flexion)

summary(hip_rotation_within)[["coefficients"]]
summary(hip_rotation_chron)[["coefficients"]]

#Eq2 vs Eq3--------------
age_eq3<- age_eq3_pqlseq %>%
  select(c(outcome, beta, fdr))
colnames(age_eq3)<- c("outcome", "eq3_beta", "eq3_fdr")
age_w<- age_w_pqlseq  %>%
  select(c(outcome, beta, fdr))
colnames(age_w)<- c("outcome", "age_w_beta", "age_w_fdr")
age_w<- age_w[age_w$outcome %in% age_eq3$outcome,]
age_w<- age_w %>%
  select(-outcome)
age_compare<- cbind(age_w, age_eq3)

age_compare %>%
  #filter(age_w_fdr < .05 & eq3_beta < .05) %>%
  ggplot(aes(age_w_beta, eq3_beta)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic(base_size=24)








