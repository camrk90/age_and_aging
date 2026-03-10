library(tidyverse)
library(ggplot2)
library(ggpubr)

#GLMER--------------------------------------------------------------------------
#Generate Import Function
import_glmer<- function(x, y){
  
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
 
}

setwd('/scratch/ckelsey4/Cayo_meth/glmer_model_compare')

#DNA Methylation----------------------------------------------------------------
#Intercepts
eq2_int_files<- 'glmer_eq2_int_'
eq2_int_glmer<- import_glmer(eq2_int_files, y = 3)

eq2_int<- eq2_int_glmer %>%
  dplyr::select(estimate_age.w, estimate_age.m, pval_age.w, pval_age.m)
colnames(eq2_int)<- paste("eq2_int", colnames(eq2_int), sep = "_")

eq3_int_files<- 'glmer_eq3_int_'
eq3_int_glmer<- import_glmer(eq3_int_files, y = 3)

eq3_int<- eq3_int_glmer %>%
  dplyr::select(estimate_age, estimate_age.m, pval_age, pval_age.m)
colnames(eq3_int)<- paste("eq3_int", colnames(eq3_int), sep = "_")

int_compare<- cbind(eq2_int, eq3_int)

nrow(int_compare[int_compare$eq2_int_pval_age.w < .05,])
nrow(int_compare[int_compare$eq3_int_pval_age < .05,])

int_compare %>%
  ggplot(aes(eq2_int_estimate_age.w, eq3_int_estimate_age)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic(base_size=24) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        #plot.margin = margin(0, 0, 0, 0, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  stat_cor()
  
#Slopes
eq2_slopes_files<- 'glmer_eq2_slopes_'
eq2_slopes_glmer<- import_glmer(eq2_slopes_files, y = 3)

eq2_slopes<- eq2_slopes_glmer %>%
  dplyr::select(estimate_age.w, estimate_age.m, pval_age.w, pval_age.m)
colnames(eq2_slopes)<- paste("eq2_slopes", colnames(eq2_slopes), sep = "_")

eq3_slopes_files<- 'glmer_eq3_slopes_'
eq3_slopes_glmer<- import_glmer(eq3_slopes_files, y = 3)

eq3_slopes<- eq3_slopes_glmer %>%
  dplyr::select(estimate_age, estimate_age.m, pval_age, pval_age.m)
colnames(eq3_slopes)<- paste("eq3_slopes", colnames(eq3_slopes), sep = "_")

slopes_compare<- cbind(eq2_slopes, eq3_slopes)

nrow(slopes_compare[slopes_compare$eq2_slopes_pval_age.w < .05,])
nrow(slopes_compare[slopes_compare$eq3_slopes_pval_age < .05,])

slopes_compare %>%
  ggplot(aes(eq2_slopes_estimate_age.w, eq3_slopes_estimate_age)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic(base_size=24) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        #plot.margin = margin(0, 0, 0, 0, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  stat_cor()


#RNA-seq------------------------------------------------------------------------
setwd('/home/ckelsey4/rna_data')
#Intercepts
rna_eq2_int_files<- 'eq2_int_glmer'
rna_eq2_int_glmer<- import_glmer(rna_eq2_int_files, y = 3)

rna_eq2_int<- rna_eq2_int_glmer %>%
  select(estimate_age.w, estimate_age.m, pval_age.w, pval_age.m)
colnames(rna_eq2_int)<- paste("eq2_int", colnames(rna_eq2_int), sep = "_")

rna_eq3_int_files<- 'eq3_int_glmer'
rna_eq3_int_glmer<- import_glmer(rna_eq3_int_files, y = 3)

rna_eq3_int<- rna_eq3_int_glmer %>%
  select(estimate_age, estimate_age.m, pval_age, pval_age.m)
colnames(rna_eq3_int)<- paste("eq3_int", colnames(rna_eq3_int), sep = "_")

int_compare<- cbind(rna_eq2_int, rna_eq3_int)

int_compare %>%
  ggplot(aes(eq2_int_estimate_age.w, eq3_int_estimate_age)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic(base_size=24) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        #plot.margin = margin(0, 0, 0, 0, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  stat_cor()

#Slopes
eq2_slopes_files<- 'eq2_slope_glmer'
eq2_slopes_glmer<- import_glmer(eq2_slopes_files, y = 3)

eq2_slopes<- eq2_slopes_glmer %>%
  select(estimate_age.w, estimate_age.m, pval_age.w, pval_age.m)
colnames(eq2_slopes)<- paste("eq2_slopes", colnames(eq2_slopes), sep = "_")

eq3_slopes_files<- 'eq3_slope_glmer'
eq3_slopes_glmer<- import_glmer(eq3_slopes_files, y = 3)

eq3_slopes<- eq3_slopes_glmer %>%
  select(estimate_age, estimate_age.m, pval_age, pval_age.m)
colnames(eq3_slopes)<- paste("eq3_slopes", colnames(eq3_slopes), sep = "_")

slopes_compare<- cbind(eq2_slopes, eq3_slopes)

slopes_compare %>%
  drop_na() %>%
  ggplot(aes(eq2_slopes_estimate_age.w, eq3_slopes_estimate_age)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic(base_size=24) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        #plot.margin = margin(0, 0, 0, 0, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  stat_cor()
 




