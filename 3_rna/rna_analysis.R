library(tidyverse)
library(ggplot2)
library(ggeffects)
library(ggcorrplot)
library(ggpubr)
library(ggrepel)
library(fgsea)
library(msigdbr)
library(Biostrings)
library(GenomicFeatures)
library(GenomicRanges)
library(biomaRt)
library(ggvenn)
library(UpSetR)

load("/home/ckelsey4/rna_data/rna_analysis.RData")

#Load data
eq2_int<- readRDS("/home/ckelsey4/rna_data/rna_eq2_int")
eq2_slopes<- readRDS("/home/ckelsey4/rna_data/rna_eq2_slopes")

eq3_int<- readRDS("/home/ckelsey4/rna_data/rna_eq3_int")
eq3_slopes<- readRDS("/home/ckelsey4/rna_data/rna_eq3_slopes")

chron_int<- readRDS("/home/ckelsey4/rna_data/rna_chron_int")
chron_slopes<- readRDS("/home/ckelsey4/rna_data/rna_chron_slopes")

#Make simplified outcome df
rna_int<- as.data.frame(cbind(chron_int$outcome, chron_int$beta_trapped_age, chron_int$pvalue_trapped_age,
                             eq2_int$beta_within_age, eq2_int$pvalue_within_age,
                             eq2_int$beta_mean_age, eq2_int$pvalue_mean_age,
                             eq3_int$beta_trapped_age, eq3_int$pvalue_trapped_age,
                             eq3_int$beta_mean_age, eq3_int$pvalue_mean_age))

c_names<- c("outcome", "beta_chron_age", "pval_chron_age",
                     "beta_eq2_w", "pval_eq2_w",
                     "beta_eq2_m", "pval_eq2_m",
                     "beta_eq3_age", "pval_eq3_age",
                     "beta_eq3_m", "pval_eq3_m")
colnames(rna_int)<- c_names

rna_int<- rna_int %>%
  mutate(across(2:11, as.numeric)) %>%
  mutate(eq2_eq3_diff = abs(beta_eq2_w) - abs(beta_eq3_age))

rna_slopes<- as.data.frame(cbind(chron_slopes$outcome, chron_slopes$beta_trapped_age, chron_slopes$pvalue_trapped_age,
                              eq2_slopes$beta_within_age, eq2_slopes$pvalue_within_age,
                              eq2_slopes$beta_mean_age, eq2_slopes$pvalue_mean_age,
                              eq3_slopes$beta_trapped_age, eq3_slopes$pvalue_trapped_age,
                              eq3_slopes$beta_mean_age, eq3_slopes$pvalue_mean_age))

colnames(rna_slopes)<- c_names

rna_slopes<- rna_slopes %>%
  mutate(across(2:11, as.numeric)) %>%
  mutate(eq2_chron_diff = abs(beta_chron_age) - abs(beta_eq2_w),
         eq2_eq3_diff = abs(beta_eq2_w) - abs(beta_eq3_age),
         eq2w_eq2m_diff = beta_eq2_m - beta_eq2_w)

#Compare slopes vs intercepts---------------------------------------------------
rna_int %>%
  ggplot(aes(beta_eq2_w, beta_eq3_age)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic(base_size=24) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        #plot.margin = margin(0, 0, 0, 0, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  stat_cor() +
  ggtitle("Intercept")

rna_slopes %>%
  ggplot(aes(beta_eq2_w, beta_eq3_age)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic(base_size=24) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  stat_cor() +
  ggtitle("Slopes")

rna_slopes %>%
  ggplot(aes(beta_eq3_m, eq2w_eq2m_diff)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic(base_size=24) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        #plot.margin = margin(0, 0, 0, 0, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  stat_cor() +
  xlab("Beta Eq3 Between Age") +
  ylab("(Beta Eq2 Within) - (Beta Eq2 Between)")

#PCA----------------------------------------------------------------------------
rna_pca<- prcomp(rna_norm, center = TRUE, scale. = TRUE)

pcs<- as.data.frame(rna_pca$x) 

#Check which pca's explain the most variance
summary(rna_pca)$importance[2, ]

pcs<- cbind(pcs[1:2], base_meta)

pcs %>% ggplot(aes(x = pid, y = age_at_sampling)) +
  geom_boxplot()

pcs %>% ggplot(aes(x = prep_year, y = age_at_sampling)) +
  geom_boxplot()

pcs %>%
  ggplot(aes(within_age, PC1)) +
  geom_point() +
  geom_smooth(method = "lm")

pc.matrix<- model.matrix(~ PC1 + PC2 + trapped_age + within_age + mean_age + sex, data = pcs)
pc.matrix %>% 
  cor(use="pairwise.complete.obs") %>%
  ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2)

# Plot Model Outcomes ----------------------------------------------------------
#Effect sizes
compare_plot<- function(df, fdr1, fdr2, var1, var2, plot_type) {
  
  v1<-deparse(substitute(var1))
  v2<-deparse(substitute(var2))
  
  df<- df %>%
    filter({{fdr1}} < .05 & {{fdr2}} < .05) %>%
    mutate(diff = abs({{var2}}) - abs({{var1}}),
           ratio = abs({{var2}})/abs({{var1}}))
  
  print(paste("Median", v2, "/", v1, "=", median(df$ratio), sep = " "))
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
      ggplot(aes({{var1}}, {{var2}}, colour = ratio)) +
      geom_point(size = 1, alpha = 0.8) +
      geom_abline() +
      geom_smooth(method = "lm", linewidth = 0.5) +
      geom_vline(xintercept=0, linetype="dashed") +
      geom_hline(yintercept=0, linetype="dashed") +
      theme_classic(base_size = 5) +
      theme(legend.key.width = unit(1, 'mm'), 
            legend.key.height = unit(5, 'mm'),
            legend.position = "right") +
      theme(panel.background = element_rect(colour = "black", linewidth=1),
            axis.line = element_line(colour = "black", linewidth = 0.5),
            plot.margin = margin(1, 1, 1, 1, "pt"),
            aspect.ratio = 1,
            panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
            panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
      xlim(-1.0, 1.0) +
      ylim(-1.0, 1.0) 
    
  } else if (plot_type == "hist") {
    
    df %>%
      ggplot(aes(ratio, fill = after_stat(x))) +
      geom_histogram(bins = 50) +
      geom_vline(xintercept=1, linetype="dashed") +
      geom_vline(xintercept=median(df$ratio), linetype="dashed", colour = 'red') +
      theme_classic(base_size = 5) +
      theme(legend.key.width = unit(1, 'mm'), 
            legend.key.height = unit(5, 'mm'),
            legend.position = "right") +
      theme(panel.background = element_rect(colour = "black", linewidth=1),
            axis.line = element_line(colour = "black", linewidth = 0.5),
            plot.margin = margin(1, 1, 1, 1, "pt"),
            aspect.ratio = 1,
            panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
            panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) 
    
  }
}

within_chron_plot_rna<- compare_plot(rna_int, pval_chron_age, pval_eq2_w, 
                                     beta_chron_age, beta_eq2_w, "scatter")

within_chron_plot_rna<- within_chron_plot_rna +
  scale_color_gradient2(low = "steelblue2", mid = "grey70", high = "green4", midpoint = 0, 
                        name = "") +
  xlab("Beta Eq.1 Age") +
  ylab("Beta Eq.2 Age Within")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/within_chron_scatter_rna.svg", 
       within_chron_plot_rna, 
       height = 55, width = 55, units = "mm")

within_chron_hist_rna<- compare_plot(rna_int, pval_chron_age, pval_eq2_w, 
                                     beta_chron_age, beta_eq2_w,"hist")
within_chron_hist_rna<- within_chron_hist_rna +
  scale_fill_gradient2(low = "steelblue2", mid = "grey70", high = "green4", midpoint = 0, name = "") +
  xlab('| Eq.2 Age Within |/| Eq.1 |') +
  ylab("Count")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/within_chron_hist_rna.svg", 
       within_chron_hist_rna, 
       height = 55, width = 55, units = "mm")

compare_plot(rna_int, pval_chron_age, pval_eq3_age,
             beta_chron_age, beta_eq3_age,
             "purple", "steelblue2", "scatter")

compare_plot(rna_int, pval_chron_age, pval_eq3_age,
             beta_chron_age, beta_eq3_age,
             "purple", "steelblue2", 
             "Chron Age", "Eq3. Age", "hist")

compare_plot(rna_int, pval_eq2_w, pval_eq3_age,
             beta_eq2_w, beta_eq3_age,
             "purple", "green4", "scatter")

compare_plot(rna_int, pval_eq2_w, pval_eq3_age,
             beta_eq2_w, beta_eq3_age,
             "purple", "green4","hist")

#Counts for significant genes
generate_counts<- function(rna) {
  
  age.chron.count<- nrow(rna[rna$pval_chron_age < 0.05,])
  age.w.count<- nrow(rna[rna$pval_eq2_w < 0.05,])
  age.eq3.count<- nrow(rna[rna$pval_eq3_age < 0.05,])
  total<- data.frame(count = c(age.w.count, age.eq3.count, age.chron.count),
                      predictor = as.factor(c('Eq2. Within Age', 'Eq3 Age', 'Chron Age')))
  
  total<- total %>%
    mutate(predictor = as.factor(predictor)) %>%
    mutate(predictor=fct_reorder(predictor, count, .desc=T),
           perc_signif = count/nrow(rna))
  
  return(total)
  
}

int_counts<- generate_counts(rna_int)
slope_counts<- generate_counts(rna_slopes)

signif_counts<- counts %>%
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
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/signif_counts_rna.svg", 
       signif_counts, 
       height = 20, width = 30, units = "mm")

beta_dist_rna<- rna_int %>%
  dplyr::select(c(beta_eq2_w, beta_chron_age)) %>%
  pivot_longer(cols = c(beta_eq2_w, beta_chron_age),
               values_to = 'beta',
               names_to = 'var') %>%
  ggplot(aes(beta, fill=var)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_fill_manual(values = c('steelblue2', 'green4'),
                    labels = c("Chron Age", "Eq2. Within Age")) +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt")) +
  ylab("Density") +
  xlab("Beta")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/beta_dist_rna.svg", 
       beta_dist_rna, 
       height = 45, width = 60, units = "mm")

## Significant regions Venn diagram
chron.age<- rna_int$outcome[rna_int$pval_chron_age < 0.05]
age.w<- rna_int$outcome[rna_int$pval_eq2_w < 0.05]
eq2.btwn<- rna_int$outcome[rna_int$pval_eq2_m < 0.05]

venn_list<- list(chron.age, age.w, eq2.btwn)
names(venn_list)<- c("chron.age", "age.w", "age.btwn")

ggvenn(venn_list,
       text_size = 8,
       show_percentage = F)

venn_all<- list(cross.age, chron.age, age.w, eq3)
names(venn_all)<- c("cross.age", "chron.age", "age.w", "eq3")

upset(fromList(venn_list), order.by = "freq", 
      text.scale = c(2, 2, 2, 1, 2, 1.5), 
      line.size = 1, point.size = 2)

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/n_samples.svg", 
       n_samples, 
       height = 16, width = 24, units = "mm")

#Plot top genes-----------------------------------------------------------------
#Collect all macaque genes
mm_mart<- useEnsembl(biomart="genes", dataset="mmulatta_gene_ensembl", mirror = 'useast')
mm_genes<- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                 mart = mm_mart)
colnames(mm_genes)<- c("anno", "gene_name")

#Replace ENSMMUG names with gene names 
rna_genes<- rna_int$outcome
rna_genes2<- rna_genes[grepl("ENSMMUG*", rna_genes)]
mm_genes2<- mm_genes[mm_genes$anno %in% rna_genes2, ]

mm_genes2 <- mm_genes2 %>%
  mutate(gene_name = ifelse(gene_name == "", anno, gene_name))

rna_int2<- rna_int

rna_int$outcome[match(mm_genes2$anno, rna_int$outcome)]<- mm_genes2$gene_name

top20<- rna_int %>%
  arrange(desc(beta_chron_age))

top20<- top20[c(1:10, (nrow(top20)-9):nrow(top20)), ]

nrow(rna_int[rna_int$beta_chron_age<0,])/nrow(rna_int)
nrow(rna_int[rna_int$beta_chron_age>0,])

eq2_volcano<- rna_int %>%
  ggplot(aes(beta_chron_age, -log10(pval_chron_age), colour = -log10(pval_chron_age) < -log10(0.05))) +
  geom_point(alpha = 0.5, size = 0.05) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_colour_manual(values = c('steelblue4', 'steelblue1')) +
  theme_classic(base_size = 7) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt")) +
  xlab("Beta Eq.1 Age") +
  ylab("-log10(P-value)")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq2_volcano.svg", 
       eq2_volcano, 
       height = 50, width = 50, units = "mm")

geom_text_repel(data = top20,
                aes(label = outcome),
                size = 3,
                max.overlaps = Inf,
                box.padding = 0.3,
                point.padding = 0.2,
                show.legend = FALSE)

nrow(rna_int[rna_int$beta_eq2_w<0,])/nrow(rna_int)
nrow(rna_int[rna_int$beta_eq2_w>0,])

eq3_volcano<- rna_int %>%
  ggplot(aes(beta_eq2_w, -log10(pval_eq2_w), colour = -log10(pval_eq2_w) < -log10(0.05))) +
  geom_point(alpha = 0.5, size = 0.05) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_colour_manual(values = c('green4', 'green1')) +
  theme_classic(base_size = 7) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt")) +
  xlab("Beta Eq.2 Age Within") +
  ylab("-log10(P-value)")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq3_volcano.svg", 
       eq3_volcano, 
       height = 50, width = 50, units = "mm")


#GSEA---------------------------------------------------------------------------
#Generate function
run_gsea<- function(g, b, df, geneset) {
  
  dat<- df %>%
    dplyr::select(c({{g}}, {{b}})) 
  
  colnames(dat)<- c("genes", "betas")
  
  dat<- dat %>% arrange(desc(betas))
  
  dat2<- dat$betas
  names(dat2) = dat$genes
  
  gsea_out<- fgsea(pathways = geneset, 
                   stats = dat2)
  
  gsea_out<- gsea_out %>%
    mutate(pathway = gsub("HALLMARK_", "", pathway),
           pathway = gsub("_", " ", pathway)) %>%
    arrange(desc(NES))
  
  return(gsea_out)
  
}

#Generate hallmark gene set
hallmark.msigdb = msigdbr(species = "Macaca mulatta", collection = "H")
hallmark_list = split(x = hallmark.msigdb$gene_symbol, f = hallmark.msigdb$gs_name)

chron_gsea<- run_gsea(outcome, beta_chron_age, rna_int, hallmark_list)
chron_gsea<- chron_gsea %>%
  mutate(pathway = gsub("_", " ", pathway))

gsea_chron_plot<- chron_gsea %>%
  ggplot(aes(x=reorder(pathway, NES), y=NES, fill = NES < 0)) +
  geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  scale_fill_manual(values = c('steelblue2', 'steelblue4')) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        axis.line = element_line(colour = "black", linewidth = 0.5)) +
  ylab("NES") +
  xlab("Hallmark Gene Set") +
  coord_flip()
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/gsea_chron.svg", 
       gsea_chron_plot, 
       height = 85, width = 85, units = "mm")

eq2_gsea<- run_gsea(outcome, beta_eq2_w, rna_int, hallmark_list)

gsea_eq2_plot<- eq2_gsea %>%
  ggplot(aes(x=reorder(pathway, NES), y=NES, fill = NES < 0)) +
  geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  scale_fill_manual(values = c('green3', 'green4')) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        axis.line = element_line(colour = "black", linewidth = 0.5)) +
  ylab("NES") +
  xlab("Hallmark Gene Set") +
  coord_flip()

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/gsea_eq2.svg", 
       gsea_eq2_plot, 
       height = 85, width = 85, units = "mm")

eq2_m_gsea<- run_gsea(outcome, beta_eq2_m, rna_int, hallmark_list)

eq2_m_gsea %>%
  ggplot(aes(x=reorder(pathway, NES), y=NES, fill = NES < 0)) +
  geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  scale_fill_manual(values = c('green3', 'green4')) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_line(colour = "black", linewidth = 0.5)) +
  ylab("NES") +
  xlab("Annotation") +
  coord_flip()

#Import promoter pqlseq files---------------------------------------------------
setwd('/scratch/ckelsey4/Cayo_meth/glmer_model_compare')
import_prom_pqlseq<- function(x){
  
  #Generate list of file names
  file_list<- list.files(pattern = x)
  file_order<- str_split_i(file_list, "_", 4)
  
  #Import glm models as list
  model_list<- lapply(file_list, readRDS)
  
  #Rename list elements
  names(model_list)<- file_order
  
  #Bind model list to df and add rownames
  model<- do.call(rbind, model_list)
  model$outcome2<- model$outcome
  model$outcome<- str_split_i(model$outcome, "\\.", 2)
  rownames(model)<- model$region
  
  #Separate region coordinates into start and end, delete the chr col, and move region col to front
  #model$outcome2<- model$outcome
  model<- model %>%
    separate_wider_delim(outcome2, ".", names = c("chrom", "gene")) %>%
    dplyr::select(-c(gene))
  model<- model %>%
    relocate(c(chrom), .before = n)
  #model<- model %>%
  #mutate(chromStart = as.numeric(chromStart))
  
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
    dplyr::select(-c(elapsed_time, converged))
  
}

#Age Within
prom_agew_files<- 'prom_pqlseq2_agew_'
prom_agew_pqlseq<- import_prom_pqlseq(prom_agew_files)

colnames(mm_genes)<- c("outcome", "gene_name")

#Assign gene names
prom_agew_pqlseq<- left_join(prom_agew_pqlseq, mm_genes, by = "outcome")
prom_agew_pqlseq<- prom_agew_pqlseq %>%
  mutate(across(4:11, as.numeric)) %>%
  relocate(gene_name, .after = outcome)
colnames(prom_agew_pqlseq)<- c("ensembl_name", "outcome", "chrom", "n", paste("dnam_", colnames(prom_agew_pqlseq[,5:12]), sep = ""))

age.w_full<- inner_join(rna_int, prom_agew_pqlseq, by = 'outcome')

age.w_full %>%
  ggplot(aes(dnam_beta, beta_eq2_w)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt")) +
  stat_cor(label.x = 0.1)

age.w_full %>%
  filter(pval_eq2_w < .05 & dnam_fdr < .05) %>%
  ggplot(aes(as.numeric(dnam_beta), beta_eq2_w)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt")) +
  stat_cor(label.x = -0.3, label.y = -0.8)

save.image("/home/ckelsey4/rna_data/rna_analysis.RData")
