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
  scale_color_gradient2(low = "purple", mid = "grey70", high = "green4", midpoint = 0, name = "") +
  theme_classic(base_size=24) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  xlab(expression(beta["Eq.2"])) +
  ylab(expression(beta["Eq.3"]))
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq2_eq3_scatterplot_rna.svg", 
       eq2_eq3_intercept, 
       height = 50, width = 50, units = "mm")

rna_slopes %>%
  ggplot(aes(beta_eq2_w, beta_eq3_age)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic(base_size=6) +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  stat_cor() +
  xlab(expression(beta["Eq.2"])) +
  ylab(expression(beta["Eq.3"]))

#PCA----------------------------------------------------------------------------
rna_pca<- prcomp(rna_norm, center = TRUE, scale. = TRUE)

pcs<- as.data.frame(rna_pca$x)

#Check which pca's explain the most variance
summary(rna_pca)$importance[2, ]

pcs<- cbind(pcs[1:3], base_meta)

pcs %>% ggplot(aes(x = pid, y = age_at_sampling)) +
  geom_boxplot()

pcs %>% ggplot(aes(x = prep_year, y = age_at_sampling)) +
  geom_boxplot()

pcs %>%
  ggplot(aes(within_age, PC1)) +
  geom_point() +
  geom_smooth(method = "lm")

pc.matrix<- model.matrix(~ PC1 + PC2 + trapped_age + within_age + mean_age + sex + Seq_batch + 
                           p_reads_trimmed + p_uniq_mapped + p_duplicates + p_gene_counts, data = pcs)
pc.matrix %>% 
  cor(use="pairwise.complete.obs") %>%
  ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2)

# Plot Model Outcomes ----------------------------------------------------------
#Effect sizes
compare_plot<- function(df, fdr1, fdr2, var1, var2, plot_type) {
  
  v1<-deparse(substitute(var1))
  v2<-deparse(substitute(var2))
  
  df<- df %>%
    filter({{fdr1}} < .05 | {{fdr2}} < .05) %>%
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
      ggplot(aes({{var1}}, {{var2}}, colour = diff)) +
      geom_point(size = 1, alpha = 0.8) +
      geom_abline() +
      geom_smooth(method = "lm", linewidth = 0.5) +
      geom_vline(xintercept=0, linetype="dashed") +
      geom_hline(yintercept=0, linetype="dashed") +
      theme_classic(base_size = 18) +
      theme(legend.position = "none") +
      theme(panel.background = element_rect(colour = "black", linewidth=1),
            axis.line = element_line(colour = "black", linewidth = 0.5),
            plot.margin = margin(1, 1, 1, 1, "pt"),
            aspect.ratio = 1,
            panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
            panel.grid.minor = element_line(color = "grey98", linewidth = 0.5))
    
  } else if (plot_type == "hist") {
    
    df %>%
      ggplot(aes(diff, fill = after_stat(x))) +
      geom_histogram(bins = 50) +
      geom_vline(xintercept=0, linetype="dashed") +
      geom_vline(xintercept=median(df$diff), linetype="dashed", colour = 'red') +
      theme_classic(base_size = 18) +
      theme(legend.position = "none") +
            #legend.key.width = unit(1, 'mm'), 
            #legend.key.height = unit(5, 'mm')) +
      theme(panel.background = element_rect(colour = "black", linewidth=1),
            axis.line = element_line(colour = "black", linewidth = 0.5),
            plot.margin = margin(1, 1, 1, 1, "pt"),
            aspect.ratio = 1,
            panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
            panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) 
    
  }
}

within_chron_plot_rna<- compare_plot(rna_int, pval_chron_age, pval_eq3_age, 
                                     beta_chron_age, beta_eq3_age, "scatter")

within_chron_plot_rna<- within_chron_plot_rna +
  scale_color_gradient2(low = "steelblue2", mid = "grey70", high = "purple", midpoint = 0, 
                        name = "") +
  xlab(expression(beta["Eq.1"])) +
  ylab(expression(beta["Eq.3"]))  +
  xlim(-1.0, 1.0) +
  ylim(-1.0, 1.0) 
within_chron_plot_rna

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/within_chron_scatter_rna.svg", 
       within_chron_plot_rna, 
       height = 50, width = 50, units = "mm")

within_chron_hist_rna<- compare_plot(rna_int, pval_chron_age, pval_eq3_age, 
                                     beta_chron_age, beta_eq3_age,"hist")
within_chron_hist_rna<- within_chron_hist_rna +
  scale_fill_gradient2(low = "steelblue2", mid = "grey70", high = "purple", midpoint = 0, name = "") +
  xlab(expression(beta["Eq.3"] - beta["Eq.1"])) +
  ylab("Count")
within_chron_hist_rna
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/within_chron_hist_rna.svg", 
       within_chron_hist_rna, 
       height = 50, width = 50, units = "mm")

compare_plot(rna_int, pval_eq2_w, pval_eq2_m,
             beta_eq2_w, beta_eq2_m, "scatter")

compare_plot(rna_int, pval_eq2_, pval_eq3_age,
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
  age.m.count<- nrow(rna[rna$pval_eq2_m < 0.05,])
  total<- data.frame(count = c(age.w.count, age.eq3.count, age.chron.count, age.m.count),
                      predictor = as.factor(c('Eq.2 W.', 'Eq.3', 'Eq.1', 'Eq.2 B.')))
  
  total$predictor<- factor(total$predictor, levels = c('Eq.3', 'Eq.2 W.','Eq.1', "Eq.2 B."))
  total<- total %>%
    mutate(perc_signif = count/nrow(rna))
  
  return(total)
  
}

int_counts<- generate_counts(rna_int)
slope_counts<- generate_counts(rna_slopes)

signif_counts<- int_counts %>%
  ggplot(aes(predictor, count, fill = predictor)) +
  geom_bar(stat = 'identity') +
  geom_text(label=int_counts$count, vjust=-0.25, size=1.5) +
  theme_classic(base_size=6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt")) +
  xlab("Model") +
  ylab("# Significant DEGs") +
  scale_fill_manual(values = c("purple", 'green4', 'steelblue2', 'grey30'))
signif_counts
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/signif_counts_rna.svg", 
       signif_counts, 
       height = 50, width = 50, units = "mm")

beta_dist_rna<- rna_int %>%
  dplyr::select(c(beta_eq3_age, beta_eq2_w, beta_chron_age, beta_eq2_m)) %>%
  pivot_longer(cols = c(beta_eq3_age, beta_eq2_w, beta_chron_age, beta_eq2_m),
               values_to = 'beta',
               names_to = 'var') %>%
  mutate(var = factor(var, levels = rev(c("beta_eq2_m", "beta_chron_age",
                                          "beta_eq2_w", "beta_eq3_age")))) %>%
  ggplot(aes(var, beta, fill=var)) +
  geom_violin() +
  geom_boxplot(width = 0.05, fill = "white", outlier.size = 0.25) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_fill_manual(values = c("purple", 'green4', 'steelblue2', 'grey30')) +
  theme_classic(base_size=6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt")) +
  scale_x_discrete(labels = c('Eq.3', 'Eq.2 W.','Eq.1', "Eq.2 B.")) +
  xlab("Model") +
  ylab("Beta")
beta_dist_rna
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/beta_dist_rna.svg", 
       beta_dist_rna, 
       height = 50, width = 50, units = "mm")

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
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1) +
  xlab(expression(beta["Eq.1"])) +
  ylab("-log10(P-value)")
eq2_volcano
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
  ggplot(aes(beta_eq3_age, -log10(pval_eq3_age), colour = -log10(pval_eq3_age) < -log10(0.05))) +
  geom_point(alpha = 0.5, size = 0.05) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_colour_manual(values = c("purple","purple4")) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1) +
  xlab(expression(beta["Eq.3"])) +
  ylab("-log10(P-value)")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq3_volcano.svg", 
       eq3_volcano, 
       height = 50, width = 50, units = "mm")

eq2_volcano_m<- rna_int %>%
  ggplot(aes(beta_eq2_m, -log10(pval_eq2_m), colour = -log10(pval_eq2_m) < -log10(0.05))) +
  geom_point(alpha = 0.5, size = 0.05) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  #scale_colour_manual(values = c('green4', 'green1')) +
  theme_classic(base_size = 7) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1) +
  xlab("Beta Eq.2 Age-Between") +
  ylab("-log10(P-value)")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq2_m_volcano_rna.svg", 
       eq2_volcano_m, 
       height = 50, width = 50, units = "mm")

eq3_volcano_m<- rna_int %>%
  ggplot(aes(beta_eq3_m, -log10(pval_eq3_m), colour = -log10(pval_eq3_m) < -log10(0.05))) +
  geom_point(alpha = 0.5, size = 0.05) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  #scale_colour_manual(values = c('green4', 'green1')) +
  theme_classic(base_size = 7) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1) +
  xlab("Beta Eq.3 Age-Between") +
  ylab("-log10(P-value)")
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq3_m_volcano_rna.svg", 
       eq3_volcano_m, 
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

eq3_gsea<- run_gsea(outcome, beta_eq3_age, rna_int, hallmark_list)

gsea_eq3_plot<- eq3_gsea %>%
  ggplot(aes(x=reorder(pathway, NES), y=NES)) +
  geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black" , fill = 'purple') +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        axis.line = element_line(colour = "black", linewidth = 0.5)) +
  ylab("NES") +
  xlab("Hallmark Gene Set") +
  coord_flip()
gsea_eq3_plot
ggsave("/home/ckelsey4/Cayo_meth/aging_plots/gsea_eq3.svg", 
       gsea_eq3_plot, 
       height = 85, width = 85, units = "mm")

eq2_m_gsea<- run_gsea(outcome, beta_eq2_m, rna_int, hallmark_list)

eq2_m_gsea %>%
  ggplot(aes(x=reorder(pathway, NES), y=NES, fill = padj > .05)) +
  geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        axis.line = element_line(colour = "black", linewidth = 0.5)) +
  ylab("NES") +
  xlab("Hallmark Gene Set") +
  coord_flip()

eq2_m_gsea %>%
  ggplot(aes(x=reorder(pathway, NES), y=NES, fill = padj > .05)) +
  geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        axis.line = element_line(colour = "black", linewidth = 0.5)) +
  ylab("NES") +
  xlab("Hallmark Gene Set") +
  coord_flip()

full_gsea<- inner_join(chron_gsea[,1:6], eq2_gsea[,1:6], suffix = c("_eq1", "_eq2_w"), by = "pathway")
full_gsea<- inner_join(full_gsea, eq2_m_gsea[,1:6], by = "pathway")
colnames(full_gsea)[12:16] <- c(paste(colnames(full_gsea)[12:16], "_eq2_m", sep = ""))

full_gsea %>%
  ggplot(aes(NES_eq2_m, NES_eq1, shape = padj_eq1 < .05, colour = padj_eq2_w < .05)) +
  geom_point(size =3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic(base_size = 7) +
  theme(
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1) +
  xlab("NES Eq.2 Between-Age") +
  ylab("NES Eq.1")

#Import promoter pqlseq files---------------------------------------------------
#Import nested pqlseq model-----------------------------------------------------
setwd('/scratch/ckelsey4/Cayo_meth/glmer_model_compare')

#Define import function
import_prom_models<- function(model_path, file_string, mod_type){
  
  file_list<- list.files(path = model_path, pattern = file_string)
  file_order<- str_split_i(file_list, "_", 5)
  
  #Import glm models as list
  model_list<- lapply(paste(model_path, file_list, sep = "/"), readRDS)
  
  #Rename list elements
  names(model_list)<- file_order
  
  if (mod_type == "eq1") {
    
    df_nms<- c(mod_type)
    
  } else if (mod_type == "eq2") {
    
    df_nms<- c("eq2.w", "eq2.m")
    
  } else if (mod_type == "eq3") {
    
    df_nms<- c("eq3_age", "eq3.m")
    
  }
  
  mod2<- lapply(model_list, function(mod){
    
    names(mod)<- df_nms
    
    mod<- lapply(names(mod), function(y){
      
      df<- mod[[y]]
      df<- df %>%
        dplyr::select(-c(h2, sigma2))
      
      colnames(df)<- c("outcome", "n", paste(colnames(df[3:length(df)]), y, sep = "_"))
      
      df
    })
    
    names(mod)<- df_nms
    
    mod
    
  })
  
  mod2 <- lapply(df_nms, function(nm) {
    do.call(rbind, lapply(mod2, function(x) x[[nm]]))
  })
  
  names(mod2)<- df_nms
  
  mod2<- lapply(mod2, function(df){
    
    df<- df %>%
      separate_wider_delim(outcome, delim = ".", names = c("chr", "outcome"))
    
  })
  
  return(mod2)
  
}

#Import female and male pqlseq output
eq1_prom_list<- import_prom_models(model_path = "/home/ckelsey4/age_and_aging/models_out",
                                   file_string = "dnam_prom_eq1_model", mod_type = "eq1")

eq2_prom_list<- import_prom_models(model_path = "/home/ckelsey4/age_and_aging/models_out",
                                   file_string = "dnam_prom_eq2_model", mod_type = "eq2")

eq3_prom_list<- import_prom_models(model_path = "/home/ckelsey4/age_and_aging/models_out",
                                   file_string = "dnam_prom_eq3_model", mod_type = "eq3")

prom_df<- left_join(eq1_prom_list[[1]], eq2_prom_list[[1]][,c(2,4:9)], by = "outcome")
prom_df<- left_join(prom_df, eq2_prom_list[[2]][,c(2,4:9)], by = "outcome")
prom_df<- left_join(prom_df, eq3_prom_list[[1]][,c(2,4:9)], by = "outcome")
prom_df<- left_join(prom_df, eq3_prom_list[[2]][,c(2,4:9)], by = "outcome")

rm(eq1_prom_list);rm(eq2_prom_list);rm(eq3_prom_list)

#Assign gene names
prom_df<- left_join(prom_df, mm_genes, by = "outcome")
prom_df<- prom_df %>%
  relocate(gene_name, .after = outcome)
colnames(prom_df)<- c("chr", "ensembl_name", "outcome", "n", 
                      paste("dnam_", colnames(prom_df[,5:length(prom_df)]), 
                            sep = ""))

prom_df<- prom_df %>%
  mutate(eq3_chron_diff = abs(dnam_beta_eq1) - abs(dnam_beta_eq3_age))

within_chron_plot_dnam<- compare_plot(prom_df, dnam_fdr_eq1, dnam_fdr_eq3_age, 
                                      dnam_beta_eq1, dnam_beta_eq3_age, "scatter")
within_chron_plot_dnam<- within_chron_plot_dnam +
  scale_colour_gradient2(low = "steelblue2", mid = "grey70", high = "purple", midpoint = 0, name = "") +
  xlim(-0.20, 0.20) +
  ylim(-0.20, 0.20)

within_chron_plot_dnam

dna_rna<- inner_join(rna_int, prom_df, by = 'outcome')

dna_rna %>%
  #filter(pval_eq2_w < .20 & dnam_fdr_eq2.w < .20) %>%
  ggplot(aes(dnam_beta_eq2.w, beta_eq2_w)) +
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

dna_rna %>%
  #filter(pval_chron_age < .20 & dnam_fdr_eq1 < .20) %>%
  ggplot(aes(dnam_beta_eq1, beta_chron_age)) +
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
  stat_cor(label.x = 0.05)

#FDR Thresholds
thresholds <- c(0.01, 0.05, 0.10, 0.15, 0.20, 1.0)

results_eq1<- sapply(thresholds, function(i) {
  
  mod<- dna_rna %>%
    filter(pval_chron_age < i & dnam_fdr_eq1 < i)
  
  res<- cor(mod$beta_chron_age, mod$dnam_beta_eq1)
  
  res
})

results_eq2<- sapply(thresholds, function(i) {
  
  mod<- dna_rna %>%
    filter(pval_eq2_w < i & dnam_fdr_eq2.w < i)
  
  res<- cor(mod$beta_eq2_w, mod$dnam_beta_eq2.w)
  
  res
})

save.image("/home/ckelsey4/rna_data/rna_analysis.RData")
