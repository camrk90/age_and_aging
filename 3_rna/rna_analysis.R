library(tidyverse)
library(ggplot2)
library(ggeffects)
library(ggcorrplot)
library(ggpubr)
library(ggrepel)
library(fgsea)
library(msigdbr)
library(Biostrings)
library(biomaRt)
library(ggvenn)
library(UpSetR)
library(lme4)
library(limma)
library(edgeR)
library(EMMREML)
library(variancePartition)

#load("/home/ckelsey4/rna_data/rna_analysis.RData")
#Use for local
parent_dir<- paste0(getwd(), "/local_data/")

#Use for remote (SOL)
parent_dir<- "/home/ckelsey4/"

#Load data
eq1<- readRDS(paste0(parent_dir, "rna_eq1"))
eq2<- readRDS(paste0(parent_dir, "rna_eq2"))
eq3<- readRDS(paste0(parent_dir, "rna_eq3"))
base_meta<- read.table(paste0(parent_dir, "base_meta.txt"))
rna_counts<- readRDS(paste0(parent_dir, "Cayo_PBMC_longLPS_counts_9Jan26.rds"))

#Make simplified outcome df
eq1_int<- eq1[["df"]]
eq2_int<- eq2[["df"]]
eq3_int<- eq3[["df"]]
rna_int<- as.data.frame(cbind(eq1_int$outcome, eq1_int$beta_trapped_age, 
                              eq1_int$pvalue_trapped_age, eq1_int$se_trapped_age,
                              eq2_int$beta_within_age, eq2_int$pvalue_within_age,
                              eq2_int$se_within_age, eq2_int$beta_mean_age, 
                              eq2_int$pvalue_mean_age, eq2_int$se_mean_age,
                              eq3_int$beta_trapped_age, eq3_int$pvalue_trapped_age,
                              eq3_int$se_trapped_age, eq3_int$beta_mean_age, 
                              eq3_int$pvalue_mean_age, eq3_int$se_mean_age))

c_names<- c("outcome", "beta_chron_age", "pval_chron_age", "se_chron_age",
            "beta_eq2_w", "pval_eq2_w", "se_eq2_w",
            "beta_eq2_m", "pval_eq2_m", "se_eq2_m",
            "beta_eq3_age", "pval_eq3_age", "se_eq3_age",
            "beta_eq3_m", "pval_eq3_m", "se_eq3_m")

colnames(rna_int)<- c_names

rna_int<- rna_int %>%
  mutate(across(2:16, as.numeric))

#Variance Partition
plotVarPart(eq1[["vp"]])
plotVarPart(eq2[["vp"]])
plotVarPart(eq3[["vp"]])

#PCA----------------------------------------------------------------------------
rna_counts<- rna_counts[, base_meta$Sample_ID]
rna_norm<- voom(calcNormFactors(DGEList(counts=rna_counts)), plot=FALSE)
rna_norm<- rna_norm[["E"]]
colnames(rna_norm)<- colnames(rna_counts)
rownames(rna_norm)<- rownames(rna_counts)
rna_norm<- t(rna_norm)

rna_norm<- rna_norm[rownames(rna_norm) %in% base_meta$Sample_ID,]

rna_pca<- prcomp(rna_norm, center = TRUE, scale. = TRUE)

pcs<- as.data.frame(rna_pca$x)

#Check which pca's explain the most variance
summary(rna_pca)$importance[2, ]

pcs<- cbind(pcs[1:5], base_meta)

pc.matrix<- model.matrix(~ PC1 + PC2 + trapped_age + within_age + mean_age + sex + Seq_batch + 
                           p_reads_trimmed + p_uniq_mapped + p_duplicates + p_gene_counts,
                         data = pcs)
pc.matrix %>% 
  cor(use="pairwise.complete.obs") %>%
  ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2)

# Plot Model Outcomes ----------------------------------------------------------
### Counts for significant genes------------------------------------------------
count_signif_regions<- function(x) {
  
  #Generate count, proportion, and percentage of significant genes
  df<- x %>%
    dplyr::select(starts_with("pval_"))
  
  vars<- gsub("pval_", "", colnames(df))
  
  counts <- colSums(df < 0.05, na.rm = TRUE)
  
  counts<- data.frame(predictor = vars,
                      count = counts)
  
  counts<- counts %>%
    arrange(counts) %>%
    mutate(predictor = factor(predictor, levels = predictor))
  
  counts<- counts %>%
    mutate(proportion_signif = count/nrow(df),
           perc_signif = proportion_signif*100)
  
  #Generate comparison with Eq.1 Age beta values
  test<- lapply(x[,startsWith(colnames(x), "beta_")], function(y){
    
    dfr<- cor.test(y, rna_int$beta_chron_age)["estimate"]
    dfm<- median(y - rna_int$beta_chron_age)
    
    dd<- data.frame(correlation = dfr, median_diff = dfm)
    return(dd)
  })
  
  test<- as.data.frame(do.call(rbind, test))
  
  test$predictor<- gsub("beta_", "", rownames(test))
  
  #Generate significance comparison
  signif<- lapply(x[,startsWith(colnames(x), "pval_")], function(y){
    
    shared<- x %>% filter(y < .05 | pval_chron_age < .05) %>% nrow()
    
  })
  
  signif<- as.data.frame(do.call(rbind, signif))
  colnames(signif)<- "shared_genes"
  
  signif$predictor<- gsub("pval_", "", rownames(signif))
  
  #Combine outputs
  counts<- left_join(counts, test, by = "predictor")
  counts<- left_join(counts, signif, by = "predictor")
  
  #Plot counts
  counts_plot<- counts %>%
    filter(!predictor %in% c("eq2_w", "eq3_m")) %>%
    ggplot(aes(reorder(predictor, count), count, fill = predictor)) +
    geom_bar(stat = 'identity') +
    geom_text(label=counts$count[!counts$predictor %in% c("eq2_w", "eq3_m")], vjust=-0.25, size = 3) +
    theme_classic(base_size = 6) +
    theme(legend.position = "none",
          panel.background = element_rect(colour = "black", linewidth=1),
          axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.title.x = element_blank(),
          plot.margin = margin(1, 1, 1, 1, "pt")) +
    ylab("# Significant Regions")
  
  return(list(plot = counts_plot, df = counts))
  
}

counts<- count_signif_regions(rna_int)
counts[["plot"]] +
  scale_fill_manual(values = c("steelblue2", "grey30", "purple")) +
  scale_x_discrete(labels = c("chron_age" = "Eq.1", 
                              "eq2_m" = "Eq.2 Btwn", "eq3_age" = "Eq.3 Within"))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/signif_counts.svg",
       height = 65, width = 65, units = "mm")

### Effect sizes----------------------------------------------------------------
# Plot distribution of effect sizes for models
## Eq2 within is removed, x-axis is limited to 0+/-0.5, 201 Eq.3 genes are lost
rna_int %>%
  dplyr::select(c(beta_eq3_age, beta_eq2_w, beta_chron_age, beta_eq2_m)) %>%
  pivot_longer(cols = c(beta_eq3_age, beta_eq2_w, beta_chron_age, beta_eq2_m),
               values_to = 'beta',
               names_to = 'var') %>%
  mutate(var = factor(var, levels = c("beta_eq2_m", "beta_chron_age",
                                      "beta_eq2_w", "beta_eq3_age"))) %>% 
  filter(var %in% c("beta_chron_age", "beta_eq2_m", "beta_eq3_age")) %>%
  ggplot(aes(beta, fill=var)) +
  #geom_violin() +
  #geom_boxplot(width = 0.05, fill = "white", outlier.size = 0.25) +
  geom_density(alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 'dashed', colour = "red") +
  scale_fill_manual(values = rev(c("purple",'steelblue2', 'grey30'))) +
  theme_classic(base_size=6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt")) +
  #scale_x_discrete(labels = c('Eq.3', 'Eq.2 W.','Eq.1', "Eq.2 B.")) +
  xlim(-0.5, 0.5) +
  xlab(expression(beta)) +
  ylab("Density")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/beta_dist_rna.svg",
       height = 50, width = 50, units = "mm")

#Generate plot function
compare_plot<- function(df, fdr1, fdr2, var1, var2, plot_type) {
  
  v1<-deparse(substitute(var1))
  v2<-deparse(substitute(var2))
  
  df<- df %>%
    filter({{fdr1}} < .05 | {{fdr2}} < .05) %>%
    mutate(diff = abs({{var2}}) - abs({{var1}}),
           ratio = abs({{var2}})/abs({{var1}}))
  
  if (plot_type == "scatter"){
    
    df %>%
      ggplot(aes({{var1}}, {{var2}}, colour = diff)) +
      geom_point(size = 1, alpha = 0.8) +
      geom_abline() +
      geom_smooth(method = "lm", linewidth = 0.5) +
      geom_vline(xintercept=0, linetype="dashed") +
      geom_hline(yintercept=0, linetype="dashed") +
      theme_classic(base_size = 6) +
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
      theme_classic(base_size = 6) +
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

#Eq.2 Age Within vs Eq.3 Age
compare_plot(rna_int, pval_eq2_w, pval_eq3_age,
             beta_eq2_w, beta_eq3_age, "scatter") +
  scale_fill_gradient2(low = "green4", mid = "grey70", high = "purple", 
                       midpoint = 0, name = "") +
  xlab(expression(beta["Eq.2"])) +
  ylab(expression(beta["Eq.3"]))  +
  xlim(-1.0, 1.0) +
  ylim(-1.0, 1.0) 

# Eq.1 Age vs Eq.3 Age
## Scatterplot
compare_plot(rna_int, pval_chron_age, pval_eq3_age, 
             beta_chron_age, beta_eq3_age, "scatter") +
  scale_color_gradient2(low = "steelblue2", mid = "grey70", high = "purple", 
                        midpoint = 0, name = "") +
  xlab(expression(beta["Eq.1"])) +
  ylab(expression(beta["Eq.3"]))  +
  xlim(-1.0, 1.0) +
  ylim(-1.0, 1.0)

ggsave(paste0(parent_dir, "plots", "within_chron_scatter_rna.svg"),
       height = 50, width = 50, units = "mm")

## Histogram
compare_plot(rna_int, pval_chron_age, pval_eq3_age, 
             beta_chron_age, beta_eq3_age,"hist") +
  scale_fill_gradient2(low = "steelblue2", mid = "grey70", high = "purple", 
                       midpoint = 0, name = "") +
  xlab(expression(beta["Eq.3"] - beta["Eq.1"])) +
  ylab("Count")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/within_chron_hist_rna.svg", 
       height = 50, width = 50, units = "mm")

# Eq.2 Within vs Eq.2 Between
compare_plot(rna_int, pval_eq2_w, pval_eq2_m, 
             beta_eq2_w, beta_eq3_age,"scatter") +
  scale_fill_gradient2(low = "steelblue2", mid = "grey70", high = "purple", 
                       midpoint = 0, name = "") +
  xlab(expression(beta["Eq.3"] - beta["Eq.1"])) +
  ylab("Count")

# Eq.1 Age vs Eq.2 Between Age
compare_plot(rna_int, pval_chron_age, pval_eq2_m, 
             beta_chron_age, beta_eq2_m, "scatter") +
  scale_colour_gradient2(low = "steelblue2", mid = "grey70", high = "grey30", 
                         midpoint = 0, name = "") +
  xlab(expression(beta["Eq.2 Between"])) +
  ylab(expression(beta["Eq.3"]))  +
  xlim(-0.1, 0.1) +
  ylim(-0.1, 0.1) 

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

rna_int$outcome[match(mm_genes2$anno, rna_int$outcome)]<- mm_genes2$gene_name

rm(rna_genes);rm(rna_genes2)

top20<- rna_int %>%
  arrange(desc(beta_chron_age))

top20<- top20[c(1:10, (nrow(top20)-9):nrow(top20)), ]

nrow(rna_int[rna_int$beta_chron_age<0,])/nrow(rna_int)
nrow(rna_int[rna_int$beta_chron_age>0,])

rna_int %>%
  ggplot(aes(beta_chron_age, -log10(pval_chron_age), colour = -log10(pval_chron_age) < -log10(0.05))) +
  geom_point(alpha = 0.5, size = 0.05) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_colour_manual(values = c('steelblue4', 'steelblue1')) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1) +
  xlab(expression(beta["Eq.1"])) +
  ylab("-log10(P-value)")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq1_volcano.svg", 
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

rna_int %>%
  ggplot(aes(beta_eq3_age, -log10(pval_eq3_age), colour = -log10(pval_eq3_age) < -log10(0.05))) +
  geom_point(alpha = 0.5, size = 0.05) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_colour_manual(values = c("purple","purple4")) +
  theme_classic(base_size = 18) +
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
       height = 50, width = 50, units = "mm")

rna_int %>%
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

chron_gsea %>%
  filter(abs(NES) > 1) %>%
  ggplot(aes(x=reorder(pathway, NES), y=NES)) +
  #geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  geom_point(aes(alpha=padj<0.05, size = size), colour = "steelblue2") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("NES") +
  xlab("Hallmark Gene Set") +
  coord_flip()

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/gsea_chron.svg", 
       height = 85, width = 85, units = "mm")

eq3_gsea<- run_gsea(outcome, beta_eq3_age, rna_int, hallmark_list)

eq3_gsea %>%
  filter(abs(NES) > 1) %>%
  ggplot(aes(x=reorder(pathway, NES), y=NES)) +
  #geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black" , fill = 'purple') +
  geom_point(aes(alpha=padj<0.05, size = size), colour = 'purple') +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("NES") +
  xlab("Hallmark Gene Set") +
  coord_flip()

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/gsea_eq3.svg", 
       gsea_eq3_plot, 
       height = 85, width = 85, units = "mm")

eq2_m_gsea<- run_gsea(outcome, beta_eq2_m, rna_int, hallmark_list)

eq2_m_gsea %>%
  filter(abs(NES) > 1) %>%
  ggplot(aes(x=reorder(pathway, NES), y=NES)) +
  #geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black", fill = "grey30") +
  geom_point(aes(alpha=padj<0.05, size = size), colour = "grey30") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
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
long_data<- read.table("/scratch/ckelsey4/Cayo_meth/long_data_adjusted.txt")

prom_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/prom_cov_filtered")
prom_cov<- do.call(rbind, prom_cov[1:21])
rownames(prom_cov)<- str_split_i(rownames(prom_cov), "\\.", 3)

prom_m<- readRDS("/scratch/ckelsey4/Cayo_meth/prom_m_filtered")
prom_m<- do.call(rbind, prom_m[1:21])
rownames(prom_m)<- str_split_i(rownames(prom_m), "\\.", 3)

prom_examples<- c("ENSMMUG00000033572", "ENSMMUG00000060797", "ENSMMUG00000044622", "ENSMMUG00000013067")

prom_cov<- as.data.frame(t(prom_cov[rownames(prom_cov) %in% prom_examples, ]))
prom_cov$lid_pid<- rownames(prom_cov)
long_data<- inner_join(long_data, prom_cov, by = "lid_pid")

prom_m<- as.data.frame(t(prom_m[rownames(prom_m) %in% prom_examples, ]))
prom_m$lid_pid<- rownames(prom_m)
long_data<- inner_join(long_data, prom_m, by = "lid_pid", suffix = c("_cov", "_meth"))

prom_glm<- lapply(setNames(prom_examples, prom_examples), function(x){
  
  cov<- paste0(x, "_cov")
  meth<- paste0(x, "_meth")
  
  mod<- as.formula(paste0("cbind(`", meth, "`, `", cov, "`) ~ ",
                            "age_at_sampling + mean.age + individual_sex + university + (1|monkey_id)"))
  
  glmer(mod, data = long_data, family = binomial(link = "logit"))
  
})



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
  mutate(abs_diff = abs(beta_eq1) - abs(beta_eq3_age),
         direction = ifelse(abs_diff < 0, "Eq.3 Steeper", "Eq.1 Steeper"))

genes_to_remove<- c("Metazoa_SRP", "U1", "U2", "U3", "U4", "U5", "U6", "U7", "U8", "", "Y_RNA")
prom_df<- prom_df[!prom_df$gene_name %in% genes_to_remove, ]
prom_df<- prom_df %>% filter(!grepl("*_rRNA", gene_name))
prom_df<- prom_df %>% filter(!grepl("mml-mir-*", gene_name))
df2<- prom_df[prom_df$gene_name %in% genes_to_remove[2:9], ]

prom_df<- prom_df %>%
  mutate(fdr_eq1 = p.adjust(pvalue_eq1, method = "fdr"),
         fdr_eq2.m = p.adjust(pvalue_eq2.m, method = "fdr"),
         fdr_eq3_age = p.adjust(pvalue_eq3_age, method = "fdr"))

df<- prom_df %>%
  dplyr::select(pvalue_eq1,pvalue_eq3_age) %>%
  pivot_longer(cols = c(pvalue_eq1, pvalue_eq3_age))

df %>%
  ggplot(aes(value, fill = name)) + 
  geom_histogram(breaks = seq(0, 1, 0.01), colour = "black", alpha = 0.5, position = "identity") + 
  theme_classic(base_size = 6) +
  theme(
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        #aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.01)) +
  ylab("Count") +
  xlab("P-value")
  
prom_df %>%
  ggplot(aes(pvalue_eq1)) +
  geom_histogram(breaks = seq(0, 1, 0.05), colour = "black", fill = "steelblue1") +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  scale_x_continuous(breaks = seq(0, 1, 0.05)) +
  ylab("Count") +
  xlab("P-value Eq.1")

prom_df %>%
  ggplot(aes(pvalue_eq3_age)) +
  geom_histogram(breaks = seq(0, 1, 0.05), colour = "black", fill = "purple") +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  scale_x_continuous(breaks = seq(0, 1, 0.05)) +
  ylab("Count") +
  xlab("P-value Eq.3")
  

## Make Upset Plot
chron.age<- prom_df$outcome[prom_df$fdr_eq1 < 0.05]
eq2.btwn<- prom_df$outcome[prom_df$fdr_eq2.m < 0.05]
eq3<- prom_df$outcome[prom_df$fdr_eq3_age < 0.05]

venn_all<- list(chron.age, eq3, eq2.btwn)
names(venn_all)<- c("Eq.1", "Eq.3 Within", "Eq.2 Btwn")

upset(fromList(venn_all), order.by = "freq", 
      line.size = 0.5, point.size = 1)

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/dnam_proms_upset.svg", 
       height = 50, width = 50, units = "mm")


prom_df %>%
  ggplot(aes(direction, fill = direction)) +
  geom_bar(colour = 'black') +
  scale_fill_manual(values = c("steelblue1", "purple1")) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("Count") +
  xlab("Direction")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/dnam_proms_larger.svg", 
       height = 50, width = 50, units = "mm")

prom_df$eq3_signif<- "Neither"
prom_df$eq3_signif[prom_df$fdr_eq1 < .05 & prom_df$fdr_eq3_age > .05]<- "Eq.1 Signif."
prom_df$eq3_signif[prom_df$fdr_eq1 > .05 & prom_df$fdr_eq3_age < .05]<- "Eq.3 Signif."
prom_df$eq3_signif[prom_df$fdr_eq1 < .05 & prom_df$fdr_eq3_age < .05]<- "Both Signif."

prom_df %>%
  filter(fdr_eq1 < .05 | fdr_eq3_age < .05) %>%
  ggplot(aes(beta_eq1, beta_eq3_age)) +
  geom_point(aes(colour = eq3_signif)) +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_abline() +
  scale_colour_manual(values = c("black", "steelblue1", "purple1")) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  scale_y_continuous(breaks = seq(-0.15, 0.15, 0.05), limits = c(-0.15, 0.15), labels = scales::comma) +
  scale_x_continuous(breaks = seq(-0.15, 0.15, 0.05), limits = c(-0.15, 0.15), labels = scales::comma) +
  xlab(expression(beta["Eq.1"])) +
  ylab(expression(beta["Eq.3"]))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/dnam_proms_scatter.svg", 
       height = 50, width = 50, units = "mm")

prom_df$eq2_signif<- "Neither"
prom_df$eq2_signif[prom_df$fdr_eq1 < .05 & prom_df$fdr_eq2.m > .05]<- "Eq.1 Signif."
prom_df$eq2_signif[prom_df$fdr_eq1 > .05 & prom_df$fdr_eq2.m < .05]<- "Eq.2 Between Signif."
prom_df$eq2_signif[prom_df$fdr_eq1 < .05 & prom_df$fdr_eq2.m < .05]<- "Both Signif."

prom_df %>%
  filter(fdr_eq1 < .05 | fdr_eq2.m < .05) %>%
  ggplot(aes(beta_eq1, beta_eq2.m)) +
  geom_point(aes(colour = eq2_signif)) +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_abline() +
  scale_colour_manual(values = c("black", "steelblue1", "grey80")) +
  theme_classic(base_size = 6) +
  theme(
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  scale_y_continuous(breaks = seq(-0.10, 0.10, 0.05), limits = c(-0.10, 0.10), labels = scales::comma) +
  scale_x_continuous(breaks = seq(-0.10, 0.10, 0.05), limits = c(-0.10, 0.10), labels = scales::comma) +
  xlab(expression(beta["Eq.1"])) +
  ylab(expression(beta["Eq.2 Between"]))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/dnam_proms_scatter_eq2.svg", 
       height = 50, width = 50, units = "mm")

prom_df %>%
  ggplot(aes(abs_diff, fill = abs_diff < 0)) +
  geom_density() +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  scale_fill_manual(values = c("steelblue1", "purple1")) +
  ylab("Density") +
  xlab(expression(abs(beta["Eq.1"]) - abs(beta["Eq.3"])))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/dnam_proms_hist.svg", 
       height = 50, width = 50, units = "mm")

prom_df %>%
  filter(fdr_eq1 < .05 | fdr_eq2.m < .05) %>%
  ggplot(aes(beta_eq1, beta_eq2.m)) +
  geom_point(aes(colour = eq2_signif)) +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_abline() +
  scale_colour_manual(values = c("black", "steelblue1", "grey80")) +
  theme_classic(base_size = 6) +
  theme(
    panel.background = element_rect(colour = "black", linewidth=1),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    plot.margin = margin(1, 1, 1, 1, "pt"),
    aspect.ratio = 1,
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  scale_y_continuous(breaks = seq(-0.10, 0.10, 0.05), limits = c(-0.10, 0.10), labels = scales::comma) +
  scale_x_continuous(breaks = seq(-0.10, 0.10, 0.05), limits = c(-0.10, 0.10), labels = scales::comma) +
  xlab(expression(beta["Eq.1"])) +
  ylab(expression(beta["Eq.2 Between"]))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/dnam_proms_scatter_eq2.svg", 
       height = 50, width = 50, units = "mm")



genes_to_remove<- c("Metazoa_SRP", "U1", "U2", "U3", "U4", "U5", "U6", "U7", "U8", "", "Y_RNA")
df<- prom_df[!prom_df$gene_name %in% genes_to_remove, ]
df<- df %>% filter(!grepl("*_rRNA", gene_name))
df<- df %>% filter(!grepl("mml-mir-*", gene_name))
df2<- prom_df[prom_df$gene_name %in% genes_to_remove[2:9], ]

prom_eq3_signif<- df %>%
  filter(fdr_eq1 > .05 & fdr_eq3_age < .05)

top10<- df %>%
  arrange(abs_diff) %>%
  dplyr::slice(c(1:10, (n() - 9):n()))

top10 %>%
  dplyr::select(beta_eq1, beta_eq3_age, se_beta_eq1, se_beta_eq3_age, gene_name, abs_diff) %>%
  dplyr::rename(se_eq1 = se_beta_eq1, se_eq3_age = se_beta_eq3_age) %>%
  pivot_longer(cols = c(beta_eq1, beta_eq3_age, se_eq1, se_eq3_age),
               names_to = c(".value", "model"),
               names_sep = "_") %>%
  ggplot(aes(x=beta, y=reorder(gene_name, abs(abs_diff)), colour = model)) +
  geom_point(aes()) +
  #scale_size(range = c(0.05, 2)) +
  geom_path(aes(group = gene_name), colour = "black") +
  scale_colour_manual(values = c("steelblue2", "purple")) +
  geom_vline(xintercept=0, linetype="dashed") +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
    panel.background = element_rect(colour = "black", linewidth=1),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    plot.margin = margin(1, 1, 1, 1, "pt"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("Promoter") +
  xlab(expression(beta))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/top10_diff.svg", 
       height = 100, width = 60, units = "mm")

# vector of model suffixes
models <- c("eq1", "eq2.m", "eq3_age")

# create a named list of results
top10_list <- lapply(models, function(model) {
  
  beta_col <- paste0("beta_", model)
  fdr_col  <- paste0("fdr_", model)
  
  df %>%
    filter(.data[[fdr_col]] < 0.05) %>%
    drop_na() %>%
    arrange(desc(.data[[beta_col]])) %>%
    dplyr::slice(c(1:10, (n() - 9):n())) %>%
    mutate(mod = model)
})

#Name the list elements
names(top10_list)<- models
top_10<- do.call(rbind, top10_list)

mod_names<- as_labeller(c("eq1" = "Eq.1", "eq3_age" = "Eq.3 Within"))

top_10 %>%
  filter(mod == "eq1" | mod == "eq3_age") %>%
  dplyr::select(beta_eq1, beta_eq3_age, se_beta_eq1, se_beta_eq3_age, gene_name, mod) %>%
  dplyr::rename(se_eq1 = se_beta_eq1, se_eq3_age = se_beta_eq3_age) %>%
  pivot_longer(cols = c(beta_eq1, beta_eq3_age, se_eq1, se_eq3_age),
               names_to = c(".value", "model"),
               names_sep = "_") %>%
ggplot(aes(x=beta, y=reorder(gene_name, beta), colour = model)) +
  geom_point(aes(size = se)) +
  scale_size(range = c(0.05, 2)) +
  geom_path(aes(group = gene_name), colour = "black") +
  scale_colour_manual(values = c("steelblue2", "purple")) +
  geom_vline(xintercept=0, linetype="dashed") +
  theme_classic(base_size = 6) +
  theme(
    panel.background = element_rect(colour = "black", linewidth=1),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    plot.margin = margin(1, 1, 1, 1, "pt"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  facet_wrap(vars(mod), ncol =2, scales = "free_y", labeller = mod_names) +
  ylab("Promoter") +
  xlab(expression(beta))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/top10_proms_dnam.svg", 
       height = 75, width = 105, units = "mm")

#GSEA
#Generate hallmark gene set
hallmark.msigdb = msigdbr(species = "Macaca mulatta", category = "H")
hallmark_list = split(x = hallmark.msigdb$ensembl_gene, f = hallmark.msigdb$gs_name)

proms_gsea<- df %>% 
  dplyr::select(c(outcome, beta_eq1)) %>% 
  arrange(desc(beta_eq1))

proms_gsea2<- proms_gsea$beta_eq1
names(proms_gsea2) = proms_gsea$outcome

#Enrichment for Hallmark set
eq1_gsea_out<- fgsea(pathways = hallmark_list, 
                   stats = proms_gsea2,
                   minSize = 15,
                   maxSize = 500,
                   eps = 0.0)

eq1_gsea_out %>%
  arrange(NES) %>%
  dplyr::slice(c(1:10, (n() - 9):n())) %>%
  ggplot(aes(x=NES, y=reorder(pathway, NES))) +
  #geom_col(aes(alpha = padj<.05)) +
  geom_point(colour = 'steelblue2') +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #scale_colour_manual(values = c("steelblue2")) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("Pathway") +
  xlab("NES")

proms_gsea<- df %>% 
  dplyr::select(c(outcome, beta_eq3_age)) %>% 
  arrange(desc(beta_eqe_age))

proms_gsea2<- proms_gsea$beta_eq3_age
names(proms_gsea2) = proms_gsea$outcome

#Enrichment for Hallmark set
eq3_gsea_out<- fgsea(pathways = hallmark_list, 
                     stats = proms_gsea2,
                     minSize = 15,
                     maxSize = 500,
                     eps = 0.0)

eq3_gsea_out %>%
  arrange(NES) %>%
  dplyr::slice(c(1:10, (n() - 9):n())) %>%
  ggplot(aes(x=NES, y=reorder(pathway, NES))) +
  #geom_col(aes(alpha = padj<.05)) +
  geom_point(colour = "purple") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #scale_colour_manual(values = c("steelblue2", "purple")) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("Pathway") +
  xlab("NES")

proms_gsea<- df %>% 
  dplyr::select(c(outcome, eq3_chron_diff)) %>% 
  arrange(desc(eq3_chron_diff))

proms_gsea2<- proms_gsea$eq3_chron_diff
names(proms_gsea2) = proms_gsea$outcome

#Enrichment for Hallmark set
diff_gsea_out<- fgsea(pathways = hallmark_list, 
                     stats = proms_gsea2,
                     minSize = 15,
                     maxSize = 500,
                     eps = 0.0)

diff_gsea_out$pathway<- gsub("HALLMARK_", "", diff_gsea_out$pathway)

diff_gsea_out %>%
  arrange(NES) %>%
  dplyr::slice(c(1:10, (n() - 9):n())) %>%
  ggplot(aes(x=NES, y=reorder(pathway, NES), colour = NES < 0)) +
  #geom_col(aes(alpha = padj<.05)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c("steelblue2", "purple")) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("Pathway") +
  xlab("NES")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/diff_gsea.svg", 
       height = 100, width = 60, units = "mm")

#DNA vs RNA---------------------------------------------------------------------
dna_rna<- inner_join(rna_int, prom_df, by = 'outcome')

dna_rna %>%
  mutate(signif = ifelse(pval_eq3_age < .20 & dnam_fdr_eq3_age < .20, "Y", "N")) %>%
  #filter(pval_eq3_age < .05 & dnam_fdr_eq3_age < .05) %>%
  ggplot(aes(dnam_beta_eq3_age, beta_eq3_age)) +
  geom_point(aes(alpha = 0.3, colour = signif), size = 0.1) +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c("purple", "purple4")) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  scale_y_continuous(breaks = seq(-1.0, 1.0, 0.5), limits = c(-1.0, 1.0)) +
  scale_x_continuous(breaks = seq(-0.2, 0.2, 0.1), limits = c(-0.2, 0.2)) +
  xlab(expression(beta["DNAm"])) +
  ylab(expression(beta["GE"]))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/dnam_rna_eq3.svg", 
       height = 50, width = 50, units = "mm")

dna_rna %>%
  mutate(signif = ifelse(pval_chron_age < .20 & dnam_fdr_eq1 < .20, "Y", "N")) %>%
  #filter(pval_chron_age < .05 & dnam_fdr_eq1 < .05) %>%
  ggplot(aes(dnam_beta_eq1, beta_chron_age)) +
  geom_point(aes(alpha = 0.3, colour = signif), size = 0.1) +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c("steelblue1", "steelblue4")) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1,
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  scale_y_continuous(breaks = seq(-0.10, 0.10, 0.05), limits = c(-0.10, 0.10)) +
  scale_x_continuous(breaks = seq(-0.10, 0.10, 0.05), limits = c(-0.10, 0.10)) +
  xlab(expression(beta["DNAm"])) +
  ylab(expression(beta["GE"]))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/dnam_rna_eq1.svg", 
       height = 50, width = 50, units = "mm")

save.image("/home/ckelsey4/rna_data/rna_analysis.RData")
