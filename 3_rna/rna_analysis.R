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

#Use for local
parent_dir<- paste0(getwd(), "/local_data/")

#Use for remote (SOL)
parent_dir<- "/home/ckelsey4/rna_data/"

load(paste0(parent_dir, "rna_analysis.RData"))

#Load model data
eq1_int<- readRDS(paste0(parent_dir, "rna_eq1"))
eq2_int<- readRDS(paste0(parent_dir, "rna_eq2"))
eq3_int<- readRDS(paste0(parent_dir, "rna_eq3"))

#Load metadata
base_meta<- read.table(paste0(parent_dir, "base_meta.txt"))
rna_counts<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_counts_9Jan26.rds")

base_meta<- base_meta %>%
  arrange(Sample_ID) %>%
  filter(Seq_batch %in% c(1, 2, 3)) %>%
  mutate(y = 1)

base_meta<- base_meta[base_meta$Sample_ID %in% colnames(rna_counts),]
rna_counts<- rna_counts[, base_meta$Sample_ID]

#Make simplified outcome df
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

rna_int<- rna_int %>%
  mutate(eq3_chron_diff = beta_eq3_age - beta_chron_age,
         eq2b_chron_diff = beta_eq2_m - beta_chron_age,
         eq3_chron_ratio = beta_eq3_age/beta_chron_age,
         eq2b_chron_ratio = beta_eq2_m/beta_chron_age)

#PCA----------------------------------------------------------------------------
#Generate normalized counts
if (all.equal(base_meta$Sample_ID, colnames(rna_counts)) == T) {
  
  rna_norm<- voom(calcNormFactors(DGEList(counts=rna_counts)), plot=FALSE)
  rna_norm<- rna_norm[["E"]]
  colnames(rna_norm)<- colnames(rna_counts)
  rownames(rna_norm)<- rownames(rna_counts)
  rna_norm<- t(rna_norm)
  
} else {
  print("Metadata Sample IDs and RNA Count cols do not match")
}

rna_norm<- rna_norm[rownames(rna_norm) %in% base_meta$Sample_ID,]

rna_pca<- prcomp(rna_norm, center = TRUE, scale. = TRUE)

pcs<- as.data.frame(rna_pca$x)

pcs<- cbind(pcs[1:5], base_meta)

pc.matrix<- model.matrix(~ PC1 + PC2 + PC3 + PC4 + PC5 + trapped_age + within_age + mean_age + sex + Seq_batch + 
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
    
    dfr<- cor.test(y, rna_int$beta_chron_age)$estimate
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
  colnames(signif)<- "genes_shared_with_eq1"
  
  signif$predictor<- gsub("pval_", "", rownames(signif))
  
  #Combine outputs
  counts<- left_join(counts, test, by = "predictor")
  counts<- left_join(counts, signif, by = "predictor")
  
  #Plot counts
  counts_plot<- counts %>%
    filter(!predictor %in% c("eq3_m")) %>%
    ggplot(aes(reorder(predictor, count), count, fill = predictor)) +
    geom_bar(stat = 'identity') +
    geom_text(label=counts$count[!counts$predictor %in% c("eq3_m")], vjust=-0.25, size = 2) +
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

counts[["df"]]

counts[["plot"]] +
  scale_fill_manual(values = c("steelblue2", "grey30", "green4", "purple")) +
  scale_x_discrete(labels = c("chron_age" = "Eq.1", "eq2_w" = "Eq.2 W.",
                              "eq2_m" = "Eq.2 B.", "eq3_age" = "Eq.3 W."))

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/signif_counts.svg",
       height = 50, width = 50, units = "mm")

### Effect sizes----------------------------------------------------------------
# Plot distribution of effect sizes for models
rna_int %>%
  dplyr::select(c(beta_eq3_age, beta_eq2_w, beta_chron_age, beta_eq2_m)) %>%
  pivot_longer(cols = c(beta_eq3_age, beta_eq2_w, beta_chron_age, beta_eq2_m),
               values_to = 'beta',
               names_to = 'var') %>%
  mutate(var = factor(var, levels = c("beta_eq2_w", "beta_eq3_age",
                                      "beta_eq2_m", "beta_chron_age"))) %>% 
  #filter(var %in% c("beta_chron_age", "beta_eq2_m", "beta_eq3_age")) %>%
  ggplot(aes(x=var, y=beta, fill=var)) +
  geom_violin() +
  geom_boxplot(width = 0.05, fill = "white", outlier.size = 0.25) +
  #geom_density(alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 'dashed', colour = "red") +
  scale_fill_manual(values = c("green4", "purple", 'grey30', 'steelblue2')) +
  theme_classic(base_size=6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt")) +
  scale_x_discrete(labels = c('Eq.2 W.', 'Eq.3 W', "Eq.2 B.", "Eq.1")) +
  ylim(-0.5, 0.5) +
  ylab(expression(beta)) +
  xlab("Var")

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
      ggplot(aes(x={{var1}}, y={{var2}}, colour = diff)) +
      geom_point(size = 1, alpha = 0.8) +
      geom_abline() +
      geom_smooth(method = "lm", linewidth = 0.5) +
      geom_vline(xintercept=0, linetype="dashed") +
      geom_hline(yintercept=0, linetype="dashed") +
      theme_classic(base_size = 6) +
      theme(legend.position = "top",
            legend.key.width = unit(5, 'mm'), 
            legend.key.height = unit(2, 'mm')) +
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
  ylab(expression(beta["Eq.3"])) +
  xlim(-0.5, 0.5) +
  ylim(-0.5, 0.5)

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq2_eq3_scatter_rna.svg",
       height = 50, width = 50, units = "mm")

# Eq.1 Age vs Eq.3 Age
## Scatterplot
compare_plot(rna_int, pval_chron_age, pval_eq3_age, 
             beta_chron_age, beta_eq3_age, "scatter") +
  scale_color_gradient2(low = "steelblue2", mid = "grey70", high = "purple", 
                        midpoint = 0, name = "") +
  xlab(expression(beta["Eq.1"])) +
  ylab(expression(beta["Eq.3"]))  +
  xlim(-0.5, 0.5) +
  ylim(-0.5, 0.5)

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/within_chron_scatter_rna.svg",
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
             beta_eq2_w, beta_eq2_m,"scatter") +
  scale_colour_gradient2(low = "green4", mid = "grey70", high = "grey30", 
                       midpoint = 0, name = "") +
  xlab(expression(beta["Eq.2 W"])) +
  ylab(expression(beta["Eq.2 B"])) +
  xlim(-0.5, 0.5) +
  ylim(-0.5, 0.5)

# Eq.1 Age vs Eq.2 Between Age
compare_plot(rna_int, pval_chron_age, pval_eq2_m, 
             beta_chron_age, beta_eq2_m, "scatter") +
  scale_colour_gradient2(low = "steelblue2", mid = "grey70", high = "grey30", 
                         midpoint = 0, name = "") +
  xlab(expression(beta["Eq.1"])) +
  ylab(expression(beta["Eq.2 B."]))  +
  xlim(-0.2, 0.2) +
  ylim(-0.2, 0.2)

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/between_chron_scatter_rna.svg", 
       height = 50, width = 50, units = "mm")

## Histogram
compare_plot(rna_int, pval_chron_age, pval_eq2_m, 
             beta_chron_age, beta_eq2_m, "hist") +
  scale_colour_gradient2(low = "steelblue2", mid = "grey70", high = "grey30", 
                         midpoint = 0, name = "") +
  xlab(expression(beta["Eq.2 B."] - beta["Eq.1"])) +
  ylab("Count")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/within_chron_hist_rna.svg", 
       height = 50, width = 50, units = "mm")

### Upset Plot------------------------------------------------------------------
## Significant regions Venn diagram
chron.age<- rna_int$outcome[rna_int$pval_chron_age < 0.05]
age.w<- rna_int$outcome[rna_int$pval_eq2_w < 0.05]
eq2.btwn<- rna_int$outcome[rna_int$pval_eq2_m < 0.05]
eq3<- rna_int$outcome[rna_int$pval_eq3_age < 0.05]

venn_all<- list(chron.age, age.w, eq3, eq2.btwn)
names(venn_all)<- c("Eq.1", "Eq.2 Within", "Eq.3 Within", "Eq.2 Btwn")

upset(fromList(venn_all), order.by = "freq", 
      text.scale = c(1, 1, 1, 1, 1, 1), 
      line.size = 1, point.size = 2)

#Plot top genes-----------------------------------------------------------------
#Collect all macaque genes
mm_mart<- useEnsembl(biomart="genes", dataset="mmulatta_gene_ensembl")
mm_genes<- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                 mart = mm_mart)
colnames(mm_genes)<- c("anno", "gene_name")

#Replace ENSMMUG names with gene names 
rna_genes<- rna_int$outcome
rna_genes2<- rna_genes[grepl("ENSMMUG*", rna_genes)]
mm_genes2<- mm_genes[mm_genes$anno %in% rna_genes2, ]

mm_genes2 <- mm_genes2 %>%
  mutate(gene_name = ifelse(gene_name == "", anno, gene_name)) %>%
  arrange(anno)

rna_int$outcome[rna_int$outcome %in% mm_genes2$anno]<- mm_genes2$gene_name

rm(rna_genes);rm(rna_genes2)

eq1_top20<- rna_int %>%
  filter(!grepl("ENSMMUG", outcome)) %>%
  filter(pval_chron_age < .05) %>%
  arrange(desc(beta_chron_age)) %>%
  dplyr::slice(c(1:10, (n() - 9):n()))

eq3_top20<- rna_int %>%
  filter(!grepl("ENSMMUG", outcome)) %>%
  filter(pval_eq3_age < .05) %>%
  arrange(desc(beta_eq3_age)) %>%
  dplyr::slice(c(1:10, (n() - 9):n()))

eq2b_top20<- rna_int %>%
  filter(!grepl("ENSMMUG", outcome)) %>%
  filter(pval_eq2_m < .05) %>%
  arrange(desc(beta_eq2_m)) %>%
  dplyr::slice(c(1:10, (n() - 9):n()))

options(scipen = 1000)

rna_int %>%
  ggplot(aes(beta_chron_age, -log10(pval_chron_age), colour = -log10(pval_chron_age) < -log10(0.05))) +
  geom_point(alpha = 0.5, size = 0.05) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_colour_manual(values = c('steelblue4', 'steelblue1')) +
  geom_label_repel(
    data = eq1_top20,
    aes(label = outcome),
    size = 1,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.size = 0.1,
    label.padding = 0.1
  ) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1) +
  xlab(expression(beta["Eq.1"])) +
  scale_x_continuous(breaks = seq(-0.15, 0.15, 0.1), limits = c(-0.15, 0.15)) +
  ylab("-log10(P-value)")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq1_volcano.svg", 
       height = 50, width = 50, units = "mm")

rna_int %>%
  ggplot(aes(beta_eq3_age, -log10(pval_eq3_age), colour = -log10(pval_eq3_age) < -log10(0.05))) +
  geom_point(alpha = 0.5, size = 0.05) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_colour_manual(values = c("purple4","purple")) +
  geom_label_repel(
    data = eq3_top20,
    aes(label = outcome),
    size = 1,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.size = 0.1,
    label.padding = 0.1
  ) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1) +
  scale_x_continuous(breaks = seq(-0.6, 0.6, 0.3), limits = c(-0.6, 0.6)) +
  xlab(expression(beta["Eq.3"])) +
  ylab("-log10(P-value)")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq3_volcano.svg", 
       height = 50, width = 50, units = "mm")

rna_int %>%
  ggplot(aes(beta_eq2_m, -log10(pval_eq2_m), colour = -log10(pval_eq2_m) < -log10(0.05))) +
  geom_point(alpha = 0.5, size = 0.05) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_colour_manual(values = c("grey30","grey80")) +
  geom_label_repel(
    data = eq2b_top20,
    aes(label = outcome),
    size = 1,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.size = 0.1,
    label.padding = 0.1
  ) +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        aspect.ratio = 1) +
  scale_x_continuous(breaks = seq(-0.2, 0.2, 0.1), limits = c(-0.2, 0.2)) +
  xlab(expression(beta["Eq.2. B."])) +
  ylab("-log10(P-value)")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/eq2_m_volcano_rna.svg", 
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
  mutate(pathway = gsub("_", " ", pathway),
         model = "Eq.1")

chron_gsea_df %>%
  dplyr::slice(c(1:10, (n() - 9):n())) %>%
  filter(abs(NES) > 1) %>%
  ggplot(aes(x=reorder(pathway, NES), y=NES)) +
  #geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  geom_point(aes(size = padj, alpha = padj < .05), colour = 'steelblue2') +
  theme_classic(base_size = 6) +
  theme(
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("NES") +
  xlab("Hallmark Gene Set") +
  coord_flip()

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/gsea_chron.svg", 
       height = 80, width = 85, units = "mm")

eq3_gsea<- run_gsea(outcome, beta_eq3_age, rna_int, hallmark_list)
eq3_gsea<- eq3_gsea %>%
  mutate(pathway = gsub("_", " ", pathway),
         model = "Eq.3")

eq3_gsea %>%
  dplyr::slice(c(1:10, (n() - 9):n())) %>%
  filter(abs(NES) > 1) %>%
  ggplot(aes(x=reorder(pathway, NES), y=NES)) +
  #geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  geom_point(aes(size = padj, alpha = padj < .05), colour = 'purple') +
  theme_classic(base_size = 6) +
  theme(
    panel.background = element_rect(colour = "black", linewidth=1),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    plot.margin = margin(1, 1, 1, 1, "pt"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("NES") +
  xlab("Hallmark Gene Set") +
  coord_flip()

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/gsea_eq3.svg",
       height = 80, width = 85, units = "mm")

topPathwaysUp <- eq3_gsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- eq3_gsea[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(hallmark_list[topPathways], exampleRanks, eq3_gsea, 
              gseaParam=0.5)

full_gsea<- rbind(chron_gsea,eq3_gsea)

full_gsea %>%
  arrange(NES) %>%
  #dplyr::slice(c(1:15, (n() - 14):n())) %>%
  ggplot(aes(x=reorder(pathway, NES), y=NES, colour = model)) +
  #geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  geom_point(aes(size = padj, alpha = padj < .05)) +
  scale_colour_manual(values = c("steelblue2", "purple")) +
  theme_classic(base_size = 6) +
  theme(
    panel.background = element_rect(colour = "black", linewidth=1),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    plot.margin = margin(1, 1, 1, 1, "pt"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("NES") +
  xlab("Hallmark Gene Set") +
  coord_flip()

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/gsea_full.svg",
       height = 100, width = 85, units = "mm")

eq2_m_gsea<- run_gsea(outcome, beta_eq2_m, rna_int, hallmark_list)
eq2_m_gsea<- eq2_m_gsea %>%
  mutate(pathway = gsub("_", " ", pathway),
         model = "Eq.2 B.")

eq2_m_gsea %>%
  dplyr::slice(c(1:10, (n() - 9):n())) %>%
  filter(abs(NES) > 1) %>%
  ggplot(aes(x=reorder(pathway, NES), y=NES)) +
  #geom_col(aes(alpha=padj<0.05), position = position_dodge(0.5), colour="black") +
  geom_point(aes(size = padj, alpha = padj < .05), colour = 'grey30') +
  theme_classic(base_size = 6) +
  theme(
    panel.background = element_rect(colour = "black", linewidth=1),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    plot.margin = margin(1, 1, 1, 1, "pt"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("NES") +
  xlab("Hallmark Gene Set") +
  coord_flip()

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/gsea_eq2b.svg",
       height = 80, width = 85, units = "mm")

full_gsea<- inner_join(chron_gsea[,1:6], eq3_gsea[,1:6], suffix = c("_eq1", "_eq3_w"), by = "pathway")
full_gsea<- inner_join(full_gsea, eq2_m_gsea[,1:6], by = "pathway")
colnames(full_gsea)[12:16] <- c(paste(colnames(full_gsea)[12:16], "_eq2_m", sep = ""))

full_gsea %>%
  ggplot(aes(NES_eq3_w, NES_eq1, shape = padj_eq1 < .05, colour = padj_eq3_w < .05)) +
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

# Compare estimates to public data----------------------------------------------
### Dataset from doi: 10.1186/s13059-019-1840-y
hg_compare<- readxl::read_xls(paste0(parent_dir, "13059_2019_1840_MOESM2_ESM.xls"))
hg_compare<- hg_compare %>%
  dplyr::rename(outcome = HGNCID)

hg_compare<- hg_compare %>%
  filter(outcome %in% rna_int$outcome)

rna_compare<- rna_int %>%
  filter(outcome %in% hg_compare$outcome)

plot(hg_compare$`Age effect size`, rna_compare$beta_eq3_age)

### Split hg df into deciles where list element 1 is the lowest-ranked age effects
### and list element 10 is the highest-ranked
hg_compare<- hg_compare %>%
  mutate(age_decile = ntile(dplyr::desc(`Age effect size`), 10))

hg_compare_list<- hg_compare %>%
  dplyr::select(outcome, age_decile)

hg_compare_list<- split(hg_compare_list$outcome, hg_compare_list$age_decile)

names(hg_compare_list)<- paste0("decile", 1:10)

rna_compare_vector<- rna_compare %>% 
  dplyr::select(outcome, beta_chron_age) %>% 
  arrange(desc(beta_chron_age))

rna_compare_vector2<- rna_compare_vector$beta_chron_age
names(rna_compare_vector2) = rna_compare_vector$outcome

options(scipen = 100000)
hg_compare_gsea<- fgsea(pathways = hg_compare_list, 
                      stats = rna_compare_vector2,
                      eps = 0.0)
hg_compare_gsea$pathway<- factor(hg_compare_gsea$pathway, levels = rev(str_sort(hg_compare_gsea$pathway, numeric = T)))

hg_compare_gsea %>%
  arrange(NES) %>%
  ggplot(aes(x=NES, y=pathway, colour = NES < 0)) +
  geom_point(aes(alpha = padj < .05)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #scale_colour_manual(values = c("steelblue2", "purple")) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("Pathway") +
  xlab("NES") +
  ggtitle("Eq.1")

### Comparison with Marina's paper
mm_compare<- readxl::read_xlsx(paste0(parent_dir, "pnas.2121663119.sd02.xlsx"))
mm_compare<- mm_compare %>%
  dplyr::rename(outcome = `common gene name`)

mm_compare<- mm_compare %>%
  filter(`FDR aging` < .05) %>%
  filter(outcome %in% rna_int$outcome)

rna_compare2<- rna_int %>%
  filter(outcome %in% mm_compare$outcome)

### Split hg df into deciles where list element 10 is the lowest-ranked age effects
### and list element 1 is the highest-ranked
mm_compare<- mm_compare %>%
  mutate(age_decile = ntile(dplyr::desc(`beta aging`), 10),
         direction = ifelse(`beta aging` < 0, "downreg", "upreg"))

rna_compare2<- rna_compare2 %>%
  mutate(age_decile = ntile(dplyr::desc(beta_eq3_age), 10))

mm_compare_list<- mm_compare %>%
  dplyr::select(outcome, direction)

mm_compare_list<- split(mm_compare_list$outcome, mm_compare_list$direction)

names(mm_compare_list)<- paste0("decile", 1:10)

rna_compare_vector<- rna_compare2 %>% 
  dplyr::select(outcome, beta_eq3_age) %>% 
  arrange(desc(beta_eq3_age))

rna_compare_vector2<- rna_compare_vector$beta_eq3_age
names(rna_compare_vector2) = rna_compare_vector$outcome

options(scipen = 100000)

mm_compare_gsea<- fgsea(pathways = mm_compare_list, 
                     stats = rna_compare_vector2,
                     eps = 0.0)

mm_compare_gsea$pathway<- factor(mm_compare_gsea$pathway, levels = rev(str_sort(mm_compare_gsea$pathway, numeric = T)))

mm_compare_gsea %>%
  arrange(NES) %>%
  ggplot(aes(x=NES, y=pathway, colour = NES < 0)) +
  geom_point(aes(alpha = padj < .05)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic(base_size = 12) +
  theme(
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("Pathway") +
  xlab("NES") +
  ggtitle("Marina - Eq.3")

### Comparison with Marquez paper
marquez_compare<- readxl::read_xlsx(paste0(parent_dir, "marquez_aging.xlsx"), 
                                    sheet = "Figure 3a 3b")
marquez_compare<- marquez_compare %>%
  dplyr::rename(outcome = GeneName)

marquez_compare<- marquez_compare %>%
  filter(outcome %in% rna_int$outcome)

rna_compare3<- rna_int %>%
  filter(outcome %in% marquez_compare$outcome)

### Split hg df into deciles where list element 10 is the lowest-ranked age effects
### and list element 1 is the highest-ranked
marquez_compare<- marquez_compare %>%
  mutate(age_decile = ntile(dplyr::desc(Males.rna), 10))

marquez_compare_list<- marquez_compare %>%
  dplyr::select(outcome, age_decile)

marquez_compare_list<- split(marquez_compare_list$outcome, marquez_compare_list$age_decile)

names(marquez_compare_list)<- paste0("decile", 1:10)

rna_compare_vector<- rna_compare3 %>% 
  dplyr::select(outcome, beta_chron_age) %>% 
  arrange(desc(beta_chron_age))

rna_compare_vector3<- rna_compare_vector$beta_chron_age
names(rna_compare_vector3) = rna_compare_vector$outcome

options(scipen = 100000)

marquez_compare_gsea<- fgsea(pathways = marquez_compare_list, 
                        stats = rna_compare_vector3,
                        eps = 0.0)

marquez_compare_gsea$pathway<- factor(marquez_compare_gsea$pathway, levels = rev(str_sort(marquez_compare_gsea$pathway, numeric = T)))

marquez_compare_gsea %>%
  arrange(NES) %>%
  ggplot(aes(x=NES, y=pathway, colour = NES < 0)) +
  geom_point(aes(alpha = padj < .05)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey98", linewidth = 0.5)) +
  ylab("Pathway") +
  xlab("NES") +
  ggtitle("Marquez Males - Eq.3")


save.image("/home/ckelsey4/rna_data/rna_analysis.RData")
