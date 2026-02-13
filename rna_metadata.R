library(tidyverse)
library(EMMREML)
library(lme4)
library(ggplot2)
library(variancePartition)
library(limma)
library(edgeR)
library(ggcorrplot)
library(fgsea)
library(msigdbr)
library(Biostrings)
library(GenomicFeatures)
library(GenomicRanges)

load("/home/ckelsey4/Cayo_meth/rna_compare.RData")

rna_counts<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_counts_9Jan26.rds")
rna_kin<- readRDS("/scratch/ckelsey4/Cayo_meth/rna_kin_matrix.rds")

#### Import full count matrices ####
#Function to import and clean
clean_mat<- function(x) {
  
  counts<- read_table(x, col_names = F)
  colnames(counts)<- counts[2,]
  counts<- counts[-c(1:2), -c((length(counts)-10):length(counts))]
  samples<- colnames(counts[,7:length(counts)])
  samples<- str_split_i(samples, "/", 7)
  samples<- str_split_i(samples, "\\.", 1)
  samples<- c("gene_id", "chr", "start", "end", "strand", "length", samples)
  colnames(counts)<- samples
  return(counts)
  
}

#Import and clean both count datasets
m_counts<- clean_mat("/home/ckelsey4/Cayo_meth/rna_seq/Mitchell_counts_dedup_all")
r_counts<- clean_mat("/home/ckelsey4/Cayo_meth/rna_seq/Rachel_counts_dedup_all")

#Import metadata
meta_vanderbilt<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_metadata_wELA_18Jul25.rds")
meta_asu<- readRDS("/home/ckelsey4/data/metadata_Ashlee.rds")

####Import sequencing stats and match to count matrix samples
seq_stat1<- read_csv("/home/ckelsey4/Cayo_meth/rna_seq/9Jan26_Mitchell_seq_stats.csv")
seq_stat2<- read_csv("/home/ckelsey4/Cayo_meth/rna_seq/9Jan26_Rachel_seq_stats.csv")

seq_stat2<- seq_stat2 %>%
  select(-dup_Sample_ID)

seq_stats<- rbind(seq_stat1, seq_stat2)

#Remove duplicate samples by minimum unique_mapped values
seq_stats<- seq_stats %>%
  filter(!c(Sample_ID == "RID_0172" & uniq_mapped == 1064558)) %>%
  filter(!c(Sample_ID == "RID_0517" & uniq_mapped == 732267))

rm(seq_stat1);rm(seq_stat2)

#Subset vanderbilt metadata to relevant variables
rna_meta<- meta_vanderbilt %>%
  select(c(animal_ID, sex, trapped_age, trapping_ID, Sample_ID, Seq_batch, Stimulation))

#Subset ASU metadata to relevant variables and capitalize sex/stimulation vars
# to match vanderbilt df
rna_meta2<- meta_asu %>%
  select(c(cayoid, sex, age, trapid, lid, batch, stimulation)) %>%
  mutate(sex = toupper(sex),
         stimulation = toupper(stimulation))

rna_meta2$stimulation[rna_meta2$stimulation == "CONT"]<- "H2O"
  
colnames(rna_meta2)<- c("animal_ID", "sex", "trapped_age", "trapping_ID", "Sample_ID", "Seq_batch", "Stimulation")

rna_meta<- rbind(rna_meta, rna_meta2)

rm(rna_meta2)

#Generate within and between age
rna_meta<- rna_meta %>%
  group_by(animal_ID) %>%
  mutate(mean_age = mean(trapped_age)) %>%
  mutate(within_age = trapped_age - mean_age) %>%
  arrange(animal_ID, trapped_age) %>%
  relocate(within_age, .after = trapped_age) %>%
  relocate(mean_age, .after = within_age)

#Make batch factor
rna_meta$Seq_batch<- gsub("batch", "", rna_meta$Seq_batch)
rna_meta$Seq_batch<- factor(rna_meta$Seq_batch, levels = c("1", "2", "3"))

#Join sequencing stats and metadata dfs
rna_meta<- inner_join(rna_meta, seq_stats, by = "Sample_ID")

rna_meta<- rna_meta[rna_meta$Sample_ID %in% colnames(rna_counts),]

rna_meta<- rna_meta %>%
  mutate(across(10:18, as.numeric))

base_meta<- rna_meta %>%
  filter(Stimulation == "H2O")

#Plot rna samples 
base_meta<- base_meta %>%
  group_by(animal_ID) %>%
  mutate(min_age = min(trapped_age))
base_meta$trapped_age<- round(base_meta$trapped_age, 0)

base_meta %>%
  ggplot(aes(x=trapped_age, y=reorder(animal_ID, min_age), colour=sex)) +
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
  #theme(legend.position = "none") +
  theme(panel.background = element_rect(colour = "black", linewidth=2))

#Age distribution
base_meta %>%
  ggplot(aes(x=trapped_age, fill=as.factor(sex))) +
  geom_bar(colour='black', position = 'dodge') +
  scale_x_continuous(breaks = seq(5, 30, by=5)) +
  coord_cartesian(xlim = c(5, 30)) +
  scale_fill_manual(values = c("red3", "pink2"), name = "Sex") +
  ylab("Count") +
  xlab("Age") +
  theme_classic(base_size = 24) +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(colour = "black", linewidth=3))

base_meta<- base_meta %>% 
  group_by(animal_ID) %>%
  mutate(n = n())

#N Samples
base_meta %>% 
  ggplot(aes(x=n, fill=as.factor(sex))) +
  geom_bar(colour='black', position = 'dodge') +
  scale_fill_manual(values = c("red3", "pink2"), name = "Sex") +
  ylab("Count") +
  xlab("N Samples") +
  theme_classic(base_size = 24) +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(colour = "black", linewidth=3))

#Variance Partition-------------------------------------------------------------
vp_model<- ~ within_age + mean_age + (1|sex) + (1|Seq_batch)

vp<- fitExtractVarPartModel(t(rna_counts), vp_model, rna_meta)

plotVarPart(vp)
plotPercentBars(vp[1:10, ])


#Normalize RNA Count Data-------------------------------------------------------
base_meta<- base_meta %>%
  arrange(Sample_ID)

rna_counts<- rna_counts[, base_meta$Sample_ID]

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
  ggplot(aes(PC1, PC2, colour=sex)) +
  geom_point()

pc.matrix<- model.matrix(~ PC1 + PC2 + trapped_age + within_age + mean_age + sex + Seq_batch + p_gene_counts, data = pcs)
pc.matrix %>% 
  cor(use="pairwise.complete.obs") %>%
  ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2)

#Run EMMA for within_age--------------------------------------------------------
if (all.equal(base_meta$Sample_ID, colnames(rna_counts)) == T) {
  
  m_matrix <- model.matrix(~ within_age + mean_age + sex + p_gene_counts, base_meta)
  rna_norm<- voom(calcNormFactors(DGEList(counts=rna_counts)), m_matrix, plot=FALSE)
  rna_norm<- rna_norm[["E"]]
  colnames(rna_norm)<- colnames(rna_counts)
  rownames(rna_norm)<- rownames(rna_counts)
  rna_norm_eq2<- t(rna_norm)
  
} else {
  print("Metadata Sample IDs and RNA Count cols do not match")
}

#Generate model matrix
within_matrix<- model.matrix(~ within_age + mean_age + sex + p_gene_counts, data = base_meta)

re_eq = "y ~ within_age + mean_age + sex + Seq_batch + (1|animal_ID)"

rna_meta$y<- 1

#This generates a random intercept matrix
re_mat <- lFormula(eval(re_eq), base_meta)
re_matZ <- as.matrix(t(re_mat$reTrms$Zt))

rna_kin<- rna_kin[colnames(re_matZ), colnames(re_matZ)]

df_eq2<- data.frame(matrix(nrow=0, ncol=16))

for (i in 1:ncol(rna_norm_eq2)) {
  
  em<- emmreml(rna_norm_eq2[, i], within_matrix, re_matZ, rna_kin, varbetahat=T,varuhat=T, PEVuhat=T, test=T)
  
  em_row<- as.data.frame(cbind(em[["betahat"]][1:4,], em[["Xsqtestbeta"]][1:4,],
                           em[["pvalbeta"]][1:4, "none"]))
  
  colnames(em_row)<- c("beta", "chi_square", "pvalue")
  
  em_row$vars<- rownames(em_row)
  
  em_row<- em_row %>%
    pivot_wider(names_from = vars, values_from = c(beta, chi_square, pvalue))
  
  em_row$Vu<- em[["Vu"]]
  
  em_row$Ve<- em[["Ve"]]
  
  em_row$loglik<- em[["loglik"]]
  
  em_row$outcome<- colnames(rna_norm_eq2)[i]
  
  em_row<- em_row %>% 
    relocate(outcome, .before=`beta_(Intercept)`)
  
  df_eq2<- rbind(df_eq2, em_row)
  
  print(paste("Gene", i, "out of", ncol(rna_norm_eq2), "done!", sep=" "))
  
}

#Run EMMA for Eq. 3-------------------------------------------------------------
if (all.equal(base_meta$Sample_ID, colnames(rna_counts)) == T) {
  
  m_matrix <- model.matrix(~ trapped_age + mean_age + sex + p_gene_counts, base_meta)
  rna_norm<- voom(calcNormFactors(DGEList(counts=rna_counts)), m_matrix, plot=FALSE)
  rna_norm<- rna_norm[["E"]]
  colnames(rna_norm)<- colnames(rna_counts)
  rownames(rna_norm)<- rownames(rna_counts)
  rna_norm_eq3<- t(rna_norm)
  
} else {
  print("Metadata Sample IDs and RNA Count cols do not match")
}

#Generate model matrix
eq3_matrix<- model.matrix(~ trapped_age + mean_age + sex + p_gene_counts, data = base_meta)

re_eq = "y ~ trapped_age + mean_age + sex + p_gene_counts + (1|animal_ID)"

rna_meta$y<- 1

#This generates a random intercept matrix with dims 309 x 114
re_mat <- lFormula(eval(re_eq), base_meta)
re_matZ <- as.matrix(t(re_mat$reTrms$Zt))

df_eq3<- data.frame(matrix(nrow=0, ncol=16))

for (i in 1:ncol(rna_norm_eq3)) {
  
  em<- emmreml(rna_norm_eq3[, i], eq3_matrix, re_matZ, rna_kin, varbetahat=T,varuhat=T, PEVuhat=T, test=T)
  
  em_row<- as.data.frame(cbind(em[["betahat"]][1:4,], em[["Xsqtestbeta"]][1:4,],
                               em[["pvalbeta"]][1:4,"none"]))
  
  colnames(em_row)<- c("beta", "chi_square", "pvalue")
  
  em_row$vars<- rownames(em_row)
  
  em_row<- em_row %>%
    pivot_wider(names_from = vars, values_from = c(beta, chi_square, pvalue))
  
  em_row$Vu<- em[["Vu"]]
  
  em_row$Ve<- em[["Ve"]]
  
  em_row$loglik<- em[["loglik"]]
  
  em_row$outcome<- colnames(rna_norm_eq3)[i]
  
  em_row<- em_row %>% 
    relocate(outcome, .before=`beta_(Intercept)`)
  
  df_eq3<- rbind(df_eq3, em_row)
  
  print(paste("Gene", i, "out of", ncol(rna_norm_eq3), "done!", sep=" "))
  
}

#Run EMMA for chron_age---------------------------------------------------------
if (all.equal(base_meta$Sample_ID, colnames(rna_counts)) == T) {
  
  m_matrix <- model.matrix(~ trapped_age + sex + p_gene_counts, base_meta)
  rna_norm<- voom(calcNormFactors(DGEList(counts=rna_counts)), m_matrix, plot=FALSE)
  rna_norm<- rna_norm[["E"]]
  colnames(rna_norm)<- colnames(rna_counts)
  rownames(rna_norm)<- rownames(rna_counts)
  rna_norm_chron<- t(rna_norm)
  
} else {
  print("Metadata Sample IDs and RNA Count cols do not match")
}

#Generate model matrix
chron_matrix<- model.matrix(~ trapped_age + sex + p_gene_counts, data = base_meta)

re_eq = "y ~ trapped_age + sex + p_gene_counts + (1|animal_ID)"

rna_meta$y<- 1

#This generates a random intercept matrix with dims 309 x 114
re_mat <- lFormula(eval(re_eq), base_meta)
re_matZ <- as.matrix(t(re_mat$reTrms$Zt))

df_chron<- data.frame(matrix(nrow=0, ncol=16))

for (i in 1:ncol(rna_norm_chron)) {
  
  em<- emmreml(rna_norm_chron[, i], chron_matrix, re_matZ, rna_kin, varbetahat=T,varuhat=T, PEVuhat=T, test=T)
  
  em_row<- as.data.frame(cbind(em[["betahat"]][1:4,], em[["Xsqtestbeta"]][1:4,],
                               em[["pvalbeta"]][1:4,"none"]))
  
  colnames(em_row)<- c("beta", "chi_square", "pvalue")
  
  em_row$vars<- rownames(em_row)
  
  em_row<- em_row %>%
    pivot_wider(names_from = vars, values_from = c(beta, chi_square, pvalue))
  
  em_row$Vu<- em[["Vu"]]
  
  em_row$Ve<- em[["Ve"]]
  
  em_row$loglik<- em[["loglik"]]
  
  em_row$outcome<- colnames(rna_norm_chron)[i]
  
  em_row<- em_row %>% 
    relocate(outcome, .before=`beta_(Intercept)`)
  
  df_chron<- rbind(df_chron, em_row)
  
  print(paste("Gene", i, "out of", ncol(rna_norm_chron), "done!", sep=" "))
  
}

#### Plot Model Outcomes ####
#P-Value Distributions
df_chron %>%
  ggplot(aes(pvalue_trapped_age)) +
  geom_histogram(bins=30, colour='black', fill = 'steelblue2') +
  theme_classic(base_size=24) +
  theme(legend.position = "none") +
  xlab("P-Value") +
  ylab("Count")

df_eq2 %>%
  ggplot(aes(pvalue_within_age)) +
  geom_histogram(bins=30, colour='black', fill = 'green4') +
  theme_classic(base_size=24) +
  theme(legend.position = "none") +
  xlab("P-Value") +
  ylab("Count")

df_eq3 %>%
  ggplot(aes(pvalue_trapped_age)) +
  geom_histogram(bins=30, colour='black', fill = 'purple') +
  theme_classic(base_size=24) +
  theme(legend.position = "none") +
  xlab("P-Value") +
  ylab("Count")

df_eq2 %>%
  select(c(pvalue_within_age, pvalue_mean_age)) %>%
  pivot_longer(cols=c(pvalue_within_age, pvalue_mean_age), 
               values_to = "pval",
               names_to = "term") %>%
  ggplot(aes(pval, fill = term)) +
  geom_histogram(bins=30, colour="black", position = position_dodge()) +
  scale_fill_manual(values = c("black", "green4")) +
  theme_classic(base_size=24) +
  theme(legend.position = "none")

df_eq3 %>%
  select(c(pvalue_trapped_age, pvalue_mean_age)) %>%
  pivot_longer(cols=c(pvalue_trapped_age, pvalue_mean_age), 
               values_to = "pval",
               names_to = "term") %>%
  ggplot(aes(pval, fill = term)) +
  geom_histogram(bins=30, colour="black", position = position_dodge()) +
  scale_fill_manual(values = c("black", "purple")) +
  theme_classic(base_size=24) +
  theme(legend.position = "none")

#Effect Sizes
rna_fx<- as.data.frame(cbind(df_chron$outcome, df_chron$beta_trapped_age, df_chron$pvalue_trapped_age,
                             df_eq2$beta_within_age, df_eq2$pvalue_within_age,
                             df_eq2$beta_mean_age, df_eq2$pvalue_mean_age,
                             df_eq3$beta_trapped_age, df_eq3$pvalue_trapped_age,
                             df_eq3$beta_mean_age, df_eq3$pvalue_mean_age))

colnames(rna_fx)<- c("outcome", "beta_chron_age", "pval_chron_age",
                     "beta_eq2_w", "pval_eq2_w",
                     "beta_eq2_m", "pval_eq2_m",
                     "beta_eq3_age", "pval_eq3_age",
                     "beta_eq3_m", "pval_eq3_m")

rna_fx<- rna_fx %>%
  mutate(across(2:11, as.numeric))

compare_plot<- function(df, fdr1, fdr2, var1, var2, c1, c2, lab1, lab2, plot_type) {
  
  df<- df %>%
    mutate(diff = abs({{var2}}) - abs({{var1}})) %>%
    filter({{fdr1}} < .05 | {{fdr2}} < .05)
    

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
      theme_classic(base_size=32) +
      theme(legend.key.width = unit(5, 'cm'), legend.position = "top") +
      theme(panel.background = element_rect(colour = "black", linewidth=3)) +
      #xlim(-0.20, 0.20) +
      #ylim(-0.20, 0.20) +
      xlab(lab1) +
      ylab(lab2)
    
  } else if (plot_type == "hist") {
    
    df %>%
      ggplot(aes(diff, fill = after_stat(x))) +
      geom_histogram(bins = 50, colour="black") +
      geom_vline(xintercept=0, linetype="dashed") +
      scale_fill_gradient2(low = c2, mid = "grey70", high = c1, midpoint = 0, name = "") +
      theme_classic(base_size=32) +
      theme(legend.position = "none") +
      theme(panel.background = element_rect(colour = "black", linewidth=3)) +
      xlab(paste(lab2, "-", lab1, sep=" "))
    
  }
}

compare_plot(rna_fx, pval_chron_age, pval_eq2_w, 
             beta_chron_age, beta_eq2_w,
             "green4", "steelblue2",
             "Chron Age", "Eq2. Within Age", "scatter")

compare_plot(rna_fx, pval_chron_age, pval_eq2_w, 
             beta_chron_age, beta_eq2_w,
             "green4", "steelblue2",
             "Chron Age", "Eq2. Within Age", "hist")

compare_plot(rna_fx, pval_chron_age, pval_eq3_age,
             beta_chron_age, beta_eq3_age,
             "purple", "steelblue2", 
             "Chron Age", "Eq3. Age", "scatter")

compare_plot(rna_fx, pval_chron_age, pval_eq3_age,
             beta_chron_age, beta_eq3_age,
             "purple", "steelblue2", 
             "Chron Age", "Eq3. Age", "hist")

compare_plot(rna_fx, pval_eq2_w, pval_eq3_age,
             beta_eq2_w, beta_eq3_age,
             "purple", "green4",
             "Eq2. Within Age", "Eq3. Age", "scatter")

compare_plot(rna_fx, pval_eq2_w, pval_eq3_age,
             beta_eq2_w, beta_eq3_age,
             "purple", "green4",
             "Eq2. Within Age", "Eq3. Age", "hist")

age.chron.count<- nrow(rna_fx[rna_fx$pval_chron_age < 0.05,])
age.w.count<- nrow(rna_fx[rna_fx$pval_eq2_w < 0.05,])
age.eq3.count<- nrow(rna_fx[rna_fx$pval_eq3_age < 0.05,])
counts<- data.frame(count = c(age.w.count, age.eq3.count, age.chron.count),
                    predictor = as.factor(c('Eq2. Within Age', 'Eq3 Age', 'Chron Age')))

counts<- counts %>%
  mutate(predictor = as.factor(predictor)) %>%
  mutate(predictor=fct_reorder(predictor, count, .desc=T))

rm(age.w.count);rm(age.chron.count);rm(age.eq3.count)

counts %>%
  ggplot(aes(predictor, count, fill = predictor)) +
  geom_bar(stat = 'identity', colour="black") +
  geom_text(label=counts$count, vjust=-1, size=5) +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_text(angle = 15, hjust=0.9),
        legend.position = "none") +
  xlab("Predictor") +
  ylab("Count") +
  scale_fill_manual(values = c("purple", 'green4', 'steelblue2'))

rat<- abs(rna_fx$beta_eq2_w)/abs(rna_fx$beta_chron_age)
rat<- rat[rat < 100]
median(rat)

#### GSEA ####
#make txdb from macaque gff
macaque_txdb<- makeTxDbFromGFF("/scratch/ckelsey4/Cayo_meth/Macaca_mulatta.Mmul_10.110.chr.gtf")

#extract unique chromosomes
unique_chrs <- unique(seqlevels(macaque_txdb))

# Select the first 22(1-20 + X/Y) unique chromosome names as a redundancy against missing chromosome levels in the txdb object
# This will not match to the bsseq object if the chromosomes levels in txdb don't match for whatever reason
selected_chrs <- unique_chrs[1:22]
macaque_txdb <- keepSeqlevels(macaque_txdb, selected_chrs)

#Generate genes and promoters
macaque_genes = genes(macaque_txdb)

macaque_genes$gene_id

mm_genes <- useEnsembl(biomart="genes", dataset="mmulatta_gene_ensembl", version=110)
res <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
             mart = mart)

#Generate hallmark gene set
hallmark.msigdb = msigdbr(species = "Macaca mulatta", category = "H")
hallmark_list = split(x = hallmark.msigdb$ensembl_gene, f = hallmark.msigdb$gs_name)

#Generate gene ontology set
go_set = msigdbr(species = "Macaca mulatta", category = "C5", subcategory = "GO:BP")
go_set = split(x = go_set$ensembl_gene, f = go_set$gs_name)

chron_betas<- rna_fx %>% 
  dplyr::select(outcome, beta_chron_age) %>% 
  arrange(desc(beta_chron_age))

proms2<- proms$diff
names(proms2) = proms$anno


save.image("/home/ckelsey4/Cayo_meth/rna_compare.RData")


em<- emmreml(rna_counts[, 1], predictor_matrix, re_matZ, rna_kin, varbetahat=T,varuhat=T, PEVuhat=T, test=T)

# find number of available cores
ncores <- as.numeric(system("nproc", intern = TRUE))
# create cluster
clus <- makeCluster(ncores)
registerDoParallel(cores = ncores)
# initiate local cluster
clusterExport(clus,
              varlist = c("rna_norm", "rna_meta", "re_matZ", "rna_kin"))


emma_ge_nested_add <- data.frame(parApply(clus, rna_norm, 1, function(y){
  .libPaths(c("~/age_and_aging/",.libPaths()))
  emma <- EMMREML::emmreml(y = y, 
                           X = model.matrix(~ rna_meta$within_age + rna_meta$sex + rna_meta$mean_age + rna_meta$Seq_batch),
                           Z = re_matZ, 
                           K = rna_kin, 
                           varbetahat = T, 
                           varuhat = T, 
                           PEVuhat = T, 
                           test = T)
  return(c(emma$betahat, emma$pvalbeta[, "none"], emma$varbetahat))  
}))



