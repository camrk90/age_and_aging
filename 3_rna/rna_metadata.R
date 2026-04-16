library(tidyverse)
library(ggplot2)

base_meta<- read.table("/home/ckelsey4/rna_data/base_meta.txt")

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
#m_counts<- clean_mat("/home/ckelsey4/Cayo_meth/rna_seq/Mitchell_counts_dedup_all")
#r_counts<- clean_mat("/home/ckelsey4/Cayo_meth/rna_seq/Rachel_counts_dedup_all")

rna_counts<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_counts_9Jan26.rds")
rna_samples<- colnames(rna_counts)

rm(rna_counts)

#Import metadata
meta_vanderbilt<- readRDS("/home/ckelsey4/Cayo_meth/rna_seq/Cayo_PBMC_longLPS_metadata_wELA_18Jul25.rds")
meta_asu<- readRDS("/home/ckelsey4/data/metadata_Ashlee.rds")

####Import sequencing stats and match to count matrix samples
seq_stat1<- read_csv("/home/ckelsey4/Cayo_meth/rna_seq/9Jan26_Mitchell_seq_stats.csv")
seq_stat2<- read_csv("/home/ckelsey4/Cayo_meth/rna_seq/9Jan26_Rachel_seq_stats.csv")

seq_stat2<- seq_stat2 %>%
  dplyr::select(-dup_Sample_ID)

seq_stats<- rbind(seq_stat1, seq_stat2)

#Remove duplicate samples by minimum unique_mapped values
seq_stats<- seq_stats %>%
  filter(!c(Sample_ID == "RID_0172" & uniq_mapped == 1064558)) %>%
  filter(!c(Sample_ID == "RID_0517" & uniq_mapped == 732267))

rm(seq_stat1);rm(seq_stat2)

#Subset vanderbilt metadata to relevant variables
rna_meta<- meta_vanderbilt %>%
  dplyr::select(c(animal_ID, sex, trapped_age, trapping_ID, Sample_ID, Seq_batch, Stimulation))

#Subset ASU metadata to relevant variables and capitalize sex/stimulation vars
# to match vanderbilt df
rna_meta2<- meta_asu %>%
  dplyr::select(c(cayoid, sex, age, trapid, lid, batch, stimulation)) %>%
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

rna_meta<- rna_meta[rna_meta$Sample_ID %in% rna_samples,]

rna_meta<- rna_meta %>%
  mutate(across(9:17, as.numeric))

base_meta<- rna_meta %>%
  filter(Stimulation == "H2O")

base_meta<- base_meta %>% 
  group_by(animal_ID) %>% 
  mutate(n = n()) %>%
  filter(n > 1) %>%
  ungroup()

write.table(base_meta, "/home/ckelsey4/rna_data/base_meta.txt", quote = F)
write.table(rna_meta, "/home/ckelsey4/rna_data/rna_meta.txt", quote = F)

#Plot Age and Sample Distributions----------------------------------------------
base_meta<- base_meta %>%
  group_by(animal_ID) %>%
  mutate(min_age = min(trapped_age))
base_meta$trapped_age<- round(base_meta$trapped_age, 0)

samples_dist_rna<- base_meta %>%
  ggplot(aes(x=trapped_age, y=reorder(animal_ID, min_age), colour=sex)) +
  geom_path(linewidth = 0.5) +
  geom_point(colour="black", size = 0.25) +
  scale_x_continuous(breaks = seq(5, 30, by=5)) +
  coord_cartesian(xlim = c(5, 30)) +
  scale_colour_manual(values = c("red3", "pink2"), name = "Sex") +
  ylab("Individual") +
  xlab("Age") +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin = margin(1, 1, 1, 1, "pt")) +
  xlab("Age") +
  ylab("Individual")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/samples_dist_rna.svg", 
       samples_dist_rna, 
       height = 85, width = 55, units = "mm")

#N Samples
n_samples_rna<- base_meta %>%
  ggplot(aes(x=n, fill=as.factor(sex))) +
  geom_bar(position = 'dodge') +
  scale_fill_manual(values = c("red3", "pink2"), name = "Sex") +
  ylab("Count") +
  xlab("N Samples") +
  theme_classic(base_size=6) +
  theme(panel.background = element_rect(colour = "black", linewidth=0.5),
        axis.line = element_line(colour = "black", linewidth = 0.25),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        legend.position = "none")

ggsave("/home/ckelsey4/Cayo_meth/aging_plots/n_samples_rna.svg", 
       n_samples_rna, 
       height = 20, width = 30, units = "mm")



