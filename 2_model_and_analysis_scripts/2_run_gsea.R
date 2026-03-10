#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mem=20G 

library(tidyverse)
library(fgsea)

#Load data
pqlseq_anno<- readRDS("/scratch/ckelsey4/Cayo_meth/pqlseq_anno.rds")

pqlseq_anno<- pqlseq_anno %>%
  mutate(eq2_diff = abs(beta_chron_age) - abs(beta_eq2_w_age),
         eq3_diff = abs(beta_chron_age) - abs(beta_eq3_age))

#Generate region list
func_region<- pqlseq_anno %>%
  select(unique_cpg, anno)

func_region_list<- split(x=func_region$unique_cpg, f=func_region$anno)

#Run gsea
within<- pqlseq_anno %>%
  dplyr::select(unique_cpg, eq2_diff) %>%
  arrange(desc(eq2_diff))

within2<- within$eq2_diff
names(within2) = within$unique_cpg

#Enrichment for Hallmark set
within_gsea<- fgseaSimple(pathways = func_region_list, 
                           stats = within2,
                           minSize = 5,
                           maxSize = length(within2) - 1,
                           nperm = 1000)

saveRDS(within_gsea, "/scratch/ckelsey4/Cayo_meth/within_gsea.rds")

#Run gsea
eq3<- pqlseq_anno %>%
  dplyr::select(unique_cpg, eq3_diff) %>%
  arrange(desc(eq3_diff))

eq3_2<- eq3$eq3_diff
names(eq3_2) = eq3$unique_cpg

#Enrichment for Hallmark set
eq3_gsea<- fgseaSimple(pathways = func_region_list, 
                          stats = eq3_2,
                          minSize = 5,
                          maxSize = length(within2) - 1,
                          nperm = 1000)

saveRDS(eq3_gsea, "/scratch/ckelsey4/Cayo_meth/eq3_gsea.rds")
