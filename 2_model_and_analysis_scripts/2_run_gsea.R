#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mem=20G 

library(tidyverse)
library(fgsea)

#Load data
pqlseq_anno<- readRDS("/scratch/ckelsey4/Cayo_meth/pqlseq_anno.rds")

#Generate region list
func_region<- pqlseq_anno %>%
  select(unique_cpg, anno)

func_region_list<- split(x=func_region$unique_cpg, f=func_region$anno)

#Run gsea
within<- pqlseq_anno %>%
  select(unique_cpg, beta_eq2_w_age) %>%
  arrange(desc(beta_eq2_w_age))

within2<- within$beta_eq2_w_age
names(within2) = within$unique_cpg

rm(pqlseq_anno)

#Enrichment for Hallmark set
within_gsea<- fgseaSimple(pathways = func_region_list, 
                           stats = within2,
                           minSize = 5,
                           maxSize = length(within2) - 1,
                           nperm = 1000)

saveRDS(within_gsea, "/scratch/ckelsey4/Cayo_meth/within_gsea.rds")

