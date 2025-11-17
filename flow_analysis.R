library(tidyverse)
library(flowCore)

#List all .fcs files from directory, the loop iterates through all directory hierarchies
fcs_files<- list.files("/data/CEM/smacklab/cayo_flow_cytometry", pattern = ".fcs", recursive = T, full.names = T)
fcs_files2<- fcs_files[1:5]

#Read in all fcs files
fcs_set<- read.flowSet(fcs_files2, emptyValue = F)

#Read in compensation matrices

