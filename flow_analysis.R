library(tidyverse)
library(flowCore)
library(ggcyto)
library(umap)

#List all .fcs files for each year 
fcs_files19<- list.files("/data/CEM/smacklab/cayo_flow_cytometry/Flow_Data_2019", 
                         pattern = ".fcs", recursive = T, full.names = T)
fcs_files21<- list.files("/data/CEM/smacklab/cayo_flow_cytometry/Flow_Data_2021", 
                         pattern = ".fcs", recursive = T, full.names = T)
fcs_files23<- list.files("/data/CEM/smacklab/cayo_flow_cytometry/Flow_Data_2023", 
                         pattern = ".fcs", recursive = T, full.names = T)
fcs_files24<- list.files("/data/CEM/smacklab/cayo_flow_cytometry/Flow_Data_2024", 
                         pattern = ".fcs", recursive = T, full.names = T)

#Read in compensation matrices
comp_files<- list.files("/home/ckelsey4/cayo_flow_data", pattern = ".csv", full.names=T)

comp_list<- lapply(comp_files, function(i){
  mtx<- read.csv(i)
  colnames(mtx)<- gsub("\\.", "-", colnames(mtx))
  rownames(mtx)<- mtx[,1]
  mtx<- mtx[,-1]
})

comp_19<- comp_list[[1]]
comp_21<- comp_list[[2]]
rm(comp_list)

#Read in, compensate, and rename .fcs files by year
fcs_in<- function(fcs_list, comp_mtx){
  
  fcs_set<- read.flowSet(fcs_list, emptyValue = F, transformation = F)
  fcs_comp<- fsApply(fcs_set, function(ff) compensate(ff, comp_mtx))
  rm(fcs_set)
  
  ids<- lapply(sampleNames(fcs_comp), function(i) {fcs_comp@frames[[i]]@description[["$CELLS"]]})
  ids<- do.call(rbind, ids)
  ids<- make.unique(ids, sep = "_")
  
  if(length(ids) != length(sampleNames(fcs_comp))) stop("new_names must match number of frames")
  sampleNames(fcs_comp) <- ids
  # optionally update phenoData column
  pData(fcs_comp)$new_name <- sampleNames(fcs_comp)
  
  return(fcs_comp)
}


#Read in all fcs files
fcs_19<- fcs_in(fcs_files19, comp_19)
fcs_21<- fcs_in(fcs_files21, comp_21)
fcs_23<- fcs_in(fcs_files23, comp_21)
fcs_24<- fcs_in(fcs_files24, comp_21)

autoplot(df, "FSC-A", "SSC-A")

sum(length(fcs_files19), length(fcs_files21), length(fcs_files23), length(fcs_files24))





