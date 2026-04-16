library(tidyverse)

#Import metadata----------------------------------------------------------------
long_data<- read.table("/scratch/ckelsey4/Cayo_meth/long_data_adjusted.txt")
pids<- long_data$lid_pid[long_data$university == 'uw']

#Import m/cov rds------------------------------------------------------------
prom_cov2<- readRDS("/scratch/ckelsey4/Cayo_meth/prom_cov_filtered")
prom_cov<- prom_cov[c(1:21)]
prom_m<- readRDS("/scratch/ckelsey4/Cayo_meth/prom_m_filtered")
prom_m<- prom_m[c(1:21)]

#Filter metadata to lids in regions list
long_data<- long_data[long_data$lid_pid %in% colnames(prom_cov[[1]]),]

prom_cov<- lapply(names(prom_cov), function(x){
  prom_cov<- subset(prom_cov[[x]], select=long_data$lid_pid)
  return(prom_cov)
})

prom_m<- lapply(names(prom_m), function(x){
  prom_m<- subset(prom_m[[x]], select=long_data$lid_pid)
  return(prom_m)
})

names(prom_cov)<- 1:21 #turn all chroms into integers (X = 21)
names(prom_m)<- 1:21 #turn all chroms into integers (X = 21)


samples_to_remove<- parallel::mclapply(names(prom_m), function(x) {
  
  df<- prom_m[[x]]
  
  res<- c()
  
  for (i in 1:nrow(df)) {
    
    row<- df[i, !(colnames(df) %in% pids)]
    
    if (sum(row) == 0) {
      loci<- rownames(row)
      res<- c(res, loci)
    }
  }
  
  return(res)
  
},mc.cores = 21)

names(samples_to_remove)<- 1:21

new_prom<- lapply(names(prom_m), function(x){
  
  prom<- prom_m[[x]]
  samples<- samples_to_remove[[x]]
  
  n_prom<- prom[!rownames(prom) %in% samples,]
  
  rows_diff<- nrow(prom) - nrow(n_prom)
  
  if (rows_diff == length(samples)){
    return(n_prom)
  } else {
    print("Length error")
  }
})

names(new_prom)<- 1:21

subset_cov <- function(original_list, filtered_list) {
  lapply(names(original_list), function(x) {
    
    orig <- original_list[[x]]
    filt <- filtered_list[[x]]
    
    orig[rownames(orig) %in% rownames(filt), , drop = FALSE]
  }) |> setNames(names(original_list))
}

new_cov<- subset_cov(prom_cov, new_prom)

saveRDS(new_prom, "/scratch/ckelsey4/Cayo_meth/prom_m_filtered2.rds")
saveRDS(new_cov, "/scratch/ckelsey4/Cayo_meth/prom_cov_filtered2.rds")

