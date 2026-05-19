library(tidyverse)

#Import metadata----------------------------------------------------------------
blood_metadata<- read.table("/scratch/ckelsey4/Cayo_meth/blood_metadata_full.txt", sep = "\t", header = T)
long_data<- read.table("/scratch/ckelsey4/Cayo_meth/long_data_adjusted.txt")

long_ids<- unique(long_data$monkey_id)

overlap_lids<- blood_metadata[blood_metadata$monkey_id %in% long_ids,]

overlap_lids<- overlap_lids %>%
  group_by(monkey_id) %>%
  sample_n(1)

lids_to_remove<- long_data[!long_data$lid_pid %in% overlap_lids$lid_pid,]
blood_metadata<- blood_metadata[!blood_metadata$lid_pid %in% lids_to_remove$lid_pid,]

blood_metadata<- blood_metadata %>%
  filter(age_at_sampling > 1)

pids<- blood_metadata$lid_pid[blood_metadata$university == 'uw']

rm(lids_to_remove);rm(long_ids);rm(overlap_lids);rm(long_data)

#Import m/cov rds------------------------------------------------------------
# load region lists that have been filtered for 5x coverage in 90% of samples
regions_cov<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_cov_filtered.rds")
regions_m<- readRDS("/scratch/ckelsey4/Cayo_meth/regions_m_filtered.rds")

#Filter metadata to lids in regions list
blood_metadata<- blood_metadata[blood_metadata$lid_pid %in% colnames(regions_cov[[1]]),]

regions_cov<- lapply(names(regions_cov), function(x){
  regions_cov<- subset(regions_cov[[x]], select=blood_metadata$lid_pid)
  return(regions_cov)
})

regions_m<- lapply(names(regions_m), function(x){
  regions_m<- subset(regions_m[[x]], select=blood_metadata$lid_pid)
  return(regions_m)
})

names(regions_cov)<- 1:21 #turn all chroms into integers (X = 21)
names(regions_m)<- 1:21 #turn all chroms into integers (X = 21)

samples_to_remove<- parallel::mclapply(names(regions_m), function(x) {
  
  df<- regions_m[[x]]
  
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

new_regions<- lapply(names(regions_m), function(x){
  
  regions<- regions_m[[x]]
  samples<- samples_to_remove[[x]]
  
  n_regions<- regions[!rownames(regions) %in% samples,]
  
  rows_diff<- nrow(regions) - nrow(n_regions)
  
  if (rows_diff == length(samples)){
    return(n_regions)
  } else {
    print("Length error")
  }
})

names(new_regions)<- 1:21

subset_cov <- function(original_list, filtered_list) {
  lapply(names(original_list), function(x) {
    
    orig <- original_list[[x]]
    filt <- filtered_list[[x]]
    
    orig[rownames(orig) %in% rownames(filt), , drop = FALSE]
  }) |> setNames(names(original_list))
}

new_cov<- subset_cov(regions_cov, new_regions)

saveRDS(new_regions, "/scratch/ckelsey4/Cayo_meth/regions_m_filtered_cs.rds")
saveRDS(new_cov, "/scratch/ckelsey4/Cayo_meth/regions_cov_filtered_cs.rds")

