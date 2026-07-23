library(tidyverse)
library(readxl)

dat_path<- "/scratch/ckelsey4/Cayo_meth/cell_counts"

flow_metadata_22_23<- readRDS(paste(dat_path, "/flow_metadata.rds", sep =""))

paths<- list.files(dat_path, full.names = T)

ids_21<- read.csv(paste(dat_path, "/cayo_ids_21.csv", sep = ""), header = F, col.names = c("file_id", "monkey_id"))

trap_ids<- readRDS("/scratch/ckelsey4/trapping_ids.rds")

#Remove paths where data has NAs in sample column as it does not join properly
#And remove paths that are not cell count data
paths<- paths[!paths %in% c("/scratch/ckelsey4/Cayo_meth/cell_counts/2023-01-19", "/scratch/ckelsey4/Cayo_meth/cell_counts/cayo_ids_21.csv", 
                            "/scratch/ckelsey4/Cayo_meth/cell_counts/cayo_ids.txt", "/scratch/ckelsey4/Cayo_meth/cell_counts/flow_metadata.rds",
                            "/scratch/ckelsey4/Cayo_meth/cell_counts/lymphocyte_proportions.rds")]

#Reads in all files and binds them by date in a list
cc_list<- lapply(paths, function(x){
  
  files<- list.files(x ,full.names = T)
  
  dd<- lapply(files, read_xls)
  
  dd<- dd %>% 
    reduce(left_join, by = "...1") %>%
    rename(monkey_id = "...1")
  
  dd<- dd[!dd$monkey_id %in% c("Mean", "SD"),]
  
})


#Splits path names, retains date, and names list elements by date
cc_dates<- str_split_i(paths, "/", 6)
names(cc_list)<- cc_dates

cc_list<- lapply(names(cc_list), function(x){
  
  df<- cc_list[[x]]
  df$date<- as.Date(x)
  df
  
})

names(cc_list)<- cc_dates

df<- bind_rows(cc_list)

# Replace 2021 incorrect ids with correct ones
### Add ".mqd" suffix to ids df to match cell counts df
ids_21$file_id<- gsub(".fcs", ".mqd", ids_21$file_id)

### Remove common suffix to non-2021 monkey ids
df$monkey_id<- ifelse(!df$monkey_id %in% ids_21$file_id, str_split_i(df$monkey_id, "\\.",1), df$monkey_id)

### Use ids_21 df to change ids to correct 3-string ids
df$monkey_id<- ifelse(df$monkey_id %in% ids_21$file_id, 
                      ids_21$monkey_id, df$monkey_id)

df<- df %>%
  select(monkey_id, date, `Lymphocytes/Single Cells/CD3+ /CD4+ | Count`,
           `Lymphocytes/Single Cells/CD3+ /CD8+ | Count`,
           `Lymphocytes/CD3-CD16+ | Count`,
           `Lymphocytes/Single Cells/CD20+  | Count`)

colnames(df)<- c("animal_ID", "date", "cd3_cd4", "cd3_cd8", "cd3_cd16", "cd20")

colnames(df)<- gsub("Lymphocytes/", "", colnames(df))
colnames(df)<- gsub("Single Cells/", "", colnames(df))
colnames(df)<- gsub(" | Count", "", colnames(df))
colnames(df)<- gsub("|", "", colnames(df))

df<- df %>%
  mutate(lymph_count = cd3_cd4 + cd3_cd8 + cd3_cd16 + cd20,
         cd3_cd4_proportion = cd3_cd4/lymph_count,
         cd3_cd8_proportion = cd3_cd8/lymph_count,
         cd3_cd16_proportion = cd3_cd16/lymph_count,
         cd20_proportion = cd20/lymph_count)

df<- df %>%
  mutate(trap_year = case_when(
    between(date, as.Date("2021-10-01"), as.Date("2022-04-30")) ~ 2021,
    between(date, as.Date("2022-10-01"), as.Date("2023-04-30")) ~ 2022,
    between(date, as.Date("2023-10-01"), as.Date("2024-04-30")) ~ 2023,
    between(date, as.Date("2024-10-01"), as.Date("2025-04-30")) ~ 2023,
  ))

df<- left_join(df, trap_ids, by = c("monkey_id", "trap_year"))

saveRDS(df, "/scratch/ckelsey4/Cayo_meth/cell_counts/lymphocyte_proportions.rds")


