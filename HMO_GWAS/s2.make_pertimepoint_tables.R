
setwd("/Users/Dasha/work/UMCG/data/NEXT/HMO_GWAS/v5/")

timepoints <- c("W2", "M1", "M2", "M3", "M6")

mother_ids <- read.delim("230818_mother_IDs_for_GWAS_milk_HMO.txt", sep = "\t", check.names = F, as.is = T)
hmo <- read.delim("230818_milk_HMO_data_for_GWAS.txt", sep = "\t", check.names = F, as.is = T)

for (tp in timepoints){
  hmo_tp_full <- hmo[hmo$mother_sample_ID %in% mother_ids[mother_ids$time_point == tp,"mother_sample_ID"],]
  hmo_tp <- hmo_tp_full[,19:ncol(hmo_tp_full)]
  row.names(hmo_tp) <- gsub("[_;].*", "", hmo_tp_full$NEXT_ID_mother)
  
  hmo_tp <- hmo_tp %>% 
    rownames_to_column(var = "IID")
  #write filtered table
  hmo_tp <- cbind(a=0, hmo_tp)
  colnames(hmo_tp)[1] <- "#FID"
  cat(tp, dim(hmo_tp), "\n")
  write.table(hmo_tp, file = paste0("per_timepoint/HMO_", tp,".txt"), sep = "\t", quote = F, row.names = F)
}