
res <- data.frame(matrix(ncol = 5, nrow = 0))
cnt <- 1
for (tp in unique(infant_metadata$Timepoint_categorical)){
  meta_tp <- infant_metadata[infant_metadata$Timepoint_categorical == tp, ]
  duplicate_ids <- remove_one_of_duplicates(infant_metadata, duplicates, tp)
  tp_ids <- setdiff(row.names(meta_tp), duplicate_ids)
  
  tp_taxa <- extract_timepoint_convert_ids(infant_taxa_SGB_all_filt_CLR, tp_ids, tp)
  tp_taxa_bin <- extract_timepoint_convert_ids(infant_taxa_SGB_filt_binary, tp_ids, tp)
  tp_taxa_quant <- extract_timepoint_convert_ids(infant_taxa_SGB_all_filt_CLR_quant, tp_ids, tp)
  
  taxa <- tp_taxa
  for (sp in colnames(taxa)[3:ncol(taxa)]){
    merged <- inner_join(taxa[,c("IID", sp)], 
                       meta_tp[tp_ids,c("NEXT_ID", "clean_reads_FQ_1", "dna_conc", "exact_age_months_at_collection")],
                       by = c("IID" = "NEXT_ID")) %>% inner_join(.,fut2_family, by = c("IID" = "NEXT_ID.x"))
    colnames(merged) <- gsub(sp, "species", colnames(merged))
    mother_multiple <- names(which(table(merged$NEXT_ID.y) > 1) )
    sibs <- merged[merged$NEXT_ID.y %in% mother_multiple,c("IID", "NEXT_ID.y")]
    sibs_to_rm <- get_siblings_to_remove(sibs)
    merged <- merged[! merged$IID %in% sibs_to_rm,]
    

    lm_fit <- lm(species~fut2_inter  + clean_reads_FQ_1 + dna_conc + exact_age_months_at_collection, data = merged)
    coef <- summary(lm_fit)$coefficients
    if (min(coef[2:4,4]) < 0.01) {
      #cat(tp, sp, "\n")
      #print(summary(lm_fit))
      res <- rbind(res,c(tp, sp, coef[2:4,4]))
      cnt <- cnt + 1
    }
  }
  
}
colnames(res) <- c("tp", "sp", "p10", "p01", "p11")



get_siblings_to_remove <- function(sibs){
  
  ids_to_rm <- c()
  for (mother in unique(sibs[,2])){
    siblings <- sibs[sibs$NEXT_ID.y == mother,1]
    ids_to_rm <- c(ids_to_rm,siblings[2:length(siblings)])
 }
  return(ids_to_rm)
}
