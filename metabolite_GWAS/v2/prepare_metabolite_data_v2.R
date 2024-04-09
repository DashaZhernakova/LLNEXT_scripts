setwd("/Users/Dasha/work/UMCG/data/NEXT/metabolites/v3/data/")

clr_data <- read.delim("20230206_clr_mtb_above_70.txt.gz", sep = "\t", check.names = F, as.is = T, row.names = 1)


samples <- read.delim("genotype_ids.txt", header = F, as.is = T)

d <- d[d$NEXT_ID_geno %in% samples[,1],]

save_dedup <- function(d, who, when){
  d_timepoint <- d[d$Plasma_Sample_type == who & d$Plasma_Timepoint == when,]
  
  counts <- table(d_timepoint$NEXT_ID_geno)
  duplicates <- names(counts[counts > 1])
  
  rows_to_remove <- which(d_timepoint$NEXT_ID_geno %in% duplicates, )[-1]
  cat("Removing ", length(rows_to_remove), "rows\n")
  if (length(rows_to_remove) > 0) {
    dedup <- d_timepoint[-rows_to_remove,]
  } else {
    dedup <- d_timepoint
  }
  dedup$NEXT_ID <- dedup$NEXT_ID_geno
  dedup <- subset(dedup, select = -c(PARENT_SAMPLE_NAME, Family_ID, NEXT_ID_geno ,Plasma_Sample_type, Plasma_Timepoint))
  colnames(dedup)[1] = "#IID"
  write.table(dedup, file = paste0("metabo_", who, "_", when, ".filtered.txt"), sep = "\t", quote = F, row.names = F)
  
  row.names(dedup) <- dedup[,"#IID"]
  dedup <- subset(dedup, select=-`#IID`)
  return (dedup)
}

table(d[,c("Plasma_Sample_type", "Plasma_Timepoint")])
dedup = NULL
dedup <- save_dedup(d, "Mother", "B")
#adjust_for_PCs_sex(dedup, "Mother", "B")

dedup = NULL
dedup <- save_dedup(d, "Mother", "P12")
#adjust_for_PCs_sex(dedup, "Mother", "P12")

dedup = NULL
dedup <- save_dedup(d, "Mother", "P28")
#adjust_for_PCs_sex(dedup,"Mother", "P28")

dedup = NULL
dedup <- save_dedup(d, "Baby", "B")
#adjust_for_PCs_sex(dedup, "Baby", "B")
