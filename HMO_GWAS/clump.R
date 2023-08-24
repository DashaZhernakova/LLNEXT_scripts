args <- commandArgs(trailingOnly = TRUE)
library(ieugwasr)

res_path <- args[1]
d <- read.delim(res_path, header = F,  sep = "\t", as.is = T, check.names = F, row.names = NULL)

colnames(d)[5] <- "rsid"
colnames(d)[15] <- "pval"

res2 <-  ld_clump(d, clump_r2 = 0.1, clump_kb = 250000)

write.table(res2, file = paste0(res_path, ".clumped_0.1.txt"), sep = "\t", quote = F, row.names = F)
