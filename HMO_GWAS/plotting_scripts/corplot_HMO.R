library(reshape2)
library(corrplot)

#d <- read.delim("~/work/UMCG/data/NEXT/HMO_GWAS/v5/results/hmo_gwas_imputed_pvals_for_corrplot.txt", sep = "\t", as.is = T, check.names = F)
d <- read.delim("~/work/UMCG/data/NEXT/HMO_GWAS/v5/results/hmo_gwas_nonimputed_pvals_for_corrplot.txt", sep = "\t", as.is = T, check.names = F)


# original x labels
d$time_locus <- factor(paste0(d$locus, ":",d$timepoint), levels = c("FUT2:W2","FUT2:M1","FUT2:M2","FUT2:M3","FUT2:M6","FUT3_FUT6:W2", "FUT3_FUT6:M1","FUT3_FUT6:M2","FUT3_FUT6:M3","FUT3_FUT6:M6","ABO:W2","ABO:M1","ABO:M2","ABO:M3","ABO:M6","ST3GAL6:M1","CADM1:M6","GRIN2B:M6"))
# No ABO M3
d$time_locus <- factor(paste0(d$locus, ":",d$timepoint), levels = c("FUT2:W2","FUT2:M1","FUT2:M2","FUT2:M3","FUT2:M6","FUT3_FUT6:W2", "FUT3_FUT6:M1","FUT3_FUT6:M2","FUT3_FUT6:M3","FUT3_FUT6:M6","ABO:W2","ABO:M1","ABO:M2","ABO:M3","ST3GAL6:M1","CADM1:M6","GRIN2B:M6"))
# All loci - tp combinations
d$time_locus <- factor(paste0(d$locus, ":",d$timepoint), levels = c("FUT2:W2","FUT2:M1","FUT2:M2","FUT2:M3","FUT2:M6","FUT3_FUT6:W2", "FUT3_FUT6:M1","FUT3_FUT6:M2","FUT3_FUT6:M3","FUT3_FUT6:M6","ABO:W2","ABO:M1","ABO:M2","ABO:M3","ABO:M6",
                                                                    "ST3GAL6:W2","ST3GAL6:M1","ST3GAL6:M2", "ST3GAL6:M3","ST3GAL6:M6", 
                                                                    "CADM1:W2","CADM1:M1","CADM1:M2", "CADM1:M3","CADM1:M6", 
                                                                    "GRIN2B:W2","GRIN2B:M1","GRIN2B:M2", "GRIN2B:M3","GRIN2B:M6"))
 
# only significant in non-imputed
d$time_locus <- factor(paste0(d$locus, ":",d$timepoint), levels = c("FUT2:W2","FUT2:M1","FUT2:M2","FUT2:M3","FUT2:M6", "FUT3_FUT6:M1","FUT3_FUT6:M2","FUT3_FUT6:M3","FUT3_FUT6:M6","ABO:M1","ABO:M2","ABO:M3","ST3GAL6:M1","CADM1:M6","GRIN2B:M6", "intergenic_9p23:M2"))


d$logp <- -1*log10(d$pval)


d_wide_logp <- dcast(d, hmo  ~ time_locus, value.var = "logp", fun.aggregate = max, na.rm = T)
d_wide_logp[d_wide_logp == -Inf] = 0
d_wide_logp[,"NA"] <- NULL
d_wide_logp$hmo<- gsub("_ugml_invr", "", d_wide_logp$hmo)
d_wide_logp$hmo<- gsub("mother_milk_HMO_", "", d_wide_logp$hmo)
hmo_order <- c("Total", #total HMOs
                 "Neut", "Fuc", "Sia", #grouped HMOs
                 "3GL","6GL","LNH","LNnT","LNT", #neutral HMOs
                 "2FL","3FL","A_tetra","DFLNHa","LDFT", #neutral+fucosylated HMOs
                 "LNDH_I","LNFP_I","LNFP_II","LNFP_III","LNFP_V",
                 "LNnDFH","LNnFP_V","MFLNH_III",
                 "3F3SL",
                 "3SL","6SL","DSLNT","LSTb","LSTc")

row.names(d_wide_logp) = d_wide_logp$hmo
d_wide_logp$hmo = NULL

bonf_cutoff <- -1*log10(1.785714e-09)
pdf("~/work/UMCG/data/NEXT/HMO_GWAS/v5/plots/HMO_corrplot_logp_v3_nonimputed_only_sign.v2.pdf", width = 8, height = 8)
corrplot(as.matrix(d_wide_logp[hmo_order,]), is.corr = F, tl.col = "black", pch.cex = 0.2,
         p.mat = as.matrix(d_wide_logp[hmo_order,]), insig = "pch", sig.level = c(bonf_cutoff), pch = 20)
dev.off()

