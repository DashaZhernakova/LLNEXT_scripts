library(mmrm)
library(car)
setwd("/Users/Dasha/work/UMCG/data/NEXT/mbQTLs/")

taxa <- infant_taxa_SGB_all_filt_CLR_quant
taxa <- infant_taxa_SGB_all_filt_CLR_2
taxa <- read.delim("data/infant_taxa_SGB_all_filt_CLR_zero_to_min.txt", header = F, check.names = F, as.is = T, sep = "\t")
qtls <- read.delim("data/quant_all_species_and_timepoints.1e-05.srt.pos_clumped.txt", header = F, check.names = F, as.is = T, sep = "\t")

geno <- as.data.frame(t(read.delim("data/snps_1e-5.genotypes.txt", as.is = T, sep = "\t", check.names = F, row.names = 1)))
ids <- row.names(geno)
geno <- as.data.frame(sapply(geno, function(x) as.factor(as.character(x))))
geno$"NEXT_ID" <- ids

pdf("plots/mmrm_quant_examples_v5.pdf")
#results <- c("W2-Escherichia_coli-17:56518534", "M1-Streptococcus_parasanguinis-5:52285098", "M1-Streptococcus_parasanguinis-4:189900854", "M2-Bifidobacterium_bifidum-8:13938718", "M3-Parabacteroides_distasonis-7:77836451", "M6-Bifidobacterium_bifidum-10:124514229", "M12-Phocaeicola_vulgatus-20:41719099")
for (i in 1:nrow(qtls)){
  #res = results[1]
  tp = qtls[i,1]
  sp = qtls[i,2]
  snp = qtls[i,4]
  
  # prepare and merge data for plotting
  meta_subs <- infant_metadata[,c("NEXT_ID", "NG_ID", "Timepoint_categorical", "clean_reads_FQ_1", "dna_conc", "exact_age_months_at_collection")]
  colnames(meta_subs) <- c("NEXT_ID", "NG_ID", "Timepoint", "num_reads", "dna_conc", "age")
  meta_subs$Timepoint<- factor(meta_subs$Timepoint, levels = c("W2", "M1", "M2", "M3", "M6", "M12"))
  taxa_subs <- as.data.frame(taxa[,sp])
  colnames(taxa_subs) <- "species"
  taxa_subs$NG_ID <- row.names(taxa)
  
  merged <- left_join(meta_subs, taxa_subs, by="NG_ID") %>%
    left_join(.,geno[,c("NEXT_ID",snp)], by = "NEXT_ID")
  colnames(merged)[ncol(merged)] <- "snp"
  merged <- na.omit(merged)
  
  
  #ggplot(merged, aes(x = snp, y = species)) + geom_boxplot(alpha = 0.1) + geom_point(color = "dodgerblue4", alpha = 0.5)+ facet_grid(. ~ Timepoint_categorical) + theme_bw()
    
  
  merged$tp_p <-NA
  
  
  ids_nodup <- c()
  tp_levels <- c()
  for (tp in levels(meta_subs$Timepoint)){
    meta_tp <- infant_metadata[infant_metadata$Timepoint_categorical == tp, ]
    duplicate_ids <- remove_one_of_duplicates(infant_metadata, duplicates, tp)
    tp_ids <- setdiff(row.names(meta_tp), duplicate_ids)
    ids_nodup <- c(ids_nodup,tp_ids)
    tp_taxa <- extract_timepoint_convert_ids(taxa, tp_ids, tp)
    
    merged_tp <- left_join(meta_tp[,c("NEXT_ID", "NG_ID", "Timepoint_categorical", "clean_reads_FQ_1", "dna_conc", "exact_age_months_at_collection")], 
                           tp_taxa[,c("IID",sp)], by=c("NEXT_ID" = "IID")) %>%
      left_join(.,geno[,c("NEXT_ID",snp)], by = "NEXT_ID")
    colnames(merged_tp)[ncol(merged_tp)] <- "snp"
    colnames(merged_tp)[ncol(merged_tp) - 1] <- "species"
    merged_tp <- na.omit(merged_tp)
    merged_tp$snp <- as.numeric(factor(merged_tp$snp))
    tryCatch({
      lm_fit <- lm(species ~ snp + clean_reads_FQ_1 + dna_conc + exact_age_months_at_collection,data = merged_tp)  
      p <- summary(lm_fit)$coefficients["snp", 4]
      merged[merged$Timepoint == tp, "tp_p"] <- paste0(tp, "\nP = ", formatC(p, digits = 2))
      tp_levels <- c(tp_levels, paste0(tp, "\nP = ", formatC(p, digits = 2)))
    }, error=function(e) {p = NA})
   }
  
  merged <- merged[merged$NG_ID %in% ids_nodup,]
  merged$snp_num <- as.numeric(factor(merged$snp))
  merged$snp <- as.factor(merged$snp)
  
  #ggplot(merged, aes(x = age, y = species, color = snp)) + geom_point() + geom_smooth(method = "lm")
  
  ll_p = 1
  tryCatch({
    cat ("fitting mmrm\n")
    #mmrm_fit <- mmrm(species ~ age + I(age^2) +I(age^3) + snp_num*age + num_reads + dna_conc + us(Timepoint|NEXT_ID), data=merged, reml = F)
    #mmrm_fit_0 <- mmrm(species ~ snp_num + age + I(age^2) +I(age^3) + num_reads + dna_conc + us(Timepoint|NEXT_ID), data=merged, reml = F)
    mmrm_fit <- mmrm(species ~ age + I(age^2) +I(age^3) + snp_num*Timepoint + num_reads + dna_conc + us(Timepoint|NEXT_ID), data=merged, reml = F)
    mmrm_fit_0 <- mmrm(species ~ snp_num + Timepoint + age + I(age^2) +I(age^3) + num_reads + dna_conc + us(Timepoint|NEXT_ID), data=merged, reml = F)
    
    ll_res <- run_likelihoodRatioTest(mmrm_fit_0, mmrm_fit)
    ll_p <- formatC(ll_res[2,5], digits = 3)
    cat("Log-Likelihood test P = ", ll_p, "\n")
  }, error=function(e) {ll_p = 1}) 
  
  if (as.numeric(ll_p) < 1e-6){
    print ("Worked!")
    merged$tp_p <- factor(merged$tp_p, levels = tp_levels)
    p <- ggplot(merged, aes(x = snp, y = species)) + 
      geom_boxplot(alpha = 0.1, outlier.shape = NA) + 
      geom_jitter(color = "dodgerblue4", alpha = 0.5, width = 0.2) + 
      facet_grid(. ~ tp_p) +
      ggtitle(paste0(snp, " - ", sp, "\nmmrm LL test P = ", ll_p)) + 
      theme_bw()
    p2 <- ggplot(merged, aes(x = Timepoint, y = species, color = snp)) + geom_boxplot(alpha = 0.4) +
      scale_color_manual(values = c("#F65C00", "#006DB3", "#00A072")) + theme_bw() +
      ggtitle(paste0("mmrm LL test P = ", ll_p)) +
      ylab(sp) + labs(color = snp)
    print(p)
    print(p2)
    #ggplot(merged, aes(x = snp, y = species, color = Timepoint)) + geom_boxplot(alpha = 0.4) + geom_smooth(aes(x = snp_num, y = species, color = Timepoint), method = 'lm', formula = y ~ x, alpha = 0.2) 
  }
}
dev.off()

run_likelihoodRatioTest = function(model1, model2) {
  class_checkpoint = ("mmrm" %in% class(model1) & "mmrm" %in% class(model2))
  reml_checkpoint = (model1$reml == F & model2$reml == F)
  if (class_checkpoint == FALSE) stop ("One or both models do not belong to MMRM class. Execution stopped!")
  if (reml_checkpoint == FALSE) stop ("One or both models use REML for calculation. The results will be unreliable. Execution stopped!")
  result.matrix = matrix(rep(NA,10),ncol = 5)
  colnames(result.matrix) <- c("#Df", "LogLik", "Df", "Chisq", "Pr(>Chisq)")
  rownames(result.matrix) <- c("Model1","Model2")
  result.matrix[1,1] <- length(model1$beta_est)
  result.matrix[2,1] <- length(model2$beta_est)
  result.matrix[1,2] <- model1$neg_log_lik
  result.matrix[2,2] = model2$neg_log_lik
  result.matrix[1,3] = 0
  result.matrix[2,3] = length(model2$beta_est) - length(model1$beta_est)
  result.matrix[1,4] = 0
  result.matrix[2,4] = 2*abs(result.matrix[2,2] - result.matrix[1,2])
  result.matrix[1,5] = NA
  result.matrix[2,5] = pchisq(result.matrix[2,4], round(abs(result.matrix[2,3])),lower.tail = F)
  result.matrix 
}


                                                              