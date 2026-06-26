###########################################################
# Colocalization v2: UGT1A1 Region (expanded)
# Using refined MoBa data + bilirubin GWAS
###########################################################

library(data.table)
library(coloc)

setwd("C:/tmp/gwas_data")

cat("========================================\n")
cat("STEP 1: Prepare bilirubin data\n")
cat("========================================\n")

bil <- fread("bilirubin_ugt1a1_region.tsv", header=TRUE)
cat("Bilirubin variants:", nrow(bil), "\n")

# Expand window to cover full UGT1A1 locus (chr2:234,200,000-234,700,000)
cat("Filtering to expanded UGT1A1 region...\n")
bil <- bil[bil$base_pair_location >= 234200000 & bil$base_pair_location <= 234700000, ]
cat("After filtering:", nrow(bil), "variants\n")

cat("\n========================================\n")
cat("STEP 2: Prepare MoBa fetal data\n")
cat("========================================\n")

# Read full MoBa and filter to UGT1A1 region
moba_full <- fread("moba-gwas-jaundice-fets.txt", sep="\t", header=TRUE,
                    select=c("CHR", "POS", "REF", "EFF", "EAF", "INFO", "N", 
                             "BETA", "SE", "LOG10P", "ID", "rsid", "nearestGene"),
                    nThread=4)

moba_region <- moba_full[moba_full$CHR == 2 & 
                          moba_full$POS >= 234200000 & 
                          moba_full$POS <= 234700000, ]
cat("MoBa variants in region:", nrow(moba_region), "\n")

# Remove duplicates - keep unique rsid+POS combinations
moba_region <- moba_region[!duplicated(moba_region[, .(rsid, POS)]), ]
cat("After dedup:", nrow(moba_region), "unique variants\n")

cat("\n========================================\n")
cat("STEP 3: Merge and QC\n")
cat("========================================\n")

# Rename for merge
setnames(bil, "base_pair_location", "POS")

merged <- merge(bil, moba_region, by="POS", suffixes=c("_bil","_moba"))
cat("Merged:", nrow(merged), "variants\n")

# QC: remove NAs and invalid values
merged$p_val <- as.numeric(merged$p_value_char)
valid <- !is.na(merged$beta) & !is.na(merged$standard_error) & 
         !is.na(merged$BETA) & !is.na(merged$SE) &
         merged$standard_error > 0 & merged$SE > 0 &
         !is.na(merged$effect_allele_frequency) & !is.na(merged$EAF) &
         !is.na(merged$p_val)

merged <- merged[valid, ]
cat("After QC:", nrow(merged), "variants\n")

if(nrow(merged) < 30) {
  cat("ERROR: Too few variants for coloc\n")
  # Try with all chr2 variants regardless of position
  cat("Falling back to all overlapping SNPs...\n")
  merged <- merge(bil, moba_region, by="POS", suffixes=c("_bil","_moba"))
  merged$p_val <- as.numeric(merged$p_value_char)
  merged <- merged[!is.na(merged$beta) & !is.na(merged$standard_error) & 
                    !is.na(merged$BETA) & !is.na(merged$SE) &
                    merged$standard_error > 0 & merged$SE > 0, ]
  cat("Fallback merged:", nrow(merged), "variants\n")
}

# Create SNP identifiers
merged$SNP <- merged$rsid

# Print stats
cat("\nRegion stats:\n")
cat(sprintf("  Total variants: %d\n", nrow(merged)))
cat(sprintf("  Bilirubin significant (p<1e-5): %d\n", sum(merged$p_val < 1e-5, na.rm=TRUE)))
cat(sprintf("  Jaundice significant (log10p>5): %d\n", sum(merged$LOG10P > 5, na.rm=TRUE)))

# Print key SNPs
cat("\nKey UGT1A1 SNPs:\n")
key_snps <- c("rs887829", "rs4148323", "rs6742078", "rs34352510", "rs11695484")
key <- merged[merged$rsid %in% key_snps, ]
if(nrow(key) > 0) {
  for(i in 1:nrow(key)) {
    cat(sprintf("  %s: bil_beta=%.4f p=%.2e | jaun_beta=%.4f log10p=%.1f\n",
                key$rsid[i], key$beta[i], key$p_val[i],
                key$BETA[i], key$LOG10P[i]))
  }
} else {
  cat("  None of the key UGT1A1 SNPs found in merged data\n")
}

cat("\n========================================\n")
cat("STEP 4: Colocalization\n")
cat("========================================\n")

# Prepare MAF
bil_eaf <- merged$effect_allele_frequency
bil_maf <- ifelse(bil_eaf < 0.5, bil_eaf, 1 - bil_eaf)
bil_maf[bil_maf < 0.01 | is.na(bil_maf)] <- 0.01

moba_eaf <- merged$EAF
moba_maf <- ifelse(moba_eaf < 0.5, moba_eaf, 1 - moba_eaf)
moba_maf[moba_maf < 0.01 | is.na(moba_maf)] <- 0.01

# Bilirubin = quantitative
coloc_d1 <- list(
  beta = merged$beta,
  varbeta = merged$standard_error^2,
  snp = merged$rsid,
  position = merged$POS,
  MAF = bil_maf,
  N = rep(436000, nrow(merged)),
  type = "quant"
)

# Jaundice = binary (case-control)
coloc_d2 <- list(
  beta = merged$BETA,
  varbeta = merged$SE^2,
  snp = merged$rsid,
  position = merged$POS,
  MAF = moba_maf,
  N = rep(max(merged$N, na.rm=TRUE), nrow(merged)),
  type = "cc",
  s = 0.25
)

cat("Running coloc.abf...\n")
coloc_result <- tryCatch({
  coloc.abf(coloc_d1, coloc_d2)
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

if(!is.null(coloc_result)) {
  cat("\n========== COLOCALIZATION RESULTS ==========\n")
  cat(sprintf("PP.H0 (no association):    %.4f\n", coloc_result$summary["PP.H0.abf"]))
  cat(sprintf("PP.H1 (bilirubin only):    %.4f\n", coloc_result$summary["PP.H1.abf"]))
  cat(sprintf("PP.H2 (jaundice only):     %.4f\n", coloc_result$summary["PP.H2.abf"]))
  cat(sprintf("PP.H3 (both, diff SNP):    %.4f\n", coloc_result$summary["PP.H3.abf"]))
  cat(sprintf("PP.H4 (shared causal SNP): %.4f\n", coloc_result$summary["PP.H4.abf"]))
  
  if(coloc_result$summary["PP.H4.abf"] > 0.8) {
    cat("\n✅ STRONG colocalization\n")
  } else if(coloc_result$summary["PP.H4.abf"] > 0.5) {
    cat("\n⚠️ Moderate colocalization\n")
  } else if(coloc_result$summary["PP.H3.abf"] > 0.5) {
    cat("\n📌 Two distinct causal variants\n")
  } else {
    cat("\n❌ Weak colocalization\n")
  }
  
  # Top SNPs
  results_df <- coloc_result$results
  top_idx <- order(results_df$SNP.PP.H4, decreasing=TRUE)[1:5]
  cat("\nTop 5 SNPs by PP.H4:\n")
  for(i in top_idx) {
    cat(sprintf("  %s PP.H4=%.4f\n", results_df$snp[i], results_df$SNP.PP.H4[i]))
  }
  
  # Save
  saveRDS(coloc_result, "colocalization_ugt1a1_final.rds")
  write.csv(data.frame(t(coloc_result$summary)), "colocalization_summary.csv", row.names=FALSE)
  write.csv(results_df, "colocalization_detail.csv", row.names=FALSE)
  cat("\n✅ Results saved.\n")
  
  # Plot
  pdf("figure_coloc_ugt1a1.pdf", width=10, height=8)
  par(mfrow=c(2,1))
  plot(coloc_result)
  dev.off()
  cat("Coloc plot saved.\n")
} else {
  cat("Coloc failed.\n")
}

cat("\n========================================\n")
cat("DONE\n")
cat("========================================\n")
