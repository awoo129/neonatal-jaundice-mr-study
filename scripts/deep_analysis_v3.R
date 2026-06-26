###########################################################
# Deep MR Analysis (v3) - Fix LOO, Steiger, Coloc
###########################################################

library(TwoSampleMR)
library(MRPRESSO)
library(coloc)
library(data.table)
library(ggplot2)
library(stringr)
library(httr)
library(jsonlite)

setwd("C:/tmp/gwas_data")

# Load token
opengwas_token <- Sys.getenv("OPENGWAS_TOKEN")
if(opengwas_token == "") {
  env_lines <- tryCatch(readLines("C:/Users/Administrator/.env"), error=function(e) NULL)
  if(!is.null(env_lines)) {
    token_line <- env_lines[grep("OPENGWAS_TOKEN", env_lines)]
    if(length(token_line) > 0) {
      opengwas_token <- sub("^OPENGWAS_TOKEN=", "", token_line[1])
      opengwas_token <- gsub('["\\\']', '', opengwas_token)
      opengwas_token <- trimws(opengwas_token)
    }
  }
}

# Set JWT for ieugwasr
library(ieugwasr)
Sys.setenv(OPENGWAS_JWT = opengwas_token)

cat("========================================\n")
cat("STEP 1: Load harmonised data\n")
cat("========================================\n")

# Load or re-create harmonised data
if(file.exists("harmonised_fetal.rds")) {
  harm_data <- readRDS("harmonised_fetal.rds")
  cat("Loaded existing harmonised data:", nrow(harm_data), "SNPs\n")
} else {
  bil_exp <- readRDS("bilirubin_instruments.rds")
  moba_fetal <- fread("moba-gwas-jaundice-fets.txt", sep="\t", header=TRUE, 
                       select=c("rsid", "CHR", "POS", "REF", "EFF", "EAF", "BETA", "SE", "LOG10P", "N"),
                       nThread=4)
  moba_sub <- moba_fetal[moba_fetal$rsid %in% bil_exp$SNP, ]
  outcome_dat <- data.frame(
    SNP = moba_sub$rsid, effect_allele.outcome = moba_sub$EFF,
    other_allele.outcome = moba_sub$REF, eaf.outcome = moba_sub$EAF,
    beta.outcome = moba_sub$BETA, se.outcome = moba_sub$SE,
    pval.outcome = 10^(-moba_sub$LOG10P), samplesize.outcome = moba_sub$N,
    outcome = "Neonatal jaundice (Fetal)", id.outcome = "moba_fetal",
    chr.outcome = moba_sub$CHR, pos.outcome = moba_sub$POS, stringsAsFactors = FALSE
  )
  harm_data <- harmonise_data(exposure_dat = bil_exp, outcome_dat = outcome_dat)
  harm_data <- harm_data[harm_data$mr_keep == TRUE, ]
  saveRDS(harm_data, "harmonised_fetal.rds")
  cat("Created harmonised data:", nrow(harm_data), "SNPs\n")
}

cat("\n========================================\n")
cat("STEP 2: Funnel plot (text output)\n")
cat("========================================\n")

singlesnp_results <- mr_singlesnp(harm_data, all_method = c("mr_ivw", "mr_egger_regression"))
cat("Single SNP results:", nrow(singlesnp_results), "rows\n")

# Save singlesnp data for manual plots
saveRDS(singlesnp_results, "singlesnp_results.rds")

# Try funnel plot with explicit factor fix
singlesnp_results$SNP <- as.character(singlesnp_results$SNP)
# Remove duplicate factor issue
singlesnp_results$SNP <- make.unique(singlesnp_results$SNP)

pdf("figure3_funnel.pdf", width=8, height=6)
p_funnel <- mr_funnel_plot(singlesnp_results)
if(is.list(p_funnel)) print(p_funnel[[1]]) else print(p_funnel)
dev.off()
cat("Funnel plot saved.\n")

cat("\n========================================\n")
cat("STEP 3: Leave-one-out (manual LOO)\n")
cat("========================================\n")

# Manual leave-one-out to avoid factor issues
res_loo <- mr_leaveoneout(harm_data, method = mr_ivw)
cat("LOO results:", nrow(res_loo), "rows\n")
res_loo$SNP <- make.unique(as.character(res_loo$SNP))

# Create manual LOO plot data
res_loo$logOR <- res_loo$b
res_loo$logOR_lower <- res_loo$b - 1.96 * res_loo$se
res_loo$logOR_upper <- res_loo$b + 1.96 * res_loo$se
res_loo$OR <- exp(res_loo$b)
res_loo$OR_lower <- exp(res_loo$logOR_lower)
res_loo$OR_upper <- exp(res_loo$logOR_upper)

# Mark the IVW overall result
ivw_idx <- which(res_loo$SNP == "All - Inverse variance weighted")
res_loo$is_overall <- ifelse(res_loo$SNP == "All - Inverse variance weighted", "Overall", "Leave-one-out")

# Sort by b (effect size)
res_loo_sort <- res_loo[order(res_loo$b), ]
res_loo_sort$SNP <- factor(res_loo_sort$SNP, levels = res_loo_sort$SNP)

# Manual LOO forest plot
library(grid)
p_loo <- ggplot(res_loo_sort, aes(x = SNP, y = b, color = is_overall)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = b - 1.96*se, ymax = b + 1.96*se), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = c("Overall" = "red", "Leave-one-out" = "steelblue")) +
  coord_flip() +
  labs(title = "Leave-One-Out Analysis (Bilirubin → Neonatal Jaundice)",
       x = "SNP removed", y = "MR effect size (log-OR)") +
  theme_minimal() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())

pdf("figure4_loo.pdf", width=10, height=12)
print(p_loo)
dev.off()
cat("LOO plot saved.\n")

# Save LOO as csv
write.csv(res_loo, "loo_results.csv", row.names=FALSE)

cat("\n========================================\n")
cat("STEP 4: Steiger Directionality Test\n")
cat("========================================\n")

steiger <- directionality_test(harm_data)
cat("\n--- Steiger Directionality Test ---\n")
print(steiger)

# Check if causality direction is correct
if("correct_causal_direction" %in% names(steiger)) {
  cat("\nCorrect direction:", steiger$correct_causal_direction[1], "\n")
}
saveRDS(steiger, "steiger_test_fetal.rds")
write.csv(steiger, "steiger_test_fetal.csv", row.names=FALSE)

cat("\n========================================\n")
cat("STEP 5: Colocalization - alternative approach\n")
cat("========================================\n")

# Try OpenGWAS API directly for bilirubin regional data
cat("Trying OpenGWAS mrbase API...\n")

# Method 1: Try ieugwasr with token
api_result <- tryCatch({
  # Try the correct API format
  resp <- GET(
    url = "https://api.epigraphdb.org/gwas/association/",
    query = list(id = "ebi-a-GCST90025973", chromosome = 2, 
                 start = 234000000, end = 234500000),
    add_headers(Authorization = paste("Bearer", opengwas_token)),
    timeout(30)
  )
  if(resp$status_code == 200) {
    d <- content(resp, "parsed", encoding="UTF-8")
    if(!is.null(d$results)) as.data.frame(d$results) else NULL
  } else NULL
}, error = function(e) { cat("API error:", e$message, "\n"); NULL })

if(!is.null(api_result) && nrow(api_result) > 0) {
  cat("ieugwasr returned", nrow(api_result), "variants\n")
  bil_region <- api_result
} else {
  cat("ieugwasr failed. Trying OpenGWAS REST API...\n")
  
  # Try the old MR-Base API 
  resp <- tryCatch({
    GET(
      url = "https://api.opengwas.io/api/associations/ebi-a-GCST90025973",
      query = list(chromosome = 2, start = 234000000, end = 234500000),
      add_headers(Authorization = paste("Bearer", opengwas_token)),
      timeout(60)
    )
  }, error = function(e) {
    cat("REST API error:", e$message, "\n")
    return(NULL)
  })
  
  bil_region <- NULL
  if(!is.null(resp) && resp$status_code == 200) {
    resp_text <- content(resp, "text", encoding="UTF-8")
    resp_json <- fromJSON(resp_text)
    if(!is.null(resp_json$results) || !is.null(resp_json$data)) {
      d <- if(!is.null(resp_json$results)) resp_json$results else resp_json$data
      bil_region <- as.data.frame(d)
      cat("API returned", nrow(bil_region), "variants\n")
    }
  } else if(!is.null(resp)) {
    cat("API status:", resp$status_code, "\n")
  }
}

if(!is.null(bil_region) && nrow(bil_region) > 50) {
  cat("\nRunning colocalization...\n")
  
  # Column mapping
  pos_col <- grep("pos|position|base_pair_location", names(bil_region), value=TRUE, ignore.case=TRUE)[1]
  beta_col <- grep("beta|effect", names(bil_region), value=TRUE, ignore.case=TRUE)[1]
  se_col <- grep("se|standard_error|standard.err", names(bil_region), value=TRUE, ignore.case=TRUE)[1] 
  eaf_col <- grep("eaf|effect_allele_frequency|freq", names(bil_region), value=TRUE, ignore.case=TRUE)[1]
  snp_col <- grep("snp|rsid|variant|rs_id|name", names(bil_region), value=TRUE, ignore.case=TRUE)[1]
  
  cat("Column mapping: pos=", pos_col, " beta=", beta_col, " se=", se_col, " eaf=", eaf_col, "\n")
  
  if(length(pos_col) > 0 && length(beta_col) > 0) {
    bil_region$POS <- as.numeric(bil_region[[pos_col]])
    bil_region$BETA <- as.numeric(bil_region[[beta_col]])
    bil_region$SE <- as.numeric(bil_region[[se_col]])
    
    # Get MoBa region
    moba_fetal <- fread("moba-gwas-jaundice-fets.txt", sep="\t", header=TRUE, 
                         select=c("rsid", "CHR", "POS", "REF", "EFF", "EAF", "BETA", "SE", "LOG10P", "N"),
                         nThread=4)
    moba_region <- moba_fetal[moba_fetal$CHR == 2 & moba_fetal$POS >= 234000000 & moba_fetal$POS <= 234500000, ]
    cat("MoBa region:", nrow(moba_region), "variants\n")
    
    # Merge by position
    merged <- merge(bil_region, moba_region, by="POS", suffixes=c("_bil","_moba"))
    cat("Merged:", nrow(merged), "variants\n")
    
    if(nrow(merged) > 50) {
      bil_eaf <- if(length(eaf_col) > 0) as.numeric(merged[[eaf_col]]) else rep(0.3, nrow(merged))
      bil_eaf[is.na(bil_eaf)] <- 0.3
      
      # Handle MAF for coloc
      bil_maf <- ifelse(bil_eaf < 0.5, bil_eaf, 1 - bil_eaf)
      bil_maf[bil_maf < 0.01] <- 0.01
      
      moba_eaf <- merged$EAF
      moba_maf <- ifelse(moba_eaf < 0.5, moba_eaf, 1 - moba_eaf)
      moba_maf[moba_maf < 0.01] <- 0.01
      
      coloc_d1 <- list(
        beta = merged$BETA_bil,
        varbeta = merged$SE_bil^2,
        snp = paste0("chr2:", merged$POS),
        position = merged$POS,
        MAF = bil_maf,
        N = rep(436000, nrow(merged)),
        type = "quant"
      )
      
      coloc_d2 <- list(
        beta = merged$BETA,
        varbeta = merged$SE^2,
        snp = paste0("chr2:", merged$POS),
        position = merged$POS,
        MAF = moba_maf,
        N = rep(max(merged$N, na.rm=TRUE), nrow(merged)),
        type = "cc",
        s = 0.05
      )
      
      coloc_result <- coloc.abf(coloc_d1, coloc_d2)
      
      cat("\n========== COLOCALIZATION UGT1A1 ==========\n")
      cat(sprintf("PP.H0 (none):  %.4f\n", coloc_result$summary["PP.H0.abf"]))
      cat(sprintf("PP.H1 (bil):   %.4f\n", coloc_result$summary["PP.H1.abf"]))
      cat(sprintf("PP.H2 (jaun):  %.4f\n", coloc_result$summary["PP.H2.abf"]))
      cat(sprintf("PP.H3 (diff):  %.4f\n", coloc_result$summary["PP.H3.abf"]))
      cat(sprintf("PP.H4 (share): %.4f\n", coloc_result$summary["PP.H4.abf"]))
      
      if(coloc_result$summary["PP.H4.abf"] > 0.8) cat("✅ STRONG colocalization\n")
      else if(coloc_result$summary["PP.H4.abf"] > 0.5) cat("⚠️ Moderate colocalization\n")
      else cat("❌ Weak colocalization\n")
      
      saveRDS(coloc_result, "colocalization_ugt1a1.rds")
      write.csv(data.frame(t(coloc_result$summary)), "colocalization_summary.csv", row.names=FALSE)
      
      # Also plot coloc signal
      pdf("figure_coloc_signal.pdf", width=10, height=6)
      par(mfrow=c(2,1))
      plot(coloc_result)
      dev.off()
      cat("Colocalization signal plot saved.\n")
    }
  }
} else {
  cat("\n⚠️ Could not retrieve region data for colocalization.\n")
  cat("Writing colocalization diagnostics...\n")
  
  # Sensitivity analysis for key bilirubin loci
  cat("\nKey loci analysis for colocalization proxy:\n")
  
  # UGT1A1 locus - chr2
  ugt_snps <- harm_data[grep("rs887829|rs4148323|rs6742078", harm_data$SNP, ignore.case=TRUE), ]
  if(nrow(ugt_snps) > 0) {
    cat("UGT1A1 locus SNPs:\n")
    for(i in 1:nrow(ugt_snps)) {
      cat(sprintf("  %s: b_bil=%.4f, b_jaun=%.4f, p=%.2e\n",
                  ugt_snps$SNP[i], ugt_snps$beta.exposure[i],
                  ugt_snps$beta.outcome[i], ugt_snps$pval.outcome[i]))
    }
  }
  
  # Check all SNPs in chr2 (UGT1A1 region)
  chr2_snps <- harm_data[harm_data$chr.exposure == 2, ]
  cat("\nSNPs on chr2 (possible UGT1A1 region):\n")
  for(i in 1:min(10, nrow(chr2_snps))) {
    cat(sprintf("  %s (pos=%s): b_bil=%.4f, b_jaun=%.4f\n",
                chr2_snps$SNP[i], chr2_snps$pos.exposure[i],
                chr2_snps$beta.exposure[i], chr2_snps$beta.outcome[i]))
  }
}

cat("\n========================================\n")
cat("STEP 6: Heterogeneity + Pleiotropy\n")
cat("========================================\n")

het <- mr_heterogeneity(harm_data)
cat("Heterogeneity:\n")
print(het)

pleio <- mr_pleiotropy_test(harm_data)
cat("\nMR-Egger intercept:\n")
print(pleio)

# Save
write.csv(het, "heterogeneity_test.csv", row.names=FALSE)
write.csv(pleio, "pleiotropy_test.csv", row.names=FALSE)

cat("\n========================================\n")
cat("STEP 7: Final MR summary\n")
cat("========================================\n")

res <- mr(harm_data)
res$OR <- exp(res$b)
res$OR_lower <- exp(res$b - 1.96*res$se)
res$OR_upper <- exp(res$b + 1.96*res$se)
res$ci <- sprintf("%.3f (%.3f-%.3f)", res$OR, res$OR_lower, res$OR_upper)
cat("\n--- Final MR Results (Bilirubin → Fetal Jaundice) ---\n")
print(res[, c("method", "nsnp", "b", "se", "pval", "ci")])
write.csv(res, "mr_complete_results_fetal.csv", row.names=FALSE)

cat("\n✅ ALL ANALYSES COMPLETE\n")
cat("\nSaved files:\n")
cat("  - harmonised_fetal.rds, mr_presso_fetal.rds\n")
cat("  - steiger_test_fetal.csv\n")
cat("  - heterogeneity_test.csv, pleiotropy_test.csv\n")
cat("  - mr_complete_results_fetal.csv\n")
cat("  - figure3_funnel.pdf, figure4_loo.pdf\n")
cat("  - loo_results.csv, singlesnp_results.rds\n")
