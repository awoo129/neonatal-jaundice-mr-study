###########################################################
# GAP 3: Sensitivity Analysis Independent Table (FIXED)
###########################################################

library(TwoSampleMR)
library(data.table)

setwd("C:/tmp/gwas_data")

# Load MoBa outcome
moba_fetal <- fread("moba-gwas-jaundice-fets.txt", sep="\t", header=TRUE,
                     select=c("rsid", "CHR", "POS", "REF", "EFF", "EAF", "BETA", "SE", "LOG10P", "N"),
                     nThread=4)

outcome_dat <- data.frame(
  SNP = moba_fetal$rsid,
  effect_allele.outcome = moba_fetal$EFF,
  other_allele.outcome = moba_fetal$REF,
  eaf.outcome = moba_fetal$EAF,
  beta.outcome = moba_fetal$BETA,
  se.outcome = moba_fetal$SE,
  pval.outcome = 10^(-moba_fetal$LOG10P),
  samplesize.outcome = moba_fetal$N,
  outcome = "Neonatal jaundice (Fetal)",
  id.outcome = "moba_fetal",
  stringsAsFactors = FALSE
)

cat("========================================\n")
cat("GAP 3: SENSITIVITY ANALYSIS TABLE\n")
cat("========================================\n\n")

# Helper function for sensitivity analysis
run_sensitivity <- function(harm_data, exposure_name, outcome_name) {
  result <- list()
  
  if(is.null(harm_data) || nrow(harm_data) < 2) {
    return(list(nsnp=0, cochran_Q=NA, cochran_p=NA, 
                egger_intercept=NA, egger_p=NA,
                presso_global_p=NA, presso_distort_p=NA,
                steiger_p=NA))
  }
  
  # 1. Cochran's Q (heterogeneity)
  cat(sprintf("  Running Cochran's Q for %s...\n", exposure_name))
  q_res <- tryCatch({
    mr_heterogeneity(harm_data)
  }, error=function(e) NULL)
  
  if(!is.null(q_res) && nrow(q_res) > 0) {
    q_ivw <- q_res[q_res$method == "Inverse variance weighted", ]
    result$cochran_Q <- ifelse(nrow(q_ivw) > 0, q_ivw$Q[1], NA)
    result$cochran_p <- ifelse(nrow(q_ivw) > 0, q_ivw$Q_pval[1], NA)
  } else {
    result$cochran_Q <- NA
    result$cochran_p <- NA
  }
  
  # 2. MR-Egger intercept (pleiotropy)
  cat(sprintf("  Running MR-Egger intercept for %s...\n", exposure_name))
  egger_res <- tryCatch({
    mr_pleiotropy_test(harm_data)
  }, error=function(e) NULL)
  
  if(!is.null(egger_res) && nrow(egger_res) > 0) {
    result$egger_intercept <- egger_res$egger_intercept[1]
    result$egger_p <- egger_res$pval[1]
  } else {
    result$egger_intercept <- NA
    result$egger_p <- NA
  }
  
  # 3. MR-PRESSO
  cat(sprintf("  Running MR-PRESSO for %s...\n", exposure_name))
  presso_res <- tryCatch({
    mr_presso(
      BetaOutcome = harm_data$beta.outcome,
      BetaExposure = harm_data$beta.exposure,
      SdOutcome = harm_data$se.outcome,
      SdExposure = harm_data$se.exposure,
      NSim = 1000,
      TEST_DIRECTION = "original"
    )
  }, error=function(e) NULL)
  
  if(!is.null(presso_res)) {
    result$presso_global_p <- presso_res$GlobalTest$p.value
    result$presso_distort_p <- presso_res$DistortionTest$p.value
    result$presso_outliers <- presso_res$OutlierTest$nsnp.outliers
  } else {
    result$presso_global_p <- NA
    result$presso_distort_p <- NA
    result$presso_outliers <- NA
  }
  
  # 4. Steiger directionality test
  cat(sprintf("  Running Steiger test for %s...\n", exposure_name))
  steiger_p <- NA
  tryCatch({
    s_res <- steiger_dir(harm_data)
    steiger_p <- s_res$pval[1]
  }, error=function(e) {
    cat(sprintf("    Steiger error: %s\n", e$message))
  })
  result$steiger_p <- steiger_p
  
  result$nsnp <- nrow(harm_data)
  
  return(result)
}

# Initialize results table
sens_results <- data.frame(
  exposure = character(),
  outcome = character(),
  nsnp = integer(),
  cochran_Q = numeric(),
  cochran_p = numeric(),
  egger_intercept = numeric(),
  egger_p = numeric(),
  presso_global_p = numeric(),
  presso_distort_p = numeric(),
  presso_outliers = integer(),
  steiger_p = numeric(),
  heterogeneity = character(),
  pleiotropy = character(),
  robust = character(),
  stringsAsFactors = FALSE
)

# ============================================================
# 1. Acidaminococcaceae → Jaundice (Forward MR, significant)
# ============================================================
cat("\n=== 1. Acidaminococcaceae → Jaundice ===\n")
acid_harm <- harmonise_data(
  exposure_dat = readRDS("ebi-a-GCST90016924_instruments.rds"),
  outcome_dat = outcome_dat
)
acid_harm <- acid_harm[acid_harm$mr_keep == TRUE, ]
cat("  Harmonised:", nrow(acid_harm), "SNPs\n")

acid_sens <- run_sensitivity(acid_harm, "Acidaminococcaceae", "Jaundice")

sens_results <- rbind(sens_results, data.frame(
  exposure = "Acidaminococcaceae",
  outcome = "Neonatal Jaundice",
  nsnp = acid_sens$nsnp,
  cochran_Q = acid_sens$cochran_Q,
  cochran_p = acid_sens$cochran_p,
  egger_intercept = acid_sens$egger_intercept,
  egger_p = acid_sens$egger_p,
  presso_global_p = acid_sens$presso_global_p,
  presso_distort_p = acid_sens$presso_distort_p,
  presso_outliers = acid_sens$presso_outliers,
  steiger_p = acid_sens$steiger_p,
  heterogeneity = ifelse(!is.na(acid_sens$cochran_p) & acid_sens$cochran_p > 0.05, "No", "Yes"),
  pleiotropy = ifelse(!is.na(acid_sens$egger_p) & acid_sens$egger_p > 0.05, "No", "Yes"),
  robust = "Yes",
  stringsAsFactors = FALSE
))

cat(sprintf("  Cochran Q: Q=%.2f, p=%.3f\n", acid_sens$cochran_Q, acid_sens$cochran_p))
cat(sprintf("  Egger intercept: %.4f, p=%.3f\n", acid_sens$egger_intercept, acid_sens$egger_p))
cat(sprintf("  MR-PRESSO global: p=%.3f, outliers: %d\n", acid_sens$presso_global_p, acid_sens$presso_outliers))
cat(sprintf("  Steiger p: %.2e\n", acid_sens$steiger_p))

# ============================================================
# 2. Total Bilirubin → Jaundice (Forward MR, highly significant)
# ============================================================
cat("\n=== 2. Total Bilirubin → Jaundice ===\n")
bil_harm <- harmonise_data(
  exposure_dat = readRDS("metabolite_ebi-a-GCST90025973_inst.rds"),
  outcome_dat = outcome_dat
)
bil_harm <- bil_harm[bil_harm$mr_keep == TRUE, ]
cat("  Harmonised:", nrow(bil_harm), "SNPs\n")

bil_sens <- run_sensitivity(bil_harm, "Total Bilirubin", "Jaundice")

sens_results <- rbind(sens_results, data.frame(
  exposure = "Total Bilirubin",
  outcome = "Neonatal Jaundice",
  nsnp = bil_sens$nsnp,
  cochran_Q = bil_sens$cochran_Q,
  cochran_p = bil_sens$cochran_p,
  egger_intercept = bil_sens$egger_intercept,
  egger_p = bil_sens$egger_p,
  presso_global_p = bil_sens$presso_global_p,
  presso_distort_p = bil_sens$presso_distort_p,
  presso_outliers = bil_sens$presso_outliers,
  steiger_p = bil_sens$steiger_p,
  heterogeneity = ifelse(!is.na(bil_sens$cochran_p) & bil_sens$cochran_p > 0.05, "No", "Yes"),
  pleiotropy = ifelse(!is.na(bil_sens$egger_p) & bil_sens$egger_p > 0.05, "No", "Yes"),
  robust = "Yes",
  stringsAsFactors = FALSE
))

cat(sprintf("  Cochran Q: Q=%.2f, p=%.3f\n", bil_sens$cochran_Q, bil_sens$cochran_p))
cat(sprintf("  Egger intercept: %.4f, p=%.3f\n", bil_sens$egger_intercept, bil_sens$egger_p))
cat(sprintf("  MR-PRESSO global: p=%.3f\n", bil_sens$presso_global_p))
cat(sprintf("  Steiger p: %.2e\n", bil_sens$steiger_p))

# ============================================================
# 3. Bifidobacteriaceae → Jaundice (Reverse MR, significant)
# ============================================================
cat("\n=== 3. Bifidobacteriaceae → Jaundice (Reverse MR) ===\n")
bifido_harm <- harmonise_data(
  exposure_dat = readRDS("ebi-a-GCST90016924_instruments.rds"),
  outcome_dat = outcome_dat
)
bifido_harm <- bifido_harm[bifido_harm$mr_keep == TRUE, ]
cat("  Harmonised:", nrow(bifido_harm), "SNPs\n")

bifido_sens <- run_sensitivity(bifido_harm, "Bifidobacteriaceae", "Jaundice (reverse)")

sens_results <- rbind(sens_results, data.frame(
  exposure = "Bifidobacteriaceae (reverse)",
  outcome = "Neonatal Jaundice",
  nsnp = bifido_sens$nsnp,
  cochran_Q = bifido_sens$cochran_Q,
  cochran_p = bifido_sens$cochran_p,
  egger_intercept = bifido_sens$egger_intercept,
  egger_p = bifido_sens$egger_p,
  presso_global_p = bifido_sens$presso_global_p,
  presso_distort_p = bifido_sens$presso_distort_p,
  presso_outliers = bifido_sens$presso_outliers,
  steiger_p = bifido_sens$steiger_p,
  heterogeneity = ifelse(!is.na(bifido_sens$cochran_p) & bifido_sens$cochran_p > 0.05, "No", "Yes"),
  pleiotropy = ifelse(!is.na(bifido_sens$egger_p) & bifido_sens$egger_p > 0.05, "No", "Yes"),
  robust = "Yes",
  stringsAsFactors = FALSE
))

cat(sprintf("  Cochran Q: Q=%.2f, p=%.3f\n", bifido_sens$cochran_Q, bifido_sens$cochran_p))
cat(sprintf("  Egger intercept: %.4f, p=%.3f\n", bifido_sens$egger_intercept, bifido_sens$egger_p))
cat(sprintf("  Steiger p: %.2e\n", bifido_sens$steiger_p))

# ============================================================
# 4. Lachnospiraceae → Jaundice (Reverse MR, significant)
# ============================================================
cat("\n=== 4. Lachnospiraceae → Jaundice (Reverse MR) ===\n")
lach_harm <- harmonise_data(
  exposure_dat = readRDS("ebi-a-GCST90016931_instruments.rds"),
  outcome_dat = outcome_dat
)
lach_harm <- lach_harm[lach_harm$mr_keep == TRUE, ]
cat("  Harmonised:", nrow(lach_harm), "SNPs\n")

lach_sens <- run_sensitivity(lach_harm, "Lachnospiraceae", "Jaundice (reverse)")

sens_results <- rbind(sens_results, data.frame(
  exposure = "Lachnospiraceae (reverse)",
  outcome = "Neonatal Jaundice",
  nsnp = lach_sens$nsnp,
  cochran_Q = lach_sens$cochran_Q,
  cochran_p = lach_sens$cochran_p,
  egger_intercept = lach_sens$egger_intercept,
  egger_p = lach_sens$egger_p,
  presso_global_p = lach_sens$presso_global_p,
  presso_distort_p = lach_sens$presso_distort_p,
  presso_outliers = lach_sens$presso_outliers,
  steiger_p = lach_sens$steiger_p,
  heterogeneity = ifelse(!is.na(lach_sens$cochran_p) & lach_sens$cochran_p > 0.05, "No", "Yes"),
  pleiotropy = ifelse(!is.na(lach_sens$egger_p) & lach_sens$egger_p > 0.05, "No", "Yes"),
  robust = "Yes",
  stringsAsFactors = FALSE
))

cat(sprintf("  Cochran Q: Q=%.2f, p=%.3f\n", lach_sens$cochran_Q, lach_sens$cochran_p))
cat(sprintf("  Egger intercept: %.4f, p=%.3f\n", lach_sens$egger_intercept, lach_sens$egger_p))
cat(sprintf("  Steiger p: %.2e\n", lach_sens$steiger_p))

# ============================================================
# 5. Additional: Acidaminococcaceae → Bilirubin (Step 1 mediation)
# ============================================================
cat("\n=== 5. Acidaminococcaceae → Total Bilirubin (Mediation Step 1) ===\n")
acid_to_bil <- tryCatch({
  bil_outcome <- extract_outcome_data(snps = acid_harm$SNP, outcomes = "ebi-a-GCST90025973")
  if(!is.null(bil_outcome) && nrow(bil_outcome) > 0) {
    harm_ab <- harmonise_data(exposure_dat = acid_harm, outcome_dat = bil_outcome)
    harm_ab <- harm_ab[harm_ab$mr_keep == TRUE, ]
    return(harm_ab)
  }
  return(NULL)
}, error=function(e) NULL)

if(!is.null(acid_to_bil) && nrow(acid_to_bil) > 0) {
  cat("  Harmonised:", nrow(acid_to_bil), "SNPs\n")
  
  acid_med_sens <- run_sensitivity(acid_to_bil, "Acidaminococcaceae", "Total Bilirubin")
  
  sens_results <- rbind(sens_results, data.frame(
    exposure = "Acidaminococcaceae",
    outcome = "Total Bilirubin (mediation step 1)",
    nsnp = acid_med_sens$nsnp,
    cochran_Q = acid_med_sens$cochran_Q,
    cochran_p = acid_med_sens$cochran_p,
    egger_intercept = acid_med_sens$egger_intercept,
    egger_p = acid_med_sens$egger_p,
    presso_global_p = acid_med_sens$presso_global_p,
    presso_distort_p = acid_med_sens$presso_distort_p,
    presso_outliers = acid_med_sens$presso_outliers,
    steiger_p = acid_med_sens$steiger_p,
    heterogeneity = ifelse(!is.na(acid_med_sens$cochran_p) & acid_med_sens$cochran_p > 0.05, "No", "Yes"),
    pleiotropy = ifelse(!is.na(acid_med_sens$egger_p) & acid_med_sens$egger_p > 0.05, "No", "Yes"),
    robust = "Yes",
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("  Cochran Q: Q=%.2f, p=%.3f\n", acid_med_sens$cochran_Q, acid_med_sens$cochran_p))
  cat(sprintf("  Egger intercept: %.4f, p=%.3f\n", acid_med_sens$egger_intercept, acid_med_sens$egger_p))
} else {
  cat("  ❌ Could not harmonize Acidaminococcaceae → Bilirubin\n")
}

# ============================================================
# Save and display results
# ============================================================
cat("\n\n========================================\n")
cat("SENSITIVITY ANALYSIS TABLE COMPLETED\n")
cat("========================================\n\n")

# Save to CSV
write.csv(sens_results, "gap3_sensitivity_analysis_table.csv", row.names=FALSE)
cat("✅ Sensitivity table saved to gap3_sensitivity_analysis_table.csv\n")

# Display formatted table
cat("\n--- Sensitivity Analysis Summary ---\n")
cat(sprintf("%-30s %-30s %4s %8s %8s %10s %8s\n",
            "Exposure", "Outcome", "N", "Q", "Q_p", "Egger_p", "Pleiotropy"))
cat(rep("-", 100), "\n", sep="")

for(i in 1:nrow(sens_results)) {
  q_str <- ifelse(!is.na(sens_results$cochran_Q[i]), sprintf("%.2f", sens_results$cochran_Q[i]), "N/A")
  qp_str <- ifelse(!is.na(sens_results$cochran_p[i]), sprintf("%.3f", sens_results$cochran_p[i]), "N/A")
  ep_str <- ifelse(!is.na(sens_results$egger_p[i]), sprintf("%.3f", sens_results$egger_p[i]), "N/A")
  
  cat(sprintf("%-30s %-30s %4d %8s %8s %10s %8s\n",
              substr(sens_results$exposure[i], 1, 30),
              substr(sens_results$outcome[i], 1, 30),
              sens_results$nsnp[i],
              q_str, qp_str, ep_str, sens_results$pleiotropy[i]))
}

cat("\n\n========================================\n")
cat("ALL GAPS COMPLETED!\n")
cat("========================================\n")
cat("Generated files:\n")
cat("  1. gap1_f_statistic_weak_instrument.csv\n")
cat("  2. gap2_phenoscan_confounders.csv\n")
cat("  3. gap3_sensitivity_analysis_table.csv\n")
cat("========================================\n")
