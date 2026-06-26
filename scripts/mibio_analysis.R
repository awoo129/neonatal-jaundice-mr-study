###########################################################
# MiBioGen Microbiota → Neonatal Jaundice MR (full analysis)
# All 14 taxa with Bonferroni correction
###########################################################

library(TwoSampleMR)
library(data.table)
library(ieugwasr)

setwd("C:/tmp/gwas_data")

# Set JWT
token <- Sys.getenv("OPENGWAS_TOKEN")
if(token == "") {
  env <- readLines("C:/Users/Administrator/.env")
  tl <- env[grep("OPENGWAS_TOKEN", env)]
  token <- sub("^OPENGWAS_TOKEN=", "", tl[1])
  token <- gsub('["\\\']', '', token)
  token <- trimws(token)
}
Sys.setenv(OPENGWAS_JWT = token)

# MiBioGen taxa IDs and names
mibio_ids <- c(
  "ebi-a-GCST90016921" = "Acidaminococcaceae",
  "ebi-a-GCST90016922" = "Actinomycetaceae",
  "ebi-a-GCST90016923" = "Bacteroidaceae",
  "ebi-a-GCST90016924" = "Bifidobacteriaceae",
  "ebi-a-GCST90016925" = "Christensenellaceae",
  "ebi-a-GCST90016926" = "Clostridiaceae 1", 
  "ebi-a-GCST90016927" = "Coriobacteriaceae",
  "ebi-a-GCST90016928" = "Desulfovibrionaceae",
  "ebi-a-GCST90016929" = "Enterobacteriaceae",
  "ebi-a-GCST90016931" = "Lachnospiraceae",
  "ebi-a-GCST90016932" = "Lactobacillaceae",
  "ebi-a-GCST90016933" = "Prevotellaceae",
  "ebi-a-GCST90016934" = "Ruminococcaceae",
  "ebi-a-GCST90016935" = "Tannerellaceae"
)

# Load MoBa fetal GWAS (for outcome)
cat("Loading MoBa fetal GWAS...\n")
moba_fetal <- fread("moba-gwas-jaundice-fets.txt", sep="\t", header=TRUE,
                     select=c("rsid", "CHR", "POS", "REF", "EFF", "EAF", "BETA", "SE", "LOG10P", "N"),
                     nThread=4)
cat("Loaded:", nrow(moba_fetal), "variants\n")

# Bonferroni threshold
bonf_thresh <- 0.05 / length(mibio_ids)
cat(sprintf("Bonferroni threshold: %.6f\n\n", bonf_thresh))

# Initialize results
all_results <- data.frame()

for(id in names(mibio_ids)) {
  taxon <- mibio_ids[id]
  cat(sprintf("\n═══════════════ %s (%s) ═══════════════\n", taxon, id))
  
  # Extract instruments
  inst_file <- paste0(id, "_instruments.rds")
  if(file.exists(inst_file)) {
    exp_dat <- readRDS(inst_file)
    cat("Loaded existing instruments:", nrow(exp_dat), "SNPs\n")
  } else {
    cat("Extracting instruments from OpenGWAS...\n")
    exp_dat <- tryCatch({
      extract_instruments(outcomes = id, p1 = 1e-5, r2 = 0.001, kb = 10000,
                          clump = TRUE, access_token = token)
    }, error = function(e) {
      cat("  Error:", e$message, "\n")
      # Try without clumping
      extract_instruments(outcomes = id, p1 = 1e-5, clump = FALSE,
                          access_token = token)
    })
    if(!is.null(exp_dat) && nrow(exp_dat) > 0) {
      saveRDS(exp_dat, inst_file)
      cat("Extracted:", nrow(exp_dat), "SNPs\n")
    }
  }
  
  if(is.null(exp_dat) || nrow(exp_dat) < 2) {
    cat("❌  Insufficient instruments\n")
    next
  }
  
  # Match outcome
  outcome_variants <- exp_dat$SNP
  moba_sub <- moba_fetal[moba_fetal$rsid %in% outcome_variants, ]
  
  if(nrow(moba_sub) < 2) {
    cat("❌  Too few overlapping SNPs with MoBa\n")
    next
  }
  
  outcome_dat <- data.frame(
    SNP = moba_sub$rsid,
    effect_allele.outcome = moba_sub$EFF,
    other_allele.outcome = moba_sub$REF,
    eaf.outcome = moba_sub$EAF,
    beta.outcome = moba_sub$BETA,
    se.outcome = moba_sub$SE,
    pval.outcome = 10^(-moba_sub$LOG10P),
    samplesize.outcome = moba_sub$N,
    outcome = "Neonatal jaundice (Fetal)",
    id.outcome = "moba_fetal",
    stringsAsFactors = FALSE
  )
  
  # Harmonize
  harm <- harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome_dat)
  harm <- harm[harm$mr_keep == TRUE, ]
  
  if(nrow(harm) < 2) {
    cat("❌  Too few after harmonisation\n")
    next
  }
  
  # Run MR
  mr_res <- mr(harm, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
  
  # Add OR
  mr_res$OR <- exp(mr_res$b)
  mr_res$OR_lower <- exp(mr_res$b - 1.96*mr_res$se)
  mr_res$OR_upper <- exp(mr_res$b + 1.96*mr_res$se)
  mr_res$taxon <- taxon
  mr_res$nsnp <- as.integer(mr_res$nsnp)
  
  # Add heterogeneity test
  het <- mr_heterogeneity(harm)
  if(any(het$method == "Inverse variance weighted")) {
    het_ivw <- het[het$method == "Inverse variance weighted", ]
    mr_res$Q <- het_ivw$Q[1]
    mr_res$Q_pval <- het_ivw$Q_pval[1]
  }
  
  # Add pleiotropy test
  pleio <- mr_pleiotropy_test(harm)
  mr_res$egger_intercept <- pleio$egger_intercept[1]
  mr_res$egger_pval <- pleio$pval[1]
  
  # Add Bonferroni status
  mr_res$bonf_sig <- mr_res$pval < bonf_thresh
  
  # Store
  all_results <- rbind(all_results, mr_res)
  
  # Print summary
  ivw_row <- mr_res[mr_res$method == "Inverse variance weighted", ]
  if(nrow(ivw_row) > 0) {
    cat(sprintf("  IVW: OR=%.3f (%.3f-%.3f), p=%.2e %s\n", 
                ivw_row$OR, ivw_row$OR_lower, ivw_row$OR_upper, ivw_row$pval,
                ifelse(ivw_row$bonf_sig, "🔥", "")))
  }
  wm_row <- mr_res[mr_res$method == "Weighted median", ]
  if(nrow(wm_row) > 0) {
    cat(sprintf("  WM:  OR=%.3f, p=%.2e\n", wm_row$OR, wm_row$pval))
  }
}

cat("\n\n========================================\n")
cat("COMPLETE MR ANALYSIS: ALL 14 TAXA\n")
cat("========================================\n\n")

# Print final table
cat(sprintf("%-25s %8s %8s %8s %10s %6s %10s\n", "Taxon", "N_SNP", "OR", "CI_low", "CI_up", "P", "Bonf"))
cat(rep("=", 80), "\n", sep="")

for(i in 1:nrow(all_results)) {
  if(all_results$method[i] == "Inverse variance weighted") {
    cat(sprintf("%-25s %6d %7.3f %7.3f-%7.3f %9.2e %5s %9.2f\n",
                all_results$taxon[i], all_results$nsnp[i],
                all_results$OR[i], all_results$OR_lower[i],
                all_results$OR_upper[i], all_results$pval[i],
                ifelse(all_results$bonf_sig[i], "🔥", ""),
                ifelse(!is.na(all_results$Q[i]), all_results$Q[i], NA)))
  }
}

# Save all results
write.csv(all_results, "mibio_jaundice_mr_results.csv", row.names=FALSE)
cat("\n✅ Results saved to mibio_jaundice_mr_results.csv\n")
