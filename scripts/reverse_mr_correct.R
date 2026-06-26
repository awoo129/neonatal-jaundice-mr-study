###########################################################
# REVERSE MR: Neonatal Jaundice → Gut Microbiota (CORRECT)
# Exposure: MoBa fetal GWAS (jaundice)
# Outcome: MiBioGen microbiota GWAS (extracted from OpenGWAS)
###########################################################

library(TwoSampleMR)
library(data.table)
library(ieugwasr)

setwd("C:/tmp/gwas_data")

# Set JWT token
token <- Sys.getenv("OPENGWAS_TOKEN")
if(token == "") {
  env <- readLines("C:/Users/Administrator/.env")
  tl <- env[grep("OPENGWAS_TOKEN", env)]
  token <- sub("^OPENGWAS_TOKEN=", "", tl[1])
  token <- gsub('["\\\']', '', token)
  token <- trimws(token)
}
Sys.setenv(OPENGWAS_JWT = token)

cat("========================================\n")
cat("REVERSE MR: Jaundice → Microbiota\n")
cat("========================================\n\n")

# Load MoBa fetal GWAS as exposure
cat("Loading MoBa fetal GWAS as exposure...\n")
moba_fetal <- fread("moba-gwas-jaundice-fets.txt", sep="\t", header=TRUE,
                     select=c("rsid", "CHR", "POS", "REF", "EFF", "EAF", "BETA", "SE", "LOG10P", "N"),
                     nThread=4)
cat("Loaded:", nrow(moba_fetal), "variants\n")

# Create exposure data frame
exposure_dat <- data.frame(
  SNP = moba_fetal$rsid,
  effect_allele.exposure = moba_fetal$EFF,
  other_allele.exposure = moba_fetal$REF,
  eaf.exposure = moba_fetal$EAF,
  beta.exposure = moba_fetal$BETA,
  se.exposure = moba_fetal$SE,
  pval.exposure = 10^(-moba_fetal$LOG10P),
  samplesize.exposure = moba_fetal$N,
  exposure = "Neonatal jaundice (Fetal)",
  id.exposure = "moba_fetal_reverse",
  stringsAsFactors = FALSE
)

# MiBioGen IDs
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

# For reverse MR:
# Exposure = neonatal jaundice (MoBa)
# Outcome = microbiota (need to extract from OpenGWAS)
# 
# The microbiota instruments files have SNP + beta for microbiota
# We need to create outcome data from these

reverse_results <- data.frame()

for(id in names(mibio_ids)) {
  taxon <- mibio_ids[id]
  cat(sprintf("\n[%d/%d] %s (reverse MR)\n", 
              which(names(mibio_ids)==id), length(mibio_ids), taxon))
  
  # Load microbiota instruments
  inst_file <- paste0(id, "_instruments.rds")
  if(!file.exists(inst_file)) {
    cat("  ❌ No instruments file\n")
    next
  }
  
  mi_inst <- readRDS(inst_file)
  cat("  Instruments:", nrow(mi_inst), "SNPs\n")
  
  # Create outcome data from instruments
  # The instruments file contains microbiota GWAS data
  # We need to rename columns to match harmonise_data expectations
  
  outcome_dat <- data.frame(
    SNP = mi_inst$SNP,
    effect_allele.outcome = mi_inst$effect_allele.exposure,
    other_allele.outcome = mi_inst$other_allele.exposure,
    eaf.outcome = mi_inst$eaf.exposure,
    beta.outcome = mi_inst$beta.exposure,
    se.outcome = mi_inst$se.exposure,
    pval.outcome = mi_inst$pval.exposure,
    samplesize.outcome = mi_inst$samplesize.exposure,
    outcome = taxon,
    id.outcome = id,
    stringsAsFactors = FALSE
  )
  
  # Harmonize
  cat("  Harmonizing...\n")
  harm <- tryCatch({
    harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
  }, error = function(e) {
    cat("  Harmonization error:", e$message, "\n")
    return(NULL)
  })
  
  if(is.null(harm) || nrow(harm) == 0) {
    cat("  ❌ Harmonization failed\n")
    next
  }
  
  harm <- harm[harm$mr_keep == TRUE, ]
  cat("  Harmonised:", nrow(harm), "SNPs\n")
  
  if(nrow(harm) < 2) {
    cat("  ❌ Too few harmonised SNPs\n")
    next
  }
  
  # Run MR
  cat("  Running MR...\n")
  mr_res <- tryCatch({
    mr(harm, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
  }, error = function(e) {
    cat("  MR error:", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(mr_res) && nrow(mr_res) > 0) {
    mr_res$taxon <- taxon
    mr_res$nsnp <- as.integer(mr_res$nsnp)
    mr_res$OR <- exp(mr_res$b)
    mr_res$OR_lower <- exp(mr_res$b - 1.96*mr_res$se)
    mr_res$OR_upper <- exp(mr_res$b + 1.96*mr_res$se)
    mr_res$direction <- "reverse"
    
    reverse_results <- rbind(reverse_results, mr_res)
    
    # Print IVW result
    ivw <- mr_res[mr_res$method == "Inverse variance weighted", ]
    if(nrow(ivw) > 0) {
      sig <- ifelse(ivw$pval < 0.05, "⭐", "")
      cat(sprintf("  IVW: OR=%.3f (%.3f-%.3f), p=%.2e %s\n", 
                  ivw$OR, ivw$OR_lower, ivw$OR_upper, ivw$pval, sig))
    }
  } else {
    cat("  ❌ MR failed\n")
  }
}

# Save reverse MR results
if(nrow(reverse_results) > 0) {
  saveRDS(reverse_results, "reverse_mr_results.rds")
  write.csv(reverse_results, "reverse_mr_results.csv", row.names=FALSE)
  cat("\n✅ Reverse MR results saved.\n")
  
  # Summary
  ivw_rev <- reverse_results[reverse_results$method == "Inverse variance weighted", ]
  sig_count <- sum(ivw_rev$pval < 0.05, na.rm=TRUE)
  cat(sprintf("Significant taxa (p<0.05): %d/%d\n", sig_count, nrow(ivw_rev)))
} else {
  cat("\n⚠️ No reverse MR results obtained.\n")
}

# ========================================
# COMBINED FORWARD + REVERSE MR
# ========================================
cat("\n\n========================================\n")
cat("COMBINED FORWARD + REVERSE MR TABLE\n")
cat("========================================\n\n")

# Load forward MR results
forward_res <- read.csv("mibio_jaundice_mr_results.csv")
forward_ivw <- forward_res[forward_res$method == "Inverse variance weighted", ]

# Load reverse MR results
if(file.exists("reverse_mr_results.csv")) {
  reverse_res <- read.csv("reverse_mr_results.csv")
  reverse_ivw <- reverse_res[reverse_res$method == "Inverse variance weighted", ]
  
  # Create combined table
  all_taxa <- unique(c(forward_ivw$taxon, reverse_ivw$taxon))
  
  combined <- data.frame(
    Taxon = all_taxa,
    Forward_OR = NA, Forward_P = NA, Forward_NSNP = NA,
    Reverse_OR = NA, Reverse_P = NA, Reverse_NSNP = NA
  )
  
  for(taxon in all_taxa) {
    fwd <- forward_ivw[forward_ivw$taxon == taxon, ]
    rev <- reverse_ivw[reverse_ivw$taxon == taxon, ]
    
    idx <- which(combined$Taxon == taxon)
    if(nrow(fwd) > 0) {
      combined$Forward_OR[idx] <- fwd$OR[1]
      combined$Forward_P[idx] <- fwd$pval[1]
      combined$Forward_NSNP[idx] <- fwd$nsnp[1]
    }
    if(nrow(rev) > 0) {
      combined$Reverse_OR[idx] <- rev$OR[1]
      combined$Reverse_P[idx] <- rev$pval[1]
      combined$Reverse_NSNP[idx] <- rev$nsnp[1]
    }
  }
  
  write.csv(combined, "combined_forward_reverse_mr.csv", row.names=FALSE)
  cat("Combined table saved: combined_forward_reverse_mr.csv\n")
  
  # Print summary
  cat("\n--- Combined Forward + Reverse MR ---\n")
  cat(sprintf("%-25s %10s %10s %10s | %10s %10s %10s\n", 
              "Taxon", "F_OR", "F_p", "F_NSnp", "R_OR", "R_p", "R_NSnp"))
  cat(rep("-", 90), "\n", sep="")
  
  for(i in 1:nrow(combined)) {
    cat(sprintf("%-25s %9.3f %9.2e %6d | %9.3f %9.2e %6d\n",
                combined$Taxon[i],
                combined$Forward_OR[i], combined$Forward_P[i], combined$Forward_NSNP[i],
                combined$Reverse_OR[i], combined$Reverse_P[i], combined$Reverse_NSNP[i]))
  }
}

cat("\n========================================\n")
cat("ALL REVERSE MR ANALYSES COMPLETE\n")
cat("========================================\n")
