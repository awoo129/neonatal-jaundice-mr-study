###########################################################
# CORRECTED BIDIRECTIONAL MR + TWO-STEP MR
# Part A: Reverse MR (Jaundice → Microbiota)
# Part B: Two-Step MR (Microbiota → Metabolite → Jaundice)
# Part C: Functional Pathway MR
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

# ========================================
# PART A: REVERSE MR (Jaundice → Microbiota)
# ========================================
cat("========================================\n")
cat("PART A: REVERSE MR\n")
cat("========================================\n\n")

# Load MiBioGen instruments
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

# Load MoBa fetal GWAS as exposure (neonatal jaundice)
cat("Loading MoBa fetal GWAS as exposure...\n")
moba_fetal <- fread("moba-gwas-jaundice-fets.txt", sep="\t", header=TRUE,
                     select=c("rsid", "CHR", "POS", "REF", "EFF", "EAF", "BETA", "SE", "LOG10P", "N"),
                     nThread=4)
cat("Loaded:", nrow(moba_fetal), "variants\n")

# Create exposure data frame with CORRECT column names
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

# Save exposure data
saveRDS(exposure_dat, "exposure_jaundice_reverse.rds")
cat("Exposure data saved.\n")

# Run reverse MR for each microbiota taxon
reverse_results <- data.frame()

for(id in names(mibio_ids)) {
  taxon <- mibio_ids[id]
  cat(sprintf("\n[%d/%d] %s (reverse MR)\n", 
              which(names(mibio_ids)==id), length(mibio_ids), taxon))
  
  # Load microbiota instruments (these are outcome data in reverse MR)
  inst_file <- paste0(id, "_instruments.rds")
  if(!file.exists(inst_file)) {
    cat("  ❌ No instruments file\n")
    next
  }
  
  mi_dat <- readRDS(inst_file)
  cat("  Instruments:", nrow(mi_dat), "SNPs\n")
  
  # Harmonize: exposure = jaundice, outcome = microbiota
  harm <- tryCatch({
    harmonise_data(exposure_dat = exposure_dat, outcome_dat = mi_dat)
  }, error = function(e) {
    cat("  Harmonization error:", e$message, "\n")
    return(NULL)
  })
  
  if(is.null(harm)) {
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
# PART B: TWO-STEP MR (Metabolite Mediation)
# ========================================
cat("\n\n========================================\n")
cat("PART B: TWO-STEP MR - Metabolite Mediation\n")
cat("========================================\n\n")

# Search for metabolite GWAS from OpenGWAS
cat("Searching for metabolite GWAS...\n")

# Try specific metabolite IDs from OpenGWAS
metabolite_ids <- c(
  "ebi-a-GCST90025973" = "Total Bilirubin",
  "ebi-a-GCST90025974" = "Direct Bilirubin",
  "ebi-a-GCST90025975" = "Indirect Bilirubin"
)

# Also search for bile acid and other metabolite GWAS
cat("Searching for metabolite traits...\n")
met_search <- tryCatch({
  ieugwasr::search("metabolite OR bile acid OR cholesterol OR lipid", limit=30)
}, error=function(e) {
  cat("Search error:", e$message, "\n")
  return(NULL)
})

if(!is.null(met_search) && nrow(met_search) > 0) {
  cat("Found", nrow(met_search), "metabolite GWAS datasets\n")
  
  # Filter for relevant ones
  relevant_met <- met_search[
    grepl("bilirubin|bile|cholesterol|lipid|metabolite|serum|plasma", 
          met_search$trait, ignore.case=TRUE), ]
  
  cat("Relevant metabolites:", nrow(relevant_met), "\n")
  
  # Run two-step MR for top metabolites
  for(idx in 1:min(10, nrow(relevant_met))) {
    met_id <- relevant_met$id[idx]
    met_trait <- relevant_met$trait[idx]
    cat(sprintf("\n--- Metabolite: %s ---\n", substr(met_trait, 1, 50)))
    
    # Step 1: Extract metabolite instruments
    met_inst <- tryCatch({
      extract_instruments(outcomes = met_id, p1 = 1e-5, r2 = 0.001, kb = 10000,
                          clump = TRUE, access_token = token)
    }, error = function(e) NULL)
    
    if(!is.null(met_inst) && nrow(met_inst) > 2) {
      cat("  Metabolite instruments:", nrow(met_inst), "SNPs\n")
      
      # Save instruments
      saveRDS(met_inst, paste0("met_", met_id, "_instruments.rds"))
      
      # Step 2: Metabolite → Jaundice (using MoBa as outcome)
      harm_met <- tryCatch({
        harmonise_data(exposure_dat = met_inst, outcome_dat = exposure_dat)
      }, error = function(e) NULL)
      
      if(!is.null(harm_met)) {
        harm_met <- harm_met[harm_met$mr_keep == TRUE, ]
        cat("  Harmonised:", nrow(harm_met), "SNPs\n")
        
        if(nrow(harm_met) >= 2) {
          mr_met <- mr(harm_met, method_list = "mr_ivw")
          cat(sprintf("  MR: OR=%.3f, p=%.2e\n", exp(mr_met$b[1]), mr_met$pval[1]))
        }
      }
    } else {
      cat("  No instruments found\n")
    }
  }
} else {
  cat("No metabolite GWAS found.\n")
}

# ========================================
# PART C: FUNCTIONAL PATHWAY MR
# ========================================
cat("\n\n========================================\n")
cat("PART C: FUNCTIONAL PATHWAY MR\n")
cat("========================================\n\n")

# Search for functional pathway GWAS
cat("Searching for functional pathway GWAS...\n")
pathway_search <- tryCatch({
  ieugwasr::search("functional pathway OR KEGG OR CAZy OR carbohydrate", limit=20)
}, error=function(e) NULL)

if(!is.null(pathway_search) && nrow(pathway_search) > 0) {
  cat("Found", nrow(pathway_search), "pathway GWAS datasets\n")
  
  # Look for relevant pathways
  relevant_pw <- pathway_search[
    grepl("pathway|KEGG|metabolism|carbohydrate|enzyme", 
          pathway_search$trait, ignore.case=TRUE), ]
  
  if(nrow(relevant_pw) > 0) {
    cat("Relevant pathways:", nrow(relevant_pw), "\n")
    for(i in 1:min(5, nrow(relevant_pw))) {
      cat(sprintf("  %s: %s\n", relevant_pw$id[i], substr(relevant_pw$trait[i], 1, 60)))
    }
  }
} else {
  cat("No pathway GWAS found.\n")
}

# ========================================
# PART D: COMBINED FORWARD + REVERSE MR TABLE
# ========================================
cat("\n\n========================================\n")
cat("PART D: COMPREHENSIVE MR SUMMARY TABLE\n")
cat("========================================\n\n")

# Load forward MR results
forward_res <- read.csv("mibio_jaundice_mr_results.csv")
forward_ivw <- forward_res[forward_res$method == "Inverse variance weighted", ]

# Load reverse MR results
if(file.exists("reverse_mr_results.csv")) {
  reverse_res <- read.csv("reverse_mr_results.csv")
  reverse_ivw <- reverse_res[reverse_res$method == "Inverse variance weighted", ]
  
  # Create combined table
  combined <- data.frame(
    Taxon = unique(c(forward_ivw$taxon, reverse_ivw$taxon)),
    Forward_OR = NA, Forward_P = NA, Forward_NSNP = NA,
    Reverse_OR = NA, Reverse_P = NA, Reverse_NSNP = NA
  )
  
  for(taxon in combined$Taxon) {
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
cat("ALL ANALYSES COMPLETE\n")
cat("========================================\n")
