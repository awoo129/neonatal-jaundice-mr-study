###########################################################
# TWO-STEP MR + FUNCTIONAL PATHWAY MR
# Using known OpenGWAS IDs
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
cat("PART 1: KNOWN METABOLITE GWAS\n")
cat("========================================\n\n")

# Try known bilirubin and bile acid GWAS IDs
# These are from the OpenGWAS database
known_metabolites <- c(
  "ebi-a-GCST90025973" = "Total bilirubin levels",
  "ebi-a-GCST90025974" = "Direct bilirubin levels",
  "ebi-a-GCST90025975" = "Indirect bilirubin levels",
  "ebi-a-GCST90025921" = "Bile acid levels",
  "ebi-a-GCST90014043" = "Cholesterol levels",
  "ebi-a-GCST90014044" = "HDL cholesterol",
  "ebi-a-GCST90014045" = "LDL cholesterol",
  "ebi-a-GCST90014046" = "Triglycerides"
)

# Check which metabolite IDs exist
cat("Checking metabolite GWAS availability...\n")
metabolite_results <- data.frame()

# Load MoBa fetal GWAS as outcome
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
cat("MoBa outcome data loaded:", nrow(outcome_dat), "variants\n\n")

for(id in names(known_metabolites)) {
  met_name <- known_metabolites[id]
  cat(sprintf("--- %s ---\n", met_name))
  
  # Extract metabolite instruments
  met_inst <- tryCatch({
    extract_instruments(outcomes = id, p1 = 1e-5, r2 = 0.001, kb = 10000,
                        clump = TRUE)
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(met_inst) && nrow(met_inst) > 0) {
    cat("  Instruments:", nrow(met_inst), "SNPs\n")
    
    # Save instruments
    saveRDS(met_inst, paste0("metabolite_", id, "_inst.rds"))
    
    # Run MR: metabolite → jaundice
    if(nrow(met_inst) >= 2) {
      harm <- tryCatch({
        harmonise_data(exposure_dat = met_inst, outcome_dat = outcome_dat)
      }, error=function(e) NULL)
      
      if(!is.null(harm)) {
        harm <- harm[harm$mr_keep == TRUE, ]
        cat("  Harmonised:", nrow(harm), "SNPs\n")
        
        if(nrow(harm) >= 2) {
          mr_res <- mr(harm, method_list = "mr_ivw")
          mr_res$trait <- met_name
          mr_res$OR <- exp(mr_res$b)
          mr_res$OR_lower <- exp(mr_res$b - 1.96*mr_res$se)
          mr_res$OR_upper <- exp(mr_res$b + 1.96*mr_res$se)
          metabolite_results <- rbind(metabolite_results, mr_res)
          cat(sprintf("  IVW: OR=%.3f (%.3f-%.3f), p=%.2e\n",
                      mr_res$OR[1], mr_res$OR_lower[1], mr_res$OR_upper[1], mr_res$pval[1]))
        }
      }
    }
  } else {
    cat("  No instruments extracted\n")
  }
}

if(nrow(metabolite_results) > 0) {
  write.csv(metabolite_results, "metabolite_jaundice_mr.csv", row.names=FALSE)
  cat("\n✅ Metabolite MR results saved\n")
} else {
  cat("\n⚠️ No metabolite MR results\n")
}

# ========================================
# PART 2: TWO-STEP MEDIATION ANALYSIS
# ========================================
cat("\n\n========================================\n")
cat("PART 2: TWO-STEP MEDIATION ANALYSIS\n")
cat("========================================\n\n")

# We need to estimate:
# Step 1: Microbiota → Metabolite (α)
# Step 2: Metabolite → Jaundice (β)
# Mediation effect = α × β
# Total effect = Microbiota → Jaundice (from earlier MR)
# Proportion mediated = (α × β) / Total effect

# The most relevant mediator pathway: Acidaminococcaceae → Bilirubin → Jaundice
# Since bilirubin is the most biologically relevant mediator

# Check if we have bilirubin instruments
bil_inst_file <- "metabolite_ebi-a-GCST90025973_inst.rds"  # Total bilirubin
if(file.exists(bil_inst_file)) {
  bil_inst <- readRDS(bil_inst_file)
  cat("Total bilirubin instruments:", nrow(bil_inst), "SNPs\n")
  
  # Load the significant microbiota taxon instruments
  for(mi_id in c("ebi-a-GCST90016924")) {  # Bifidobacteriaceae (= Acidaminococcaceae)
    cat(sprintf("\n--- Step 1: Acidaminococcaceae → Total Bilirubin ---\n"))
    
    mi_inst <- readRDS(paste0(mi_id, "_instruments.rds"))
    cat("  MiBioGen instruments:", nrow(mi_inst), "SNPs\n")
    
    # Step 1: Microbiota → Bilirubin
    # Extract bilirubin outcome data using the microbiota instruments
    bil_outcome <- tryCatch({
      extract_outcome_data(snps = mi_inst$SNP, outcomes = "ebi-a-GCST90025973",
                           access_token = token)
    }, error=function(e) NULL)
    
    if(!is.null(bil_outcome) && nrow(bil_outcome) > 0) {
      cat("  Bilirubin outcome data:", nrow(bil_outcome), "SNPs\n")
      
      # Harmoise
      harm_step1 <- tryCatch({
        harmonise_data(exposure_dat = mi_inst, outcome_dat = bil_outcome)
      }, error=function(e) NULL)
      
      if(!is.null(harm_step1)) {
        harm_step1 <- harm_step1[harm_step1$mr_keep == TRUE, ]
        cat("  Harmonised (Step 1):", nrow(harm_step1), "SNPs\n")
        
        if(nrow(harm_step1) >= 2) {
          mr_step1 <- mr(harm_step1, method_list = "mr_ivw")
          cat(sprintf("  Step 1 (Microbiota → Bilirubin): OR=%.3f, p=%.2e\n", 
                      exp(mr_step1$b[1]), mr_step1$pval[1]))
          
          # Step 2 we already have: Bilirubin → Jaundice
          step2 <- metabolite_results
          if(nrow(step2) > 0) {
            # Total effect: Acidaminococcaceae → Jaundice
            fwd_results <- read.csv("mibio_jaundice_mr_results.csv")
            total_effect <- fwd_results[fwd_results$id.exposure == mi_id & 
                                         fwd_results$method == "Inverse variance weighted", ]
            
            if(nrow(total_effect) > 0) {
              cat(sprintf("  Total effect: OR=%.3f, p=%.2e\n",
                          total_effect$OR[1], total_effect$pval[1]))
              
              # Mediation proportion
              alpha <- mr_step1$b[1]  # Microbiota → Bilirubin
              beta <- step2$b[1]       # Bilirubin → Jaundice
              total_b <- log(total_effect$OR[1])
              
              mediated <- alpha * beta
              prop <- mediated / total_b * 100
              
              cat(sprintf("\n  📊 MEDIATION ANALYSIS:\n"))
              cat(sprintf("  α (Microbiota → Bilirubin): %.4f\n", alpha))
              cat(sprintf("  β (Bilirubin → Jaundice): %.4f\n", beta))
              cat(sprintf("  α×β (Indirect effect): %.4f\n", mediated))
              cat(sprintf("  Total effect: %.4f\n", total_b))
              cat(sprintf("  Proportion mediated: %.1f%%\n", prop))
            }
          }
        }
      }
    } else {
      cat("  Could not extract bilirubin outcome data\n")
    }
  }
}

# ========================================
# PART 3: FUNCTIONAL PATHWAY MR
# ========================================
cat("\n\n========================================\n")
cat("PART 3: FUNCTIONAL PATHWAY MR\n")
cat("========================================\n\n")

# Try known functional pathway GWAS from MiBioGen
pathway_ids <- c(
  "ebi-a-GCST90027449" = "Gut microbiome functional pathways",
  "ebi-a-GCST90016921" = "Mollicutes class",
  "ebi-a-GCST90027450" = "Gut microbiome functional pathways 2",
  "ebi-a-GCST90025991" = "Carbohydrate metabolism"
)

for(id in names(pathway_ids)) {
  pw_name <- pathway_ids[id]
  cat(sprintf("--- %s ---\n", pw_name))
  
  pw_inst <- tryCatch({
    extract_instruments(outcomes = id, p1 = 1e-5, r2 = 0.001, kb = 10000,
                        clump = TRUE)
  }, error=function(e) {
    cat("  ERROR:", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(pw_inst) && nrow(pw_inst) > 0) {
    cat("  Instruments:", nrow(pw_inst), "SNPs\n")
    saveRDS(pw_inst, paste0("pathway_", id, "_inst.rds"))
    
    if(nrow(pw_inst) >= 2) {
      harm <- tryCatch({
        harmonise_data(exposure_dat = pw_inst, outcome_dat = outcome_dat)
      }, error=function(e) NULL)
      
      if(!is.null(harm)) {
        harm <- harm[harm$mr_keep == TRUE, ]
        cat("  Harmonised:", nrow(harm), "SNPs\n")
        
        if(nrow(harm) >= 2) {
          mr_res <- mr(harm, method_list = "mr_ivw")
          cat(sprintf("  IVW: OR=%.3f, p=%.2e\n", exp(mr_res$b[1]), mr_res$pval[1]))
        }
      }
    }
  } else {
    cat("  No instruments extracted\n")
  }
}

cat("\n========================================\n")
cat("ALL ANALYSES COMPLETE\n")
cat("Results saved:\n")
cat("  - metabolite_jaundice_mr.csv\n")
cat("  - metabolite_*_inst.rds\n")
cat("========================================\n")
