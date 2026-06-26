# ======================================================
# Gut Microbiota and Neonatal Jaundice MR Study
# Master Run Script — Reproducibility Package
# ======================================================
# Before running:
#   1. Install R packages (see README.md)
#   2. Set OPENGWAS_JWT environment variable
#   3. Download MoBa GWAS files to ../data/
#   4. Set working directory to this script's location
# ======================================================

cat(rep("=", 60), sep="")
cat("\nReproducing: Gut Microbiota → Neonatal Jaundice MR Study\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("R version:", R.version.string, "\n")
cat(rep("=", 60), sep="\n")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # RStudio
# Or manually: setwd("path/to/reproducibility")

data_dir <- "data"
script_dir <- "scripts"
output_dir <- "output"
dir.create(output_dir, showWarnings=FALSE)

# Save session info
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
cat("\n[0/9] Session info saved\n")

# Step 1: Bidirectional MR (Forward + Reverse + Family-Trio)
cat("\n[1/9] Running bidirectional MR...\n")
source(file.path(script_dir, "bidirectional_mr.R"), echo=FALSE)
cat("  ✓ Done\n")

# Step 2: Deep analysis (MR-PRESSO, LOO, Steiger, Coloc)
cat("\n[2/9] Running sensitivity analyses...\n")
source(file.path(script_dir, "deep_analysis_v3.R"), echo=FALSE)
cat("  ✓ Done\n")

# Step 3: Microbiota taxa MR
cat("\n[3/9] Running 14 taxa MR analysis...\n")
source(file.path(script_dir, "mibio_analysis.R"), echo=FALSE)
cat("  ✓ Done\n")

# Step 4: MVMR Mediation
cat("\n[4/9] Running MVMR mediation...\n")
source(file.path(script_dir, "mediation_and_mvmr.R"), echo=FALSE)
cat("  ✓ Done\n")

# Step 5: Reverse MR
cat("\n[5/9] Running reverse MR...\n")
source(file.path(script_dir, "reverse_mr_correct.R"), echo=FALSE)
cat("  ✓ Done\n")

# Step 6: Colocalization
cat("\n[6/9] Running UGT1A1 colocalization...\n")
source(file.path(script_dir, "run_coloc_v2.R"), echo=FALSE)
cat("  ✓ Done\n")

# Step 7: Two-step MR / pathway
cat("\n[7/9] Running two-step MR and pathway analysis...\n")
source(file.path(script_dir, "two_step_and_pathway_mr.R"), echo=FALSE)
cat("  ✓ Done\n")

# Step 8: Supplementary tables
cat("\n[8/9] Generating supplementary tables...\n")
source(file.path(script_dir, "fill_gaps_fixed.R"), echo=FALSE)
cat("  ✓ Done\n")

# Step 9: Figures
cat("\n[9/9] Generating figures...\n")
source(file.path(script_dir, "generate_figures.R"), echo=FALSE)
cat("  ✓ Done\n")

cat(rep("=", 60), sep="")
cat("\n✅ Reproducibility run complete.\n")
cat("Results saved to:", output_dir, "\n")
cat(rep("=", 60), sep="")
