###########################################################
# SUMMARY FIGURES: Multi-panel results visualization
###########################################################

library(ggplot2)
library(data.table)

setwd("C:/tmp/gwas_data")

# Load all results
forward_res <- read.csv("mibio_jaundice_mr_results.csv")
reverse_res <- read.csv("reverse_mr_results.csv")
metabolite_res <- read.csv("metabolite_jaundice_mr.csv")

# FIGURE 1: Forest plot - Forward MR all taxa
forward_ivw <- forward_res[forward_res$method == "Inverse variance weighted", ]
forward_ivw$OR <- as.numeric(forward_ivw$OR)
forward_ivw$OR_lower <- as.numeric(forward_ivw$OR_lower)
forward_ivw$OR_upper <- as.numeric(forward_ivw$OR_upper)
forward_ivw$pval <- as.numeric(forward_ivw$pval)

forward_ivw$label <- forward_ivw$taxon
forward_ivw$sig <- ifelse(forward_ivw$bonf_sig == "TRUE", "Bonferroni sig", 
                          ifelse(forward_ivw$pval < 0.05, "Nominal sig", "NS"))

# Sort by OR
forward_ivw <- forward_ivw[order(forward_ivw$OR), ]
forward_ivw$label <- factor(forward_ivw$label, levels = forward_ivw$label)

p1 <- ggplot(forward_ivw, aes(x = OR, y = label, color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper), height = 0.2) +
  scale_x_log10(breaks = c(0.3, 0.5, 1, 2, 3, 5)) +
  scale_color_manual(values = c("Bonferroni sig" = "red", 
                                 "Nominal sig" = "orange",
                                 "NS" = "gray50")) +
  labs(title = "A. Forward MR: Gut Microbiota → Neonatal Jaundice",
       x = "OR (95% CI)", y = "", color = "Significance") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")
ggsave("figure_forward_forest.png", p1, width = 8, height = 6, dpi = 300)
cat("✅ Figure 1: Forward forest plot saved\n")

# FIGURE 2: Forest plot - Reverse MR
if(nrow(reverse_res) > 0) {
  reverse_ivw <- reverse_res[reverse_res$method == "Inverse variance weighted", ]
  reverse_ivw$OR <- as.numeric(reverse_ivw$OR)
  reverse_ivw$OR_lower <- as.numeric(reverse_ivw$OR_lower)
  reverse_ivw$OR_upper <- as.numeric(reverse_ivw$OR_upper)
  reverse_ivw$pval <- as.numeric(reverse_ivw$pval)
  
  # Cap extreme CIs
  reverse_ivw$OR_lower[reverse_ivw$OR_lower < 0.1] <- 0.1
  reverse_ivw$OR_upper[reverse_ivw$OR_upper > 10] <- 10
  
  reverse_ivw$label <- reverse_ivw$taxon
  reverse_ivw$sig <- ifelse(reverse_ivw$pval < 0.05, "p<0.05", "NS")
  
  reverse_ivw <- reverse_ivw[order(reverse_ivw$OR), ]
  reverse_ivw$label <- factor(reverse_ivw$label, levels = reverse_ivw$label)
  
  p2 <- ggplot(reverse_ivw, aes(x = OR, y = label, color = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper), height = 0.2) +
    scale_x_log10() +
    scale_color_manual(values = c("p<0.05" = "orange", "NS" = "gray50")) +
    labs(title = "B. Reverse MR: Jaundice → Gut Microbiota",
         x = "OR (95% CI)", y = "", color = "Significance") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave("figure_reverse_forest.png", p2, width = 8, height = 6, dpi = 300)
  cat("✅ Figure 2: Reverse forest plot saved\n")
}

# FIGURE 3: Metabolite forest
if(nrow(metabolite_res) > 0) {
  met_ivw <- metabolite_res[metabolite_res$method == "Inverse variance weighted", ]
  if(nrow(met_ivw) > 0) {
    met_ivw$OR <- as.numeric(met_ivw$OR)
    met_ivw$OR_lower <- as.numeric(met_ivw$OR_lower)
    met_ivw$OR_upper <- as.numeric(met_ivw$OR_upper)
    met_ivw$pval <- as.numeric(met_ivw$pval)
    met_ivw$sig <- ifelse(met_ivw$pval < 0.05, "p<0.05", "NS")
    
    met_ivw$label <- met_ivw$trait
    met_ivw <- met_ivw[order(met_ivw$OR), ]
    met_ivw$label <- factor(met_ivw$label, levels = met_ivw$label)
    
    p3 <- ggplot(met_ivw, aes(x = OR, y = label, color = sig)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      geom_point(size = 4) +
      geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper), height = 0.2) +
      scale_x_log10() +
      scale_color_manual(values = c("p<0.05" = "red", "NS" = "gray50")) +
      labs(title = "C. Metabolite → Neonatal Jaundice MR",
           x = "OR (95% CI)", y = "", color = "Significance") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "bottom")
    ggsave("figure_metabolite_forest.png", p3, width = 8, height = 4, dpi = 300)
    cat("✅ Figure 3: Metabolite forest plot saved\n")
  }
}

# FIGURE 4: Combined Multi-panel summary figure
cat("\nCreating summary visualization...\n")

# Create manual LOO plot
loo_res <- readRDS("sensitivity_acidaminococcaceae_loo.rds")
if(inherits(loo_res, "data.frame") && nrow(loo_res) > 0) {
  names(loo_res) <- make.names(names(loo_res), allow_=TRUE)
  if("b" %in% names(loo_res) && "se" %in% names(loo_res) && "SNP" %in% names(loo_res)) {
  # Get IVW result without each SNP
  p4 <- ggplot(loo_res, aes(x = b, y = SNP)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 2.5, color = "#2166AC") +
    geom_errorbarh(aes(xmin = b - 1.96*se, xmax = b + 1.96*se), 
                   height = 0.2, color = "#2166AC") +
    labs(title = "Leave-One-Out: Acidaminococcaceae → Jaundice",
         x = "MR effect size (log-OR)", y = "SNP removed") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))
  ggsave("figure_loo_acidaminococcaceae.png", p4, width = 8, height = 6, dpi = 300)
  cat("✅ Figure 4: LOO plot saved\n")
  }
}

# FIGURE 5: Mediation diagram
cat("\nGenerating mediation diagram...\n")

mediation_data <- data.frame(
  path = c("Total Effect", 
           "Direct Effect", 
           "Indirect via Bilirubin"),
  estimate = c(0.7738, 0.7807, -0.0069),
  prop = c(100, 100.9, -0.9)
)

p5 <- ggplot(mediation_data, aes(x = path, y = estimate, fill = path)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", prop)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Total Effect" = "#2166AC", 
                                "Direct Effect" = "#35978F", 
                                "Indirect via Bilirubin" = "#BF812D")) +
  labs(title = "D. Mediation Analysis: Acidaminococcaceae → Jaundice",
       x = "", y = "Effect Size (log-OR)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none")
ggsave("figure_mediation.png", p5, width = 6, height = 4, dpi = 300)
cat("✅ Figure 5: Mediation diagram saved\n")

cat("\n\n========================================\n")
cat("ALL FIGURES GENERATED:\n")
cat("========================================\n")
cat("1. figure_forward_forest.png - Forward MR results\n")
cat("2. figure_reverse_forest.png - Reverse MR results\n")
cat("3. figure_metabolite_forest.png - Metabolite MR\n")
cat("4. figure_loo_acidaminococcaceae.png - LOO analysis\n")
cat("5. figure_mediation.png - Mediation analysis\n")
cat("========================================\n")
