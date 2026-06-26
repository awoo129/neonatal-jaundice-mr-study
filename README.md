# Gut Microbiota and Neonatal Jaundice — MR Study
## Reproducibility Package

### Manuscript
**Title**: Gut Microbiota and Neonatal Jaundice: A Bidirectional Mendelian Randomization Study with Family-Trio Design
**Authors**: [Author Names]
**Journal**: Frontiers in Genetics (submitted)

### Contents

```
reproducibility/
├── README.md              ← This file
├── scripts/               ← R analysis scripts (numbered in order)
├── data/                  ← Instrument files + intermediate results
└── output/                ← Figures + result tables (generated)
```

### Data Sources (public, requires download)

| Data | Source | Notes |
|:-----|:-------|:------|
| MoBa Neonatal Jaundice GWAS (Fetal) | https://www.fhi.no/en/ch/studies/moba/for-forskere-artikler/gwas-data-from-moba/ | `moba-gwas-jaundice-fets.txt.gz` (~266 MB) |
| MoBa Neonatal Jaundice GWAS (Maternal) | Same URL | `moba-gwas-jaundice-moms.txt.gz` |
| MoBa Neonatal Jaundice GWAS (Paternal) | Same URL | `moba-gwas-jaundice-dads.txt.gz` |
| Gut Microbiota GWAS (MiBioGen) | IEU OpenGWAS (ebi-a-GCST90016921–90016935) | Extracted via ieugwasr API |
| Bilirubin GWAS | IEU OpenGWAS (ebi-a-GCST90025973) | N=436,748 |

### Required R / Python Environment

**R 4.5.3** + packages:

| Package | Version | Purpose |
|:--------|:--------|:--------|
| TwoSampleMR | 0.7.5 | Core MR engine |
| MRPRESSO | 1.0 | Pleiotropy outlier detection |
| coloc | 5.2.3 | Colocalization analysis |
| MRMix | 0.1.0 | Robust mixture-model MR |
| RadialMR | 1.2.3 | Radial plot IVW |
| ieugwasr | 1.1.0 | OpenGWAS API client |
| ggplot2 | 4.0.3 | Figures |
| data.table | 1.18.4 | Fast data I/O |

Install in R:
```r
install.packages(c("TwoSampleMR", "MRPRESSO", "coloc", "RadialMR", "MRMix",
                   "ggplot2", "data.table", "ieugwasr", "stringr", "httr", "jsonlite"))
```

### Setup

1. **Register for OpenGWAS**: Get a JWT token at https://api.opengwas.io/
2. **Set environment variable**: `export OPENGWAS_JWT="your_token_here"`
3. **Download MoBa GWAS files** from FHI website → place in `data/`
4. **Set working directory** in each script to `reproducibility/`

### Analysis Pipeline (run in order)

```
1. bidirectional_mr.R           ← Main MR (forward + reverse + trio)
2. deep_analysis_v3.R           ← Sensitivity analyses (MR-PRESSO, LOO, Steiger, Coloc)
3. mibio_analysis.R             ← Microbiota taxa MR (14 taxa)
4. mediation_and_mvmr.R         ← MVMR mediation decomposition
5. reverse_mr_correct.R         ← Reverse MR details
6. run_coloc_v2.R               ← UGT1A1 colocalization
7. two_step_and_pathway_mr.R    ← Supplementary pathway analysis
8. fill_gaps_fixed.R            ← F-stats, PhenoScanner, supp tables
9. generate_figures.R           ← All 7 figures
```

### Key Results (pre-computed, in `data/`)

| File | Contents |
|:-----|:---------|
| `mibio_jaundice_mr_results.csv` | 14 taxa × 3 MR methods × 3 cohorts |
| `reverse_mr_results.csv` | Jaundice → microbiota reverse MR |
| `colocalization_summary.csv` | UGT1A1 colocalization (PP.H3=0.939) |
| `mr_presso_fetal.rds` | MR-PRESSO raw results |
| `mr_results_summary.csv` | Summary of primary findings |

### Contact

Corresponding author: [Name, Email]

---

*Reproducibility package generated June 26, 2026*
