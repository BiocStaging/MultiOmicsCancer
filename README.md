# MultiOmicsCancer

An R package for integrative multi-omics analysis in cancer research, with a novel Immune Methylation Score (IMS) framework and built-in validation tools.

## Key Features

### Standard analyses
- Differential methylation (limma)
- Differential expression (DESeq2)
- Multi-omics integration
- Survival analysis (KM + Cox)
- Visualization

### Novel Q1-level tools
- `build_ims()` — Construct an Immune Methylation Score
- `bootstrap_validate_ims()` — Nested bootstrap validation
- `compare_ims_clinical()` — Compare IMS with clinical variables
- `estimate_purity()` — ESTIMATE tumor purity and immune score
- `validate_proteins()` — Protein-level validation via UniProt

## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Maryam-Emami/MultiOmicsCancer")
