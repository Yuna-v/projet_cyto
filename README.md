# projet_cyto

# Myeloid Flow Cytometry Analysis Pipeline

A complete R-based analysis pipeline for high-dimensional mass cytometry (CyTOF) data, comparing immune cell responses between **Post-Prime** and **Post-Boost** vaccination conditions over time.

---

## Overview

This pipeline processes FCS files through a full analysis workflow:

- Arcsinh transformation and data normalization
- Dimensionality reduction (UMAP)
- Unsupervised clustering (k-means, k=50)
- Metaclustering and phenotypic annotation
- Differential abundance testing (Wilcoxon, FDR-corrected)
- Visualization (PCA, volcano plots, kinetics, heatmaps)

---

## Project Structure

```
projet_cyto/
├── Myeloid_data/          # FCS files (not included, see Data section)
├── analysis.Rmd           # Main analysis R Markdown file
├── README.md
└── .gitignore
```

---

## Data

This pipeline was developed to reanalyze mass cytometry data from the following publication:

> Palgen, J.-L., Tchitchek, N., Elhmouzi-Younes, J., Delandre, S., Namet, I., Rosenbaum, P., Dereuddre-Bosquet, N., Martinon, F., Cosma, A., Lévy, Y., et al. (2018). **Prime and Boost Vaccination Elicit a Distinct Innate Myeloid Cell Immune Response.** *Scientific Reports* 8, 3087. https://doi.org/10.1038/s41598-018-21222-2

The study examines how prime and boost vaccinations differentially activate myeloid innate immune cells, with samples collected at multiple timepoints (Day 0, 1, 3, and 14) post-vaccination.

---

FCS files are **not included** in this repository due to file size constraints.

Place your `.fcs` files in a folder named `Myeloid_data/` at the root of the project before running the analysis.

Expected file naming convention:
```
PPD00H00_indiv1.fcs   # Post-Prime, Day 0
PPD01H00_indiv1.fcs   # Post-Prime, Day 1
PBD03H00_indiv1.fcs   # Post-Boost, Day 3
PBD14H00_indiv1.fcs   # Post-Boost, Day 14
```

| Prefix | Condition   | Description                        |
|--------|-------------|------------------------------------|
| `PPD`  | Post-Prime  | Samples collected after first vaccination  |
| `PBD`  | Post-Boost  | Samples collected after booster vaccination |

Timepoints available: **Day 0, 1, 3, 14** post-vaccination.
---

## Requirements

### R version
R >= 4.1.0

### R packages

```r
install.packages(c(
  "dplyr", "tidyr", "tibble", "ggplot2",
  "gridExtra", "ggrepel", "pheatmap",
  "cluster", "mclust", "lme4", "lmerTest",
  "here", "tidytext", "ggpubr"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# GitHub packages
devtools::install_github("tchitchek-lab/CellVizR")
```

---

## How to Run

1. Clone the repository:
```bash
git clone https://github.com/votre-username/projet_cyto.git
cd projet_cyto
```

2. Place your FCS files in `Myeloid_data/`

3. Open `analysis.Rmd` in RStudio

4. Knit the document or run chunks sequentially

> ⚠️ The pipeline uses `here::here()` for relative paths — make sure to open the project via the `.Rproj` file so the working directory is set correctly.

---

## Pipeline Steps

| Step | Description |
|------|-------------|
| 1 | Library loading and seed setting |
| 2 | FCS file import |
| 3 | Arcsinh transformation |
| 4 | Data consolidation into a single data frame |
| 5 | Marker definition (clustering vs cytokine panels) |
| 6 | CellVizR import with marker exclusion |
| 7 | Metadata preparation (condition, timepoint, individual) |
| 8 | UMAP manifold generation |
| 9 | K-means clustering (k=50) |
| 10 | PCA on cluster abundances |
| 11 | Global population kinetics |
| 12 | Hierarchical metaclustering (6 metaclusters) |
| 13 | Differential expression testing (Wilcoxon + BH correction) |
| 14 | Volcano plot (abundance Boost vs Prime) |
| 15 | APC-focused analysis (top HLADR⁺ CD86⁺ clusters) |
| 16 | Metacluster 1 deep-dive (abundance, MFI, kinetics) |

---

## Key Outputs

- **PCA plot** — sample-level overview colored by condition and timepoint
- **UMAP** — single-cell manifold with cluster overlays
- **Heatmap** — marker expression profiles per cluster
- **Volcano plot** — differential cluster abundance (Boost vs Prime, FDR < 0.05)
- **Kinetics plots** — temporal dynamics of top clusters and metaclusters
- **Boxplots** — abundance and MFI comparisons per timepoint
- **Individual trajectories** — spaghetti plots per donor

---

## Statistical Methods

- **Differential abundance**: Wilcoxon rank-sum test, FDR correction (Benjamini-Hochberg)
- **Differential expression**: Wilcoxon rank-sum test per cluster × marker combination
- **Metaclustering**: Hierarchical clustering (Ward's method) on median marker expression profiles

---

## Authors

Fanny Xie & Guillaume Masson
