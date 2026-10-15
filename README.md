# Proteomic Profiling of *S. cerevisiae* Perturbed in RNA Polymerase III

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Interactive Jupyter notebook and visualization pipeline accompanying the Master's thesis:
> *"Proteomic profiling of Saccharomyces cerevisiae strains with perturbations in the activity of RNA Polymerase III"*  
> **Author:** Ignacy Makowski  
> **Affiliation:** Faculty of Mathematics, Informatics and Mechanics, University of Warsaw  

---

## Overview

This repository provides the code and visualization suite reproducing all figures presented in the thesis:
1. **Phase I (Mutant Profiling):** High-throughput SWATH-MS proteomic profiling of *maf1Δ* (constitutive RNAP III derepression) and *rpc128-1007* (assembly defect and phenotypic suppressor) grown across fermentative (YPD Glucose at 30°C) and respiratory conditions (YPGlycerol at 30°C and 37°C).
2. **Phase II (Multi-Dataset Harmonization):** Cross-study integration with standard laboratory strains (W303 and BY4742 in exponential and stationary phases) using empirical Bayes **ComBat** batch correction to eliminate technical inter-study artifacts while preserving strain-specific variance.

All plots adhere strictly to the color-blind friendly **Okabe-Ito palette** and publication-ready typography standards ($\ge 7.5\text{ pt}$).

---

## Thesis Figures

The primary interactive notebook `Minimalist_Plot_Generator.ipynb` generates vector visualizations for every figure:

| Figure | Scope | Description |
|---|---|---|
| **Figure 4.1** | SWATH-MS Dataset | Hierarchical clustering heatmap and PCA projections (by condition and strain). |
| **Figure 4.2** | Technical Reproducibility | Coefficient of Variation (CV) boxplot distributions across replicates. |
| **Figure 4.3** | Differential Abundance (*limma*) | 5-panel volcano plots comparing mutants vs WT with moderated t-statistics. |
| **Figure 4.4** | Codon Usage Bias | Scatter plots of protein abundance Z-scores vs target codon frequencies (`AAA`, `CAA`, `GAA`, `AAG`). |
| **Figure 4.5** | Codon Usage Bias | Categorical heatmap comparing significant regulation against codon usage thresholds. |
| **Figures 4.6 – 4.10** | Functional Enrichment | Top enriched KEGG metabolic pathways (10% overlap threshold) across carbon sources. |
| **Figure 4.17** | Stress Response | Targeted abundance profiles for molecular chaperones. |
| **Figure 4.21** | Batch Effect Correction | Evaluation of ComBat adjustment: PCA and intensity boxplots before vs after correction. |
| **Figure 4.22** | Harmonized Proteome | Hierarchical clustered heatmap of the merged, batch-corrected dataset (1,589 common proteins). |
| **Figure 4.23** | Cross-Strain Architecture | High-dimensional UMAP embedding across all genotypes, carbon sources, and phases. |

---

## Quick Start

### 1. Environment Setup
```bash
git clone https://github.com/makowskiignacy/rnap3-proteomics-pipeline.git
cd rnap3-proteomics-pipeline
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### 2. Launch Notebook
```bash
jupyter notebook Minimalist_Plot_Generator.ipynb
```
*Note: All figures are already pre-rendered in the notebook for immediate inspection without running any code.*

---

## Data Availability

The raw and processed proteomic datasets used in this study are available from the author upon reasonable request. Once obtained, data tables can be placed in the `data/` directory to re-execute the statistical pipeline from scratch.

---

## Citation

```bibtex
@mastersthesis{makowski2026rnap3,
  author  = {Ignacy Makowski},
  title   = {Proteomic profiling of Saccharomyces cerevisiae strains with perturbations in the activity of RNA Polymerase III},
  school  = {University of Warsaw},
  year    = {2026}
}
```
