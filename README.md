# GHG_emission
 💨 Smart Chamber GHG Flux Pipeline

This repository contains a working pipeline for the calculation of greenhouse gas (GHG) fluxes
— particularly **CH₄** and **CO₂** — from static or floating chamber measurements in aquatic systems using data collected from **LI-8200 smart chambers**.

The project includes:
- An **R script** using the `aquaGHG` and `goFlux` packages for automatic CH₄ flux calculation and ebullition/diffusion separation
- A **Python script** for custom CH₄ and CO₂ flux regression and validation against Smart Chamber output

---

## 🚀 What This Project Does

### ✅ R Pipeline (with `aquaGHG` and `goFlux`) (`CH4_flux_aqua_GHG.R`)
- Parses LI-8200 JSON data
- Detects valid measurement windows
- Calculates CH₄ fluxes using **linear (LM)** and **non-linear (HM)** models
- Separates **diffusion** and **ebullition** fluxes for CH₄
- Exports results to `.csv`
- Generates diagnostic plots and saves them as `.pdf`

You can read more about the packages here:
- [`aquaGHG`](https://github.com/camilleminaudo/aquaGHG): Wrapper for CH₄ flux processing and ebullition separation
- [`goFlux`](https://qepanna.quarto.pub/goflux/): Core library for GHG flux regression and visualization

### 🐍 Python Scripts (`flux_calculation_CH4_CO2.py`)
- Computes CH₄ and CO₂ fluxes using linear regression
- Validates calculated fluxes against Smart Chamber device outputs (if available)
- Creates visualizations (raw signal + fitted regression)

## 📊 Output Files

The **R script** produces:
- `CH4_combined_flux_results.csv`: raw CH₄ flux per event
- `CH4_diffusion_vs_ebullition.csv`: CH₄ breakdown (diffusion, ebullition)
- `CH4_aquaGHG_flux_plots.pdf`: flux model diagnostics
- (optional) `CH4_combined_flux_with_ebullition.csv`: merged results with quality metrics

The **Python script** produces:
- `*_flux_results.csv`: calculated fluxes per gas and event
- `.png` or `.pdf` plots comparing model vs. device output

---

## 🔧 How to Use

### ✅ R Instructions

1. Install required R packages:
```r
install.packages(c("devtools", "goFlux", "ggplot2", "dplyr", "pbapply", 
                   "openxlsx", "ggnewscale", "tidyr", "purrr", "readxl"))
devtools::install_github("camilleminaudo/aquaGHG")


