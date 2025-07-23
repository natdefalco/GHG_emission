# GHG_emission
 💨 Smart Chamber GHG Flux Pipeline

This repository contains a working pipeline for the **calculation of greenhouse gas (GHG) fluxes** 
— particularly **CH₄ and CO₂** — from static or floating chamber measurements in aquatic systems using data collected from **LI-8200 smart chambers**.


## 🚀 What This Project Does

### ✅ R Pipeline (with `aquaGHG` and `goFlux`)
- Parses LI-8200 JSON data
- Detects valid measurement windows
- Calculates CH₄ fluxes using **linear (LM)** and **non-linear (HM)** models
- Separates **diffusion** and **ebullition** fluxes for CH₄
- Exports results to `.csv`
- Generates diagnostic plots and saves them as `.pdf`

You can read more about the packages here:
- [`aquaGHG`](https://github.com/camilleminaudo/aquaGHG): Wrapper for CH₄ flux processing and ebullition separation
- [`goFlux`](https://qepanna.quarto.pub/goflux/): Core library for GHG flux regression and visualization

### 🐍 Python Scripts
- Allow custom CH₄/CO₂ flux calculation
- Provide an option to **validate flux results against Smart Chamber device outputs**
- Use regression and plotting via `scipy`, `sklearn`, `matplotlib`, and `seaborn`

---

## 🔧 How to Use

### ✅ R Instructions

1. Install required R packages:
```r
install.packages(c("devtools", "goFlux", "ggplot2", "dplyr", "pbapply", 
                   "openxlsx", "ggnewscale", "tidyr", "purrr", "readxl"))
devtools::install_github("camilleminaudo/aquaGHG")
