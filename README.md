# Urinary Albumin-to-Creatinine Ratio and Vascular Risk in Type 2 Diabetes

## 📌 Project Status
**⚠️ Under Review** - This research is currently under peer review. Results and specific findings are not yet publicly available.

## 🔍 Overview
A population-based retrospective cohort study examining the prognostic value of baseline urinary albumin-to-creatinine ratio (uACR) for cardiovascular outcomes in newly diagnosed type 2 diabetes using linked primary and secondary care data from England (CPRD-GOLD, HES, ONS).

## 💡 Research Questions
- Does baseline uACR predict multi-domain vascular outcomes in incident diabetes?
- What is the shape of the exposure-response relationship?
- How do associations vary across different vascular disease domains?

## 🛠️ Methodology

### Study Design
- **Type:** Retrospective cohort study
- **Data Source:** CPRD-GOLD linked to HES and ONS mortality data
- **Population:** Adults with incident type 2 diabetes (pre-existing CVD/diabetes excluded)
- **Exposure:** Baseline uACR (categorical: KDIGO A1-A3; continuous: restricted cubic splines)

### Outcomes
**Primary Domains:**
- Atherosclerotic cardiovascular disease (ASCVD)
- Heart failure (HF)
- Peripheral artery disease (PAD)

**Secondary Outcomes:**
- All-cause mortality
- Atrial fibrillation
- Composite vascular endpoint

### Statistical Analysis
- **Competing Risks:** Cause-specific Cox models with g-computation
- **Absolute Risks:** Population-standardized risk estimates
- **Missing Data:** Multiple imputation by chained equations (MICE) with biomarker-specific predictor matrices
- **Sensitivity Analyses:** Fine-Gray subdistribution hazards models
- **Effect Modification:** Sex-stratified analyses and interaction testing
- **Model Performance:** Time-dependent AUC, Brier scores, calibration plots

## 📁 Repository Structure
- pipeline analysis : generate all analysis required
- extract data from output bundles, generate tables and plots
- Baseline table generators: to generate the baseline table and event table
- Phenotype: the READ Codes and ICD-10 code-lists used in this study 
- Results : Not publicly available yet