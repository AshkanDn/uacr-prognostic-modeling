# uACR and Multi-Domain Vascular Risk in Incident Type 2 Diabetes

![R](https://img.shields.io/badge/R-4.3+-blue)
![Status](https://img.shields.io/badge/Status-Under%20Review-orange)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## 📌 Project Overview

A reproducible epidemiologic analysis pipeline evaluating urinary albumin-to-creatinine ratio (uACR) as a prognostic marker for multi-domain vascular outcomes in adults with newly diagnosed type 2 diabetes using linked electronic health records from England.

This repository includes advanced survival analysis workflows, competing-risk modeling, multiple imputation, and reproducible clinical research engineering using R.


**⚠️ Under Review** - This research is currently under peer review. Results and specific findings are not yet publicly available.

---

## 🔬 Scientific Background

Albuminuria is an established marker of renal dysfunction and cardiovascular risk. However, its prognostic role across multiple vascular disease domains in incident type 2 diabetes remains incompletely characterised.

This study investigates:
- Whether baseline uACR predicts future vascular outcomes
- The shape of the exposure-response relationship
- Differences across cardiovascular disease domains
- Absolute and relative risk trajectories over time

---

## 💡 Research Questions

- Does baseline uACR predict multi-domain vascular outcomes in incident type 2 diabetes?
- What is the shape of the exposure-response relationship?
- How do associations vary across different vascular disease domains?
- Does albuminuria improve vascular risk stratification?

---
# 🛠️ Methodology

## Study Design

- **Design:** Population-based retrospective cohort study
- **Data Source:** CPRD-GOLD linked to HES and ONS mortality records
- **Population:** Adults with incident type 2 diabetes
- **Exclusions:** Pre-existing diabetes and cardiovascular disease
- **Exposure:** Baseline urinary albumin-to-creatinine ratio (uACR)

### Exposure Modeling
- KDIGO categorical classification (A1–A3)
- Continuous modeling using restricted cubic splines

---

## Outcomes

### Primary Outcomes
- Atherosclerotic cardiovascular disease (ASCVD)
- Heart failure (HF)
- Peripheral artery disease (PAD)

### Secondary Outcomes
- All-cause mortality
- Atrial fibrillation
- Composite vascular outcomes

---

## Statistical Framework

### Absolute Risk Estimation
- G-computation
- Population-standardised risk estimation
- Risk differences and risk ratios

### Survival Modeling
- Cause-specific Cox proportional hazards models
- Competing-risk analysis framework

### Missing Data
- Multiple imputation by chained equations (MICE)
- Biomarker-specific predictor matrices

### Sensitivity Analyses
- Fine-Gray subdistribution hazards models
- Sex-stratified analyses
- Interaction testing

### Model Evaluation
- Time-dependent AUC
- Brier scores
- Calibration assessment
- Restricted cubic spline modeling

---

# ⚙️ Analysis Workflow

The analytical workflow is fully modular and reproducible.

1. Cohort assembly and eligibility filtering
2. Baseline variable harmonisation
3. Missing data imputation using MICE
4. Endpoint-specific competing-risk modeling
5. Absolute risk estimation via g-computation
6. Calibration and discrimination analysis
7. Automated table and figure generation
8. Sensitivity and subgroup analyses

---

---

# 📁 Repository Structure

```text
uacr-prognostic-modeling/
│
├── README.md
├── LICENSE
├── .gitignore
│
├── scripts/
│   ├── 01_baseline_tables.R
│   ├── 02_analysis_pipeline.R
│   ├── 03_results_pipeline.R
│   └── utils.R
│
├── phenotypes/
│   ├── READ2/ read_codes.csv
│   └── ICD10/ icd10_codes.csv
│   └── OPCS/ OPCS_codes.csv
│   └── BNF/ BNF_codes.csv
│
├── docs/
│   ├── consort_diagram.png
│   └── manuscript.docx
│
├── outputs/
│   ├── figures/
│   ├── tables/
│   └── logs/
│
└── manuscript/
│   └── manuscript.docx
│   └── supplementary_methods.docx
    
```

---

  

