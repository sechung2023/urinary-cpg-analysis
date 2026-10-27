# urinary-cpg-analysis
R scripts for statistical analyses and figures examining urinary CPG as a biomarker of non-small cell lung cancer in never-smokers across two independent case–control studies.

# Urinary CPG Analysis
R scripts supporting the manuscript:

**Urinary 27-nor-5β-cholestane-3α,7α,12α,24R,25S-pentol glucuronide as a biomarker of non-small cell lung cancer in never-smokers: two independent case–control studies**

## Overview

This repository contains scripts for statistical analyses and figure generation using the NCI-Maryland exploratory cohort and Colorado validation cohort.

Analyses include:
- Comparisons of urinary CPG across case/control and smoking groups
- Receiver-operating-characteristic (ROC) curves
- Kaplan–Meier survival analyses
- Logistic and Cox regression analyses
- SHAP feature-attribution plots

## Script organization

Scripts are named according to their corresponding manuscript figures, supplementary figures, and supplementary tables.

- **Figure 1:** CPG comparisons by case/control and smoking status
- **Figure 2 and Supplementary Figure S3:** ROC analyses
- **Figure 3:** Kaplan–Meier analyses in never-smoker cases
- **Supplementary Figure S2:** CPG comparisons by detailed smoking category
- **Supplementary Figure S4:** Survival analyses using a combined CR + NANA + CPG score
- **Supplementary Figure S5:** Logistic regression forest plots
- **Supplementary Figure S6:** Cox regression forest plots
- **Supplementary Figure S7:** SHAP summary plots
- **Supplementary Table S1:** Logistic regression tables, with separate scripts for each cohort
- **Supplementary Table S2:** Cox regression tables, with separate scripts for each cohort

This repository covers the analyses listed above; it does not contain code for every analysis in the manuscript.

## Data requirements

The underlying datasets are not included.

Scripts expect the appropriate prepared data frames to be loaded in the R environment:
- `NIHMD`: NCI-Maryland exploratory cohort
- `COSPFINAL`: Colorado validation cohort

Variable names and categorical coding must match the expectations of each script. Coding conventions differ between prepared datasets and analyses; consult the comments and recoding steps within each script.

Depending on the script, results are displayed in R or saved as figure files or CSV tables. Final manuscript panel assembly and formatting may have been performed separately.
