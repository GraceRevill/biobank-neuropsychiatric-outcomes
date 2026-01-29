# Traumatic Brain Injury Increases Risk of Psychiatric Disorders Following Neurological Disease: Evidence from the UK Biobank Cohort

R code for analysis of UK Biobank data investigating impact of TBI on psychiatric outcomes in people with neurological conditions

This archive contains the R code for the analysis reported in the above study.

This repository contains:

1. [Revill_et_al_preprocessing.R](https://github.com/GraceRevill/biobank-neuropsychiatric-outcomes/blob/main/Revill_et_al_preprocessing.R) - R script to complete data cleaning and preprocessing 
2. [Revill_et_al_tbi_neuro_mh_analysis.R](https://github.com/GraceRevill/biobank-neuropsychiatric-outcomes/blob/main/Revill_et_al_tbi_neuro_mh_analysis.R) - R script to complete mixed effects models
3. [Revill_et_al_KM_final.R](https://github.com/GraceRevill/biobank-neuropsychiatric-outcomes/blob/main/Revill_et_al_KM_final.R) - R script to complete survival modelling (Kaplan-Meier and Cox analyses) 

# Dataset
UK Biobank (UKB) is a large, prospective cohort of over 500 000 generally healthy adults aged 40–70 years at recruitment (2006–2010), registered with the UK National Health Service and assessed at one of 22 recruitment centres https://www.ukbiobank.ac.uk/. 
Approval for the current study was granted by the University College London Ethics Committee 

# Platform and package versions
R language version, and package versions used to generate the results are:<br>

R Version 4.5.1<br>

Package version for readr is 2.1.4<br>
Package version for dplyr is 1.1.2<br>
Package version for tidyr is 1.3.0<br>
Package version for ggplot2 is 3.5.2<br>
Package version for Hmisc is 4.7.2<br>
Package version for patchwork is 1.1.3<br>
Package version for lme4 is 1.1.3<br>
Package version for stringr is 1.5.0<br>
Package version for reshape2 is 1.1.3<br>
Package version for gtsummary is 1.7.1<br>
Package version for gt is 0.9.0<br>
Package version for afex is 1.3.0<br>
Package version for lubridate is 1.9.2<br>
Package version for sjPlot is 2.8.14<br>
Package version for jtools is 2.2.1<br>
Package version for purrr is 1.0.1<br>
Package version for glmmTMB is 1.1.9<br>
Package version for survminer is 0.5.0<br>
Package version for gridExtra is 2.3.0<br>
