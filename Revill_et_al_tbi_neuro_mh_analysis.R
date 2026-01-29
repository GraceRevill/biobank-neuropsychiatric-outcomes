################################################################

# Analysis of Biobank dataset

# Revill et al

################################################################

# Load libraries 
library(purrr)
library(dplyr)
library(tidyr)
library(lubridate)
library(gt)
library(glmmTMB)
library(ggplot2)
library(gtsummary)  

# Remove everything from memory if needed to start afresh
# rm(list = ls())


# Set data and output directories
data_dir   <- "S:/Track_TBI/Grace/Test/"
output_dir <- "S:/Track_TBI/Grace/Test/Output/"

setwd(data_dir)


################################################################
# 1. NEURO DATASET + TBI DEFINITIONS
################################################################

# Load dataset 
neuro    <- read.csv(file.path(output_dir, "neuro.csv"))


# Ensure all *date variables are Date 
date_cols <- grep("date", names(neuro), value = TRUE, ignore.case = TRUE)
neuro[date_cols] <- lapply(neuro[date_cols], as.Date)

################################################################
# GLOBAL EXCLUSION: remove ICD-coded TBI occurring AFTER MH (different year)
################################################################

neuro <- neuro %>%
  mutate(
    mh_year  = year(any_mh_date),
    icd_year = year(injury_date_B)
  ) %>%
  filter(
    # KEEP anyone who has no MH diagnosis
    is.na(any_mh_date) |
      # KEEP if they have no ICD TBI date
      is.na(injury_date_B) |
      # KEEP if ICD TBI occurs before MH
      injury_date_B <= any_mh_date |
      # KEEP if ICD TBI occurs in the same year as MH
      icd_year == mh_year
  )

# Check this worked: should be 0
sum(
  neuro$injury_date_B > neuro$any_mh_date &
    year(neuro$injury_date_B) != year(neuro$any_mh_date),
  na.rm = TRUE
)


# Construct final TBI variable 
neuro <- neuro %>%
  mutate(
    mh_year  = year(any_mh_date),
    icd_year = year(injury_date_B),
    
    # ICD-coded TBI kept if:
    # - has at least one ICD code and date
    # - And (no MH OR TBI before MH OR same year as MH)
    keep_icd_TBI = n_ST_codes > 0 & !is.na(injury_date_B) &
      (is.na(any_mh_date) |
         injury_date_B < any_mh_date |
         icd_year == mh_year),
    
    # Self-reported TBI kept if:
    # - self-report Yes
    # - and (no MH date OR MH in/after 2010)
    keep_self_TBI = tbi_self == "Yes" &
      (is.na(any_mh_date) | mh_year >= 2010),
    
    # Final combined TBI exposure
    TBI_final = if_else(keep_self_TBI | keep_icd_TBI, "Yes", "No")
  )

# Make TBI factor with "No" as reference
neuro$TBI_final <- factor(neuro$TBI_final, levels = c("No", "Yes"))
neuro$tbi       <- neuro$TBI_final
neuro$tbi       <- relevel(neuro$tbi, ref = "No")

table(neuro$tbi)

# Check: only ICD TBIs BEFORE/same-year MH or no MH kept as TBI 
neuro %>%
  filter(n_ST_codes > 0, !is.na(injury_date_B)) %>%
  summarise(
    invalid = sum(
      injury_date_B > any_mh_date &
        year(injury_date_B) != year(any_mh_date),
      na.rm = TRUE
    )
  )
xtabs( ~ tbi_report, data = neuro)

# Factorise outcomes 
neuro$dementia        <- factor(neuro$dementia,       levels = c("Yes", "No"))
neuro$stroke          <- factor(neuro$stroke,         levels = c("Yes", "No"))
neuro$parkinsons      <- factor(neuro$parkinsons,     levels = c("Yes", "No"))
neuro$mood_disorder   <- factor(neuro$mood_disorder,  levels = c("Yes", "No"))
neuro$anx_disorder    <- factor(neuro$anx_disorder,   levels = c("Yes", "No"))
neuro$psyc_disorder   <- factor(neuro$psyc_disorder,  levels = c("Yes", "No"))
neuro$any_mh_disorder <- factor(neuro$any_mh_disorder,levels = c("Yes", "No"))


################################################################
# 2. CREATE "AFTER" VARIABLES & MODEL TBI vs MH OUTCOMES BY NEURO CONDITION
################################################################

# Function to create *_after_* variables 
make_after_var <- function(data, outcome_col, condition_col) {
  outcome_col   <- as.character(outcome_col)
  condition_col <- as.character(condition_col)
  
  new_name <- paste0(
    gsub("_date", "", outcome_col),
    "_after_",
    gsub("_date", "", condition_col)
  )
  
  data[[new_name]] <- ifelse(
    is.na(data[[outcome_col]]), "No",
    ifelse(data[[outcome_col]] >= data[[condition_col]], "Yes", "No")
  )
  
  data[[new_name]] <- factor(data[[new_name]], levels = c("No", "Yes"))
  data
}

# Outcomes and conditions (date vars) 
outcomes   <- c("mood_disorder_date", "anx_disorder_date",
                "psyc_disorder_date", "any_mh_date")
conditions <- c("dementia_date", "stroke_date", "parkinsons_date")

# Create all *_after_* variables 
for (outcome in outcomes) {
  for (condition in conditions) {
    neuro <- make_after_var(neuro, outcome, condition)
  }
}

# List of *_after_* outcome names 
all_after_outcomes <- c(
  outer(
    c("mood_disorder", "anx_disorder", "psyc_disorder", "any_mh"),
    c("dementia", "stroke", "parkinsons"),
    paste, sep = "_after_"
  )
)

# Function to extract OR + 95% CI from glmmTMB 
extract_OR_CI <- function(model, term = "tbiYes") {
  coef_table <- summary(model)$coefficients$cond
  beta  <- coef_table[term, "Estimate"]
  se    <- coef_table[term, "Std. Error"]
  OR    <- exp(beta)
  lower <- exp(beta - 1.96 * se)
  upper <- exp(beta + 1.96 * se)
  OR_CI <- paste0(round(OR, 2), " (", round(lower, 2), "-", round(upper, 2), ")")
  pval  <- coef_table[term, "Pr(>|z|)"]
  
  data.frame(
    OR_CI = OR_CI,
    p.value = pval,
    stringsAsFactors = FALSE
  )
}

# Unadjusted models function
run_model_unadj <- function(outcome) {
  f     <- as.formula(paste0(outcome, " ~ tbi + (1 | site)"))
  model <- glmmTMB(f, data = neuro, family = binomial)
  df    <- extract_OR_CI(model)
  df$outcome <- outcome
  df
}

# Adjusted models function
run_model_adj <- function(outcome) {
  f <- as.formula(
    paste0(
      outcome,
      " ~ tbi + age + sex + ethnicity + deprivation_index + alcohol_freq + (1 | site)"
    )
  )
  model <- glmmTMB(f, data = neuro, family = binomial)
  df    <- extract_OR_CI(model)
  df$outcome <- outcome
  df
}

# Run all models & combine 
results_unadj <- map_dfr(all_after_outcomes, run_model_unadj)
results_adj   <- map_dfr(all_after_outcomes, run_model_adj)

results_combined <- results_unadj %>%
  left_join(
    results_adj %>%
      select(outcome, OR_CI_adj = OR_CI, p.value_adj = p.value),
    by = "outcome"
  ) %>%
  tidyr::separate(outcome, into = c("disorder", "condition"), sep = "_after_") %>%
  select(disorder, condition, OR_CI_unadj = OR_CI, OR_CI_adj)

# GT table helper 
make_gt_table <- function(cond) {
  results_combined %>%
    filter(condition == cond) %>%
    mutate(
      disorder = recode(
        disorder,
        "mood_disorder" = "Mood disorders",
        "anx_disorder"  = "Anxiety disorders",
        "psyc_disorder" = "Psychotic disorders",
        "any_mh"        = "Any mental health disorder"
      )
    ) %>%
    select(-condition) %>%
    gt() %>%
    tab_header(title = paste("TBI and mental health outcomes after", cond)) %>%
    cols_label(
      disorder    = "Mental Health Disorder",
      OR_CI_unadj = "Unadjusted OR (95% CI)",
      OR_CI_adj   = "Adjusted OR (95% CI)"
    )
}

# Save GT tables 
gt_dementia <- make_gt_table("dementia")
gtsave(gt_dementia, file.path(output_dir, "dementia_before_mh.html"))

gt_stroke <- make_gt_table("stroke")
gtsave(gt_stroke, file.path(output_dir, "stroke_before_mh.html"))

gt_parkinsons <- make_gt_table("parkinsons")
gtsave(gt_parkinsons, file.path(output_dir, "parkinsons_before_mh.html"))


###############################################################
# BREAKDOWN TABLE: How TBI cases were included/excluded
###############################################################

# 1) FULL ICD cohort with TBI codes
icd_all <- neuro %>%
  filter(n_ST_codes > 0, !is.na(injury_date_B))

# 2) ICD kept by rule
icd_kept <- neuro %>%
  filter(keep_icd_TBI == TRUE)

# 3) ICD excluded by temporal rule
icd_removed <- icd_all %>%
  filter(!(keep_icd_TBI))

# 4) SELF-REPORTED cohort
self_all <- neuro %>%
  filter(tbi_self == "Yes")

# 5) Self-report kept
self_kept <- neuro %>%
  filter(keep_self_TBI == TRUE)

# 6) Self-report removed
self_removed <- self_all %>%
  filter(!(keep_self_TBI))

# Combine into summary table
tbi_summary <- tibble(
  Source = c("ICD-coded TBI", "Self-reported TBI"),
  Total_reported = c(nrow(icd_all), nrow(self_all)),
  Kept_in_analysis = c(sum(neuro$keep_icd_TBI, na.rm = TRUE),
                       sum(neuro$keep_self_TBI, na.rm = TRUE)),
  Removed_by_rules = c(nrow(icd_all) - sum(neuro$keep_icd_TBI, na.rm = TRUE),
                       nrow(self_all) - sum(neuro$keep_self_TBI, na.rm = TRUE))
) %>%
  mutate(
    Percent_kept    = round(100 * Kept_in_analysis / Total_reported, 1),
    Percent_removed = round(100 * Removed_by_rules / Total_reported, 1)
  )

tbi_summary_gt <- tbi_summary %>%
  gt() %>%
  tab_header(title = "Supplementary Table. Inclusion and Exclusion of TBI Cases by Classification Rules") %>%
  cols_label(
    Source           = "TBI Source",
    Total_reported   = "Total Reported",
    Kept_in_analysis = "Kept",
    Removed_by_rules = "Removed",
    Percent_kept     = "% Kept",
    Percent_removed  = "% Removed"
  )

gtsave(tbi_summary_gt, file.path(output_dir, "TBI_inclusion_exclusion_summary.html"))
write.csv(tbi_summary, file.path(output_dir, "TBI_inclusion_exclusion_summary.csv"), row.names = FALSE)


################################################################
# 3. CASE TABLES: NEURO CONDITION AFTER MH, IN FULL NEURO DATA
################################################################

make_cases_table <- function(cond) {
  after_cols <- grep(paste0("_after_", cond, "$"), names(neuro), value = TRUE)
  
  disorder_order <- c(
    "Mood disorder",
    "Anxiety disorder",
    "Psychotic disorder",
    "Any mental health disorder"
  )
  
  neuro %>%
    filter(.data[[cond]] == "Yes") %>%
    pivot_longer(
      cols      = all_of(after_cols),
      names_to  = "disorder",
      values_to = "mh_status"
    ) %>%
    group_by(disorder, tbi) %>%
    summarise(
      n_yes   = sum(as.character(mh_status) == "Yes", na.rm = TRUE),
      total   = n(),
      pct_yes = 100 * n_yes / total,
      .groups = "drop"
    ) %>%
    mutate(
      disorder = gsub(paste0("_after_", cond), "", disorder),
      disorder = recode(
        disorder,
        "mood_disorder" = "Mood disorder",
        "anx_disorder"  = "Anxiety disorder",
        "psyc_disorder" = "Psychotic disorder",
        "any_mh"        = "Any mental health disorder"
      ),
      disorder = factor(disorder, levels = disorder_order),
      summary  = paste0(n_yes, "/", total, " (", round(pct_yes, 1), "%)")
    ) %>%
    arrange(disorder) %>%
    select(disorder, tbi, summary) %>%
    pivot_wider(names_from = tbi, values_from = summary) %>%
    rename(
      TBI    = "Yes",
      No_TBI = "No"
    ) %>%
    gt() %>%
    tab_header(title = paste("Cases of TBI and mental health outcomes in", cond)) %>%
    cols_label(
      disorder = "Mental health disorder",
      TBI      = "TBI",
      No_TBI   = "No TBI"
    )
}

cases_dementia   <- make_cases_table("dementia")
gtsave(cases_dementia, file.path(output_dir, "dementia_after_mh_cases.html"))

cases_stroke     <- make_cases_table("stroke")
gtsave(cases_stroke, file.path(output_dir, "stroke_after_mh_cases.html"))

cases_parkinsons <- make_cases_table("parkinsons")
gtsave(cases_parkinsons, file.path(output_dir, "parkinsons_after_mh_cases.html"))


################################################################
# 4. ANALYSIS FOR TBI IN PARTICIPANTS WITH NO NEUROLOGICAL DISORDERS
################################################################

# Load full dataset 
full_data <- read.csv(file.path(output_dir, "data.csv"))

# Standardise neuro vars + create no_neuro_disorder 
full_data <- full_data %>%
  mutate(
    dementia   = toupper(as.character(dementia)),
    stroke     = toupper(as.character(stroke)),
    parkinsons = toupper(as.character(parkinsons))
  ) %>%
  mutate(
    no_neuro_disorder = ifelse(
      dementia == "NO" & stroke == "NO" & parkinsons == "NO",
      "Yes", "No"
    ),
    no_neuro_disorder = factor(no_neuro_disorder, levels = c("No", "Yes"))
  )

# Define MH outcomes 
mh_outcomes <- c("mood_disorder", "anx_disorder", "psyc_disorder", "any_mh_disorder")

# Subset to those with no neuro disorders 
full_data_noneuro <- full_data %>%
  filter(no_neuro_disorder == "Yes")

n_noneuro <- nrow(full_data_noneuro)
if (n_noneuro == 0) {
  stop("No participants with no neurological disorder found. Check variable coding!")
}

# Create numeric 0/1 versions of MH outcomes 
full_data_noneuro <- full_data_noneuro %>%
  mutate(
    mood_disorder_num    = as.integer(mood_disorder == "Yes"),
    anx_disorder_num     = as.integer(anx_disorder  == "Yes"),
    psyc_disorder_num    = as.integer(psyc_disorder == "Yes"),
    any_mh_disorder_num  = as.integer(any_mh_disorder == "Yes")
  )

# Unadjusted logistic regressions (no neuro disorders)
results_unadj_noneuro <- map_dfr(mh_outcomes, function(outcome) {
  f     <- as.formula(paste0(outcome, "_num ~ tbi + (1 | site)"))
  model <- glmmTMB(f, data = full_data_noneuro, family = binomial)
  df    <- extract_OR_CI(model)
  df$outcome <- outcome
  df
})

# Adjusted logistic regressions (no neuro disorders) 
results_adj_noneuro <- map_dfr(mh_outcomes, function(outcome) {
  f <- as.formula(
    paste0(
      outcome,
      "_num ~ tbi + age + sex + ethnicity + deprivation_index + alcohol_freq + (1 | site)"
    )
  )
  model <- glmmTMB(f, data = full_data_noneuro, family = binomial)
  df    <- extract_OR_CI(model)
  df$outcome <- outcome
  df
})

# Combine and GT table 
results_noneuro_combined <- results_unadj_noneuro %>%
  left_join(
    results_adj_noneuro %>%
      select(outcome, OR_CI_adj = OR_CI, p.value_adj = p.value),
    by = "outcome"
  ) %>%
  mutate(condition = "no_neuro_disorder") %>%
  select(outcome, condition, OR_CI_unadj = OR_CI, OR_CI_adj)

gt_noneuro <- results_noneuro_combined %>%
  mutate(
    outcome = recode(
      outcome,
      "mood_disorder"    = "Mood disorders",
      "anx_disorder"     = "Anxiety disorders",
      "psyc_disorder"    = "Psychotic disorders",
      "any_mh_disorder"  = "Any mental health disorder"
    )
  ) %>%
  gt() %>%
  tab_header(
    title = "TBI and mental health outcomes in participants with no neurological disorders"
  ) %>%
  cols_label(
    outcome     = "Mental Health Disorder",
    OR_CI_unadj = "Unadjusted OR (95% CI)",
    OR_CI_adj   = "Adjusted OR (95% CI)"
  )

gtsave(gt_noneuro, file.path(output_dir, "no_neuro_disorder_before_mh.html"))

# Case table: no neuro disorders 
cases_noneuro <- full_data_noneuro %>%
  pivot_longer(
    cols      = all_of(mh_outcomes),
    names_to  = "disorder",
    values_to = "mh_status"
  ) %>%
  group_by(disorder, tbi) %>%
  summarise(
    n_yes   = sum(mh_status == "Yes", na.rm = TRUE),
    total   = n(),
    pct_yes = 100 * n_yes / total,
    .groups = "drop"
  ) %>%
  mutate(
    summary = paste0(n_yes, "/", total, " (", round(pct_yes, 1), "%)"),
    disorder = recode(
      disorder,
      "mood_disorder"    = "Mood disorder",
      "anx_disorder"     = "Anxiety disorder",
      "psyc_disorder"    = "Psychotic disorder",
      "any_mh_disorder"  = "Any mental health disorder"
    )
  ) %>%
  select(disorder, tbi, summary) %>%
  pivot_wider(names_from = tbi, values_from = summary) %>%
  rename(
    TBI    = "Yes",
    No_TBI = "No"
  ) %>%
  gt() %>%
  tab_header(
    title = "Cases of TBI and mental health outcomes in participants with no neurological disorders"
  ) %>%
  cols_label(
    disorder = "Mental health disorder",
    TBI      = "TBI",
    No_TBI   = "No TBI"
  )

gtsave(cases_noneuro, file.path(output_dir, "no_neuro_disorder_cases.html"))


################################################################
# 5. TBI × NEURO 4-GROUP VARIABLE (ALL AGES)
################################################################

# Create 4-group TBI × neuro variable 
full_data <- full_data %>%
  mutate(
    TBI_neuro_group = case_when(
      tbi == "Yes" & no_neuro_disorder == "No"  ~ "TBI + neuro",
      tbi == "Yes" & no_neuro_disorder == "Yes" ~ "TBI - neuro",
      tbi == "No"  & no_neuro_disorder == "No"  ~ "No TBI + neuro",
      tbi == "No"  & no_neuro_disorder == "Yes" ~ "No TBI - neuro"
    ),
    TBI_neuro_group = factor(
      TBI_neuro_group,
      levels = c("No TBI - neuro", "TBI - neuro", "No TBI + neuro", "TBI + neuro")
    )
  )




# Unadjusted ORs, only TBI_neuro_group contrasts 
run_unadj_OR_TBIgroups <- function(outcome, data = full_data) {
  df <- data %>%
    mutate(out_num = as.integer(.data[[outcome]] == "Yes"))
  
  # "Unadjusted" = exposure only (optionally keep random intercept for site)
  f <- as.formula("out_num ~ TBI_neuro_group + (1 | site)")
  model <- glmmTMB(f, data = df, family = binomial)
  
  coef_table <- summary(model)$coefficients$cond
  keep_rows  <- grep("^TBI_neuro_group", rownames(coef_table))
  if (length(keep_rows) == 0) return(NULL)
  
  coef_table <- coef_table[keep_rows, , drop = FALSE]
  beta  <- coef_table[, "Estimate"]
  se    <- coef_table[, "Std. Error"]
  OR    <- exp(beta)
  lower <- exp(beta - 1.96 * se)
  upper <- exp(beta + 1.96 * se)
  pval  <- coef_table[, "Pr(>|z|)"]
  
  tibble(
    outcome = outcome,
    group   = gsub("TBI_neuro_group", "", rownames(coef_table)),
    OR      = round(OR, 2),
    CI      = paste0(round(lower, 2), "-", round(upper, 2)),
    p.value = round(pval, 3),
    model   = "Unadjusted"
  )
}

# Adjusted ORs, only TBI_neuro_group contrasts 
run_adj_OR_TBIgroups <- function(outcome, data = full_data) {
  df <- data %>%
    mutate(out_num = as.integer(.data[[outcome]] == "Yes"))
  
  f     <- as.formula(
    "out_num ~ TBI_neuro_group + age + sex + ethnicity + deprivation_index + alcohol_freq + (1 | site)"
  )
  model <- glmmTMB(f, data = df, family = binomial)
  
  coef_table <- summary(model)$coefficients$cond
  keep_rows  <- grep("^TBI_neuro_group", rownames(coef_table))
  
  if (length(keep_rows) == 0) return(NULL)
  
  coef_table <- coef_table[keep_rows, , drop = FALSE]
  beta  <- coef_table[, "Estimate"]
  se    <- coef_table[, "Std. Error"]
  OR    <- exp(beta)
  lower <- exp(beta - 1.96 * se)
  upper <- exp(beta + 1.96 * se)
  pval  <- coef_table[, "Pr(>|z|)"]
  
  tibble(
    outcome = outcome,
    group   = gsub("TBI_neuro_group", "", rownames(coef_table)),
    OR      = round(OR, 2),
    CI      = paste0(round(lower, 2), "-", round(upper, 2)),
    p.value = round(pval, 3), 
    model = "Adjusted"
  )
}

results_unadj_OR <- map_dfr(mh_outcomes, run_unadj_OR_TBIgroups)
results_adj_OR <- map_dfr(mh_outcomes, run_adj_OR_TBIgroups)

results_both <- bind_rows(results_unadj_OR, results_adj_OR) %>%
  mutate(
    outcome = recode(
      outcome,
      "mood_disorder"    = "Mood disorders",
      "anx_disorder"     = "Anxiety disorders",
      "psyc_disorder"    = "Psychotic disorders",
      "any_mh_disorder"  = "Any mental health disorder"
    ),
    group = trimws(group), 
    model = factor(model, levels = c("Unadjusted", "Adjusted"))
  )
write.csv(results_both, 
          file.path(output_dir, "TBI_neuro_groups_unadj_adj_ORs_clean.csv"),
row.names = FALSE)



gt_adj_OR <- results_adj_OR %>%
  gt() %>%
  tab_header(title = "Adjusted ORs for Mental Health Outcomes by TBI + Neuro Groups") %>%
  cols_label(
    outcome = "Mental Health Disorder",
    group   = "Group (vs No TBI - neuro)",
    OR      = "Adjusted OR",
    CI      = "95% CI",
    p.value = "p-value"
  )

gtsave(gt_adj_OR, file.path(output_dir, "TBI_neuro_groups_adjusted_ORs_clean.html"))
write.csv(results_adj_OR,
          file.path(output_dir, "TBI_neuro_groups_adjusted_ORs_clean.csv"),
          row.names = FALSE)

# Descriptive prevalence table 
desc_table <- full_data %>%
  pivot_longer(
    cols      = all_of(mh_outcomes),
    names_to  = "disorder",
    values_to = "mh_status"
  ) %>%
  group_by(disorder, TBI_neuro_group) %>%
  summarise(
    n_yes   = sum(mh_status == "Yes", na.rm = TRUE),
    total   = n(),
    pct_yes = 100 * n_yes / total,
    .groups = "drop"
  ) %>%
  mutate(
    summary = paste0(n_yes, "/", total, " (", round(pct_yes, 1), "%)"),
    disorder = recode(
      disorder,
      "mood_disorder"    = "Mood disorders",
      "anx_disorder"     = "Anxiety disorders",
      "psyc_disorder"    = "Psychotic disorders",
      "any_mh_disorder"  = "Any mental health disorder"
    )
  ) %>%
  select(disorder, TBI_neuro_group, summary) %>%
  pivot_wider(names_from = TBI_neuro_group, values_from = summary)

gt_desc <- desc_table %>%
  gt() %>%
  tab_header(title = "Prevalence of Mental Health Disorders by TBI + Neuro Group") %>%
  cols_label(
    disorder         = "Mental Health Disorder",
    `No TBI - neuro` = "No TBI - neuro",
    `TBI - neuro`    = "TBI - neuro",
    `No TBI + neuro` = "No TBI + neuro",
    `TBI + neuro`    = "TBI + neuro"
  )

gtsave(gt_desc, file.path(output_dir, "TBI_neuro_groups_percentages.html"))
write.csv(desc_table,
          file.path(output_dir, "TBI_neuro_groups_percentages.csv"),
          row.names = FALSE)


################################################################
# 6. ANALYTIC COHORT (SUBSET OF NEURO) + CASE TABLES
################################################################

analytic <- neuro %>%
  filter(
    # keep everyone with no MH date
    is.na(any_mh_date) |
      # TBI before or same year as MH 
      is.na(injury_date_B) | injury_date_B <= any_mh_date
  ) %>%
  filter(
    # Neurological conditions must occur before or same year as MH (or not occur at all) 
    is.na(any_mh_date) | is.na(dementia_date)   | dementia_date   <= any_mh_date, 
    is.na(any_mh_date) | is.na(stroke_date)     | stroke_date     <= any_mh_date, 
    is.na(any_mh_date) | is.na(parkinsons_date) | parkinsons_date <= any_mh_date
  )

nrow(analytic)  
table(analytic$tbi)
with(analytic, table(tbi, is.na(any_mh_date)))
xtabs(~ dementia, data = analytic) 

# Demographic and descriptive tables 
descriptive_table_full <- analytic %>%
  select(age, sex, deprivation_index, ethnicity, alcohol_freq, drug_addict, alcohol_addict) %>%
  tbl_summary() %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_caption("**Table 1. Demographics of sample**")

descriptive_table_tbi <- analytic %>%
  select(tbi, age, sex, deprivation_index, ethnicity, alcohol_freq, drug_addict, alcohol_addict) %>%
  tbl_summary(by = tbi) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**TBI**") 

descriptive_table_dementia <- analytic %>%
  select(dementia, age, sex, deprivation_index, ethnicity, alcohol_freq, drug_addict, alcohol_addict) %>%
  tbl_summary(by = dementia) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**Dementia**") 

descriptive_table_stroke <- analytic %>%
  select(stroke, age, sex, deprivation_index, ethnicity, alcohol_freq, drug_addict, alcohol_addict) %>%
  tbl_summary(by = stroke) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**Stroke**") 

descriptive_table_parkinsons <- analytic %>%
  select(parkinsons, age, sex, deprivation_index, ethnicity, alcohol_freq, drug_addict, alcohol_addict) %>%
  tbl_summary(by = parkinsons) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**Parkinsons**")

tbl <- tbl_merge(
  tbls = list(
    descriptive_table_full,
    descriptive_table_tbi,
    descriptive_table_dementia,
    descriptive_table_stroke,
    descriptive_table_parkinsons
  ), 
  tab_spanner = c("**Overall**","**TBI**", "**Dementia**" , "**Stroke**", "**Parkinsons**")
)

table_1_filename <- file.path(output_dir, "final_cohort_descriptive_table.html")
gt::gtsave(as_gt(tbl), file = table_1_filename)


# MH disorder table by TBI 
psyc_table_tbi <- analytic %>%
  select(tbi, mood_disorder, anx_disorder, psyc_disorder, any_mh_disorder, tbi_report) %>%
  tbl_summary(by = tbi) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_spanning_header(all_stat_cols() ~ "**TBI**")

gt::gtsave(as_gt(psyc_table_tbi),
           file = file.path(output_dir, "psyc_table_tbi_analytic.html"))


# Case tables in analytic cohort by neuro condition 
make_cases_table_analytic <- function(cond) {
  after_cols <- grep(paste0("_after_", cond, "$"), names(analytic), value = TRUE)
  
  disorder_order <- c(
    "Mood disorder",
    "Anxiety disorder",
    "Psychotic disorder",
    "Any mental health disorder"
  )
  
  analytic %>%
    filter(.data[[cond]] == "Yes") %>%
    pivot_longer(
      cols      = all_of(after_cols),
      names_to  = "disorder",
      values_to = "mh_status"
    ) %>%
    group_by(disorder, tbi) %>%
    summarise(
      n_yes   = sum(as.character(mh_status) == "Yes", na.rm = TRUE),
      total   = n(),
      pct_yes = 100 * n_yes / total,
      .groups = "drop"
    ) %>%
    mutate(
      disorder = gsub(paste0("_after_", cond), "", disorder),
      disorder = recode(
        disorder,
        "mood_disorder" = "Mood disorder",
        "anx_disorder"  = "Anxiety disorder",
        "psyc_disorder" = "Psychotic disorder",
        "any_mh"        = "Any mental health disorder"
      ),
      disorder = factor(disorder, levels = disorder_order),
      summary  = paste0(n_yes, "/", total, " (", round(pct_yes, 1), "%)")
    ) %>%
    arrange(disorder) %>%
    select(disorder, tbi, summary) %>%
    pivot_wider(names_from = tbi, values_from = summary) %>%
    rename(
      TBI    = "Yes",
      No_TBI = "No"
    ) %>%
    gt() %>%
    tab_header(
      title = paste("Cases of TBI and mental health outcomes in", cond, "(analytic cohort)")
    ) %>%
    cols_label(
      disorder = "Mental health disorder",
      TBI      = "TBI",
      No_TBI   = "No TBI"
    )
}

cases_dementia_analytic   <- make_cases_table_analytic("dementia")
gtsave(cases_dementia_analytic,
       file.path(output_dir, "dementia_after_mh_cases_analytic.html"))

cases_stroke_analytic     <- make_cases_table_analytic("stroke")
gtsave(cases_stroke_analytic,
       file.path(output_dir, "stroke_after_mh_cases_analytic.html"))

cases_parkinsons_analytic <- make_cases_table_analytic("parkinsons")
gtsave(cases_parkinsons_analytic,
       file.path(output_dir, "parkinsons_after_mh_cases_analytic.html"))


################################################################
# 7. FOREST PLOT FOR TBI_NEURO_GROUP ADJUSTED ORs
################################################################

results_adj_OR_plot <- read.csv(
  file.path(output_dir, "TBI_neuro_groups_adjusted_ORs_clean.csv")
)

forest_dat <- results_adj_OR_plot %>%
  # Exclude psychotic disorders from plot 
  filter(outcome != "Psychotic disorders") %>%
  separate(CI, into = c("CI_low", "CI_high"), sep = "-", remove = FALSE) %>%
  mutate(
    OR      = as.numeric(OR),
    CI_low  = as.numeric(CI_low),
    CI_high = as.numeric(CI_high),
    group   = factor(
      group,
      levels = c("TBI + neuro", "TBI - neuro", "No TBI + neuro")
    ),
    outcome = factor(
      outcome,
      levels = c(
        "Mood disorders",
        "Anxiety disorders",
        "Any mental health disorder"
      )
    )
  )

p <- ggplot(forest_dat, aes(x = OR, y = group)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2) +
  scale_x_log10() +
  scale_y_discrete(
    labels = c(
      "TBI - neuro"    = "TBI only",
      "No TBI + neuro" = "Neuro only",
      "TBI + neuro"    = "TBI + Neuro"
    )
  ) +
  facet_wrap(~ outcome, ncol = 1) +
  labs(
    x     = "Adjusted odds ratio (log scale)",
    y     = "Group (Reference = No TBI or Neuro)",
    title = "Adjusted ORs for mental health outcomes by TBI + neuro group"
  ) +
  theme_bw() +
  theme(
    strip.text        = element_text(face = "bold"),
    panel.grid.minor  = element_blank()
  )

ggsave(file.path(output_dir, "forest_plot_tbi.pdf"),
       p, width = 6, height = 10)


################################################################
# 8. TBI CASE TABLE IN ANALYTIC COHORT
################################################################

xtabs(~ tbi + tbi_report, data = analytic)

tbi_cases_table <- analytic %>%
  pivot_longer(
    cols      = all_of(mh_outcomes),
    names_to  = "disorder",
    values_to = "mh_status"
  ) %>%
  group_by(disorder, tbi) %>%
  summarise(
    n_yes   = sum(mh_status == "Yes", na.rm = TRUE),
    total   = n(),
    pct_yes = 100 * n_yes / total,
    .groups = "drop"
  ) %>%
  mutate(
    disorder = recode(
      disorder,
      "mood_disorder"    = "Mood disorder",
      "anx_disorder"     = "Anxiety disorder",
      "psyc_disorder"    = "Psychotic disorder",
      "any_mh_disorder"  = "Any mental health disorder"
    ),
    summary = paste0(n_yes, "/", total, " (", round(pct_yes, 1), "%)")
  ) %>%
  select(disorder, tbi, summary) %>%
  pivot_wider(names_from = tbi, values_from = summary) %>%
  rename(
    TBI    = "Yes",
    No_TBI = "No"
  ) %>%
  gt() %>%
  tab_header(title = "Mental health outcomes by TBI status (analytic cohort)") %>%
  cols_label(
    disorder = "Mental health disorder",
    TBI      = "TBI",
    No_TBI   = "No TBI"
  )

tbi_cases_table

message("All tables and plots saved to: ", output_dir)






















