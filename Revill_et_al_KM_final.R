###############################################################

# Analysis of Biobank dataset

# Revill et al

################################################################

# Load libraries 
library(survminer)
library(gridExtra)
library(dplyr)
library(purrr)
library(survival)
library(tidyr)
library(lubridate)


# Remove everything from memory if needed to start afresh
#rm(list = ls())

# Set data and output directories
data_dir   <- "S:/Track_TBI/Grace/Test/"
output_dir <- "S:/Track_TBI/Grace/Test/Output/"
setwd(data_dir)

# Load dataset 
neuro <- read.csv(file.path(output_dir, "neuro.csv"))

# Convert all "*date" columns back to Date 
date_cols <- grep("date", names(neuro), value = TRUE, ignore.case = TRUE)
neuro[date_cols] <- lapply(neuro[date_cols], as.Date)

if (!"injury_date_B" %in% names(neuro)) {
  stop("injury_date_B not found in neuro.csv - this script expects it.")
}

################################################################
# 1. RECONSTRUCT FINAL TBI VARIABLE (SAME AS MAIN ANALYSIS)
################################################################

neuro <- neuro %>%
  mutate(
    mh_year  = year(any_mh_date),
    icd_year = year(injury_date_B),
    
    # ICD-coded TBI kept if:
    # - S/T codes present AND injury_date_B known
    # AND (no MH date OR TBI date < MH date OR same calendar year)
    keep_icd_TBI = n_ST_codes > 0 & !is.na(injury_date_B) &
      (is.na(any_mh_date) |
         injury_date_B < any_mh_date |
         icd_year == mh_year),
    
    # Self-reported TBI kept if:
    # - self-report "Yes"
    # AND (no MH date OR MH in/after 2010, when self-report was collected)
    keep_self_TBI = tbi_self == "Yes" &
      (is.na(any_mh_date) | mh_year >= 2010),
    
    # Final combined exposure
    TBI_final = if_else(keep_self_TBI | keep_icd_TBI, "Yes", "No"),
    TBI_final = factor(TBI_final, levels = c("No", "Yes")),
    tbi       = TBI_final
  )

# Check: no ICD TBIs clearly after MH in later year 
bad_icd <- sum(
  neuro$tbi == "Yes" &
    !is.na(neuro$injury_date_B) &
    !is.na(neuro$any_mh_date) &
    neuro$injury_date_B > neuro$any_mh_date &
    year(neuro$injury_date_B) != year(neuro$any_mh_date),
  na.rm = TRUE
)
message("ICD TBIs after MH (different year): ", bad_icd)

# Ensure tbi factor with "No" as reference 
neuro$tbi <- factor(neuro$tbi, levels = c("No", "Yes"))
neuro$tbi <- relevel(neuro$tbi, ref = "No")

# Drop rows missing TBI classification (just in case)
neuro <- neuro %>% filter(!is.na(tbi))

################################################################
# 2. DEFINE CONDITIONS, OUTCOMES
################################################################

conditions <- c("dementia", "stroke", "parkinsons")
# use same stems as main script: *_date exist for each of these
outcomes   <- c("mood_disorder", "anx_disorder", "psyc_disorder", "any_mh")

study_end <- as.Date("2022-12-31")
max_time  <- 10  # years

# Label outcomes 
label_outcome <- function(x) {
  recode(x, 
         "mood_disorder" = "Mood disorder", 
         "anx_disorder" = "Anxiety disorder", 
         "psyc_disorder" = "Psychotic disorder", 
         "any_mh" = "Any mental health disorder")
}

label_condition <- function(x) {
  recode(x, 
         "dementia" = "dementia", 
         "stroke" = "stroke", 
         "parkinsons" = "Parkinson's disease")
}

################################################################
# 3. BUILD SURVIVAL DATA
#    - keep only people with the neuro condition
#    - EXCLUDE those with MH outcome *before* that condition
################################################################

make_surv_data <- function(data, condition, outcome) {
  
  cond_date    <- paste0(condition, "_date")
  outcome_date <- paste0(outcome, "_date")
  
  df <- data %>%
    # must have the neuro diagnosis date
    filter(!is.na(.data[[cond_date]])) %>%
    # EXCLUDE people who already had the MH outcome before the neuro condition
    filter(is.na(.data[[outcome_date]]) |
             .data[[outcome_date]] >= .data[[cond_date]])
  
  if (nrow(df) == 0) return(NULL)
  
  df <- df %>%
    mutate(
      # time from neuro condition to MH outcome (for events)
      time = as.numeric(
        difftime(.data[[outcome_date]], .data[[cond_date]], units = "days")
      ) / 365.25,
      
      # event if outcome occurs on/after condition date
      event = ifelse(!is.na(.data[[outcome_date]]) & time >= 0, 1, 0),
      
      # if no event, censor at study_end
      time = ifelse(
        event == 0,
        as.numeric(difftime(study_end, .data[[cond_date]], units = "days")) / 365.25,
        time
      ),
      
      condition = condition,
      outcome   = outcome
    )
  
  if (sum(df$event) == 0) return(NULL)
  df
}

safe_make_surv <- safely(make_surv_data)

surv_data <- purrr::map_dfr(conditions, function(cond) {
  purrr::map(outcomes, function(out) {
    res <- safe_make_surv(neuro, cond, out)
    if (!is.null(res$result)) res$result else NULL
  }) %>% compact()
})

message("Survival dataset built. Rows: ", nrow(surv_data))
print(table(surv_data$condition, surv_data$outcome))

################################################################
# 4. KAPLAN-MEIER CURVES + NUMBER AT RISK
################################################################

cond_out_pairs <- surv_data %>% distinct(condition, outcome)

km_fits <- purrr::map(seq_len(nrow(cond_out_pairs)), function(i) {
  cond <- cond_out_pairs$condition[i]
  out  <- cond_out_pairs$outcome[i]
  
  subset_df <- surv_data %>%
    filter(condition == cond, outcome == out)
  
  if (nrow(subset_df) == 0 || sum(subset_df$event) == 0) {
    return(NULL)
  } else {
    survfit(Surv(time, event) ~ tbi, data = subset_df)
  }
})

names(km_fits) <- paste(cond_out_pairs$condition,
                        cond_out_pairs$outcome,
                        sep = "_")

km_dir <- file.path(output_dir, "km_plots")
if (!dir.exists(km_dir)) {
  dir.create(km_dir, recursive = TRUE)
}

risk_times <- seq(0, max_time, by = 2)

purrr::walk(names(km_fits), function(name) {
  
  fit <- km_fits[[name]]
  if (is.null(fit)) return(NULL)
  
  parts <- strsplit(name, "_")[[1]]
  cond  <- parts[1]
  out   <- paste(parts[-1], collapse = "_")
  
  surv_subset <- surv_data %>%
    filter(condition == cond, outcome == out)
  
  if (nrow(surv_subset) == 0 || sum(surv_subset$event) == 0) return(NULL)
  
  pretty_out  <- label_outcome(out)
  pretty_cond <- label_condition(cond)
  ylab_text <- paste0(pretty_out, " in ", pretty_cond)
  
  
  
  gp <- ggsurvplot(
    fit,
    data              = surv_subset,
    fun               = "event",       # cumulative incidence
    conf.int          = FALSE,
    pval              = TRUE,
    pval.size         = 7,
    xlim              = c(0, max_time),
    break.time.by     = 2,
    censor            = TRUE,
    ylab              = ylab_text,
    xlab              = "Time (years)",
    legend.labs       = c("No TBI", "TBI"),
    legend.title      = "TBI status",
    ggtheme           = theme_minimal(base_size =18) + theme(
      axis.text   = element_text(size = 18),
      axis.title  = element_text(size = 18),
      legend.position = "bottom"
    ),
    size              = 1,
    risk.table        = TRUE,          # number at risk panel
    risk.table.fontsize = 7,
    risk.table.height = 0.25,
    risk.table.y.text = TRUE
  )
  
  # y-axis as percentages
  gp$plot <- gp$plot +
    scale_y_continuous(labels = function(x) paste0(round(100 * x, 0), "%"))
  
  png(file.path(km_dir, paste0("KM_", name, ".png")),
      width = 1600, height = 1200, res = 150)
  print(gp)
  dev.off()
})

message("KM plots saved in: ", km_dir)

################################################################
# 5. COX PROPORTIONAL HAZARDS MODELS (USING SAME surv_data)
################################################################

# candidate adjustment covariates
candidate_covars <- c("age", "sex", "ethnicity", "deprivation_index", "alcohol_freq")
adj_covars <- intersect(candidate_covars, names(neuro))

# helper to get HR/CI/p for tbi term
extract_tbi_hr <- function(fit) {
  s  <- summary(fit)
  rn <- rownames(s$coefficients)
  tbi_row <- grep("^tbi", rn)   # e.g. "tbiYes"
  if (length(tbi_row) != 1) return(NULL)
  
  beta  <- s$coefficients[tbi_row, "coef"]
  se    <- s$coefficients[tbi_row, "se(coef)"]
  hr    <- exp(beta)
  lower <- exp(beta - 1.96 * se)
  upper <- exp(beta + 1.96 * se)
  pval  <- s$coefficients[tbi_row, "Pr(>|z|)"]
  
  data.frame(
    HR      = hr,
    CI_low  = lower,
    CI_high = upper,
    p_value = pval,
    rowname = rn[tbi_row],
    stringsAsFactors = FALSE
  )
}

cox_results <- map_dfr(seq_len(nrow(cond_out_pairs)), function(i) {
  cond <- cond_out_pairs$condition[i]
  out  <- cond_out_pairs$outcome[i]
  
  surv_subset <- surv_data %>%
    filter(condition == cond, outcome == out)
  
  # need events and both TBI levels
  if (nrow(surv_subset) == 0 ||
      sum(surv_subset$event) == 0 ||
      length(unique(surv_subset$tbi)) < 2) return(NULL)
  
  ## Unadjusted Cox
  cox_unadj <- coxph(Surv(time, event) ~ tbi, data = surv_subset)
  hr_unadj  <- extract_tbi_hr(cox_unadj)
  ph_unadj  <- try(cox.zph(cox_unadj), silent = TRUE)
  ph_p_unadj <- if (inherits(ph_unadj, "try-error")) NA_real_
  else ph_unadj$table["GLOBAL", "p"]
  
  ## Adjusted Cox
  if (length(adj_covars) > 0) {
    f_adj <- as.formula(
      paste("Surv(time, event) ~ tbi +", paste(adj_covars, collapse = " + "))
    )
  } else {
    f_adj <- as.formula("Surv(time, event) ~ tbi")
  }
  
  cox_adj <- coxph(f_adj, data = surv_subset)
  hr_adj  <- extract_tbi_hr(cox_adj)
  ph_adj  <- try(cox.zph(cox_adj), silent = TRUE)
  ph_p_adj <- if (inherits(ph_adj, "try-error")) NA_real_
  else ph_adj$table["GLOBAL", "p"]
  
  tibble(
    condition = cond,
    outcome   = out,
    
    HR_unadj  = round(hr_unadj$HR, 2),
    CI_unadj  = paste0("(", round(hr_unadj$CI_low, 2), " - ",
                       round(hr_unadj$CI_high, 2), ")"),
    PH_p_unadj = signif(ph_p_unadj, 3),
    PH_violation_unadj = ifelse(!is.na(ph_p_unadj) & ph_p_unadj < 0.05, "Yes", "No"),
    
    HR_adj    = round(hr_adj$HR, 2),
    CI_adj    = paste0("(", round(hr_adj$CI_low, 2), " - ",
                       round(hr_adj$CI_high, 2), ")"),
    PH_p_adj  = signif(ph_p_adj, 3),
    PH_violation_adj = ifelse(!is.na(ph_p_adj) & ph_p_adj < 0.05, "Yes", "No")
  )
})

# Save Cox results 
cox_file <- file.path(output_dir, "cox_results_km_subset.csv")
write.csv(cox_results, cox_file, row.names = FALSE)
message("Cox results saved to: ", cox_file)

cox_results






