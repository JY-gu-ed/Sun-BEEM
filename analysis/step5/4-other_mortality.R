if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(gtsummary)
library(dplyr)
library(survival)
library(survival)
library(broom)
setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step5/"))

getwd()

morta <- read.csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "washed_data", "4_2_2_ana_for_mortality.csv"))

morta_backup <- morta

morta <- morta %>%
  mutate(
    bi_outdoor     = factor(bi_outdoor, levels = c(0,1), labels = c("No","Yes")),
    sex            = factor(sex, levels = c("female","male")),
    age            = factor(age, levels = c("37-47","48-54","55-59","60-63","64-73")),
    bmi            = factor(bmi, levels = c("underweight","healthy weight","overweight","obese")),
    sleep_score    = factor(sleep_score, levels = c("poor","intermediate","healthy")),
    drinking       = factor(drinking, levels = c("unhealthy","healthy")),
    phys_activity  = factor(phys_activity, levels = c("low","moderate","high")),
    highest_edu    = factor(highest_edu, levels = c("less than secondary","secondary school","university or above")),
    smoking_status = factor(smoking_status, levels = c("current","previous","never")),
    deprivation_quintile = factor(deprivation_quintile, levels = 1:5, ordered = TRUE)
  )

mortality_vars <- c(
  "all_cause_mortality",
  "cvd_mortality_narrow",
  "cancer_mortality",
  "skin_cancer_mortality",
  "melanoma_skin_cancer_mortality",
  "non_melanoma_mortality",
  "non_skin_cancer_mortality"
)

morta <- morta %>%
  mutate(across(all_of(mortality_vars), ~ ifelse(. == "yes", 1, 0)))

library(survival)
library(broom)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

covariates <- c("sex","age","bmi","deprivation_quintile",
                "highest_edu","phys_activity","sleep_score",
                "smoking_status","drinking")

safe_fit <- function(formula, data) {
  tryCatch(coxph(formula, data = data, na.action = na.exclude), error = function(e) NULL)
}

run_one <- function(outcome) {
  f_crude <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ bi_outdoor"))
  fit_crude <- safe_fit(f_crude, morta)
  
  crude_df <- if (!is.null(fit_crude)) {
    tidy(fit_crude, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "bi_outdoorYes") %>%
      mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
             p = sprintf("%.4f", p.value)) %>%
      select(HR_CI, p)
  } else tibble(HR_CI = NA, p = NA)
  
  rhs <- paste(c("bi_outdoor", covariates), collapse = " + ")
  f_adj <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ ", rhs))
  fit_adj <- safe_fit(f_adj, morta)
  
  adj_df <- if (!is.null(fit_adj)) {
    tidy(fit_adj, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "bi_outdoorYes") %>%
      mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
             p = sprintf("%.4f", p.value)) %>%
      select(HR_CI, p)
  } else tibble(HR_CI = NA, p = NA)
  
  tibble(
    outcome         = outcome,
    crude_HR_Yes_vs_No = crude_df$HR_CI,
    crude_p         = crude_df$p,
    adj_HR_Yes_vs_No   = adj_df$HR_CI,
    adj_p           = adj_df$p
  )
}

final_tbl <- map_dfr(mortality_vars, run_one)
final_tbl

write.csv(final_tbl, "../result/5_4_1_cox_bi_outdoor_crude_adjusted.csv", row.names = FALSE)

morta <- morta_backup

morta <- morta %>%
  mutate(
    bi_uv_protection = factor(bi_uv_protection, levels = c(0,1), labels = c("No","Yes")),
    sex            = factor(sex, levels = c("female","male")),
    age            = factor(age, levels = c("37-47","48-54","55-59","60-63","64-73")),
    bmi            = factor(bmi, levels = c("underweight","healthy weight","overweight","obese")),
    sleep_score    = factor(sleep_score, levels = c("poor","intermediate","healthy")),
    drinking       = factor(drinking, levels = c("unhealthy","healthy")),
    phys_activity  = factor(phys_activity, levels = c("low","moderate","high")),
    highest_edu    = factor(highest_edu, levels = c("less than secondary","secondary school","university or above")),
    smoking_status = factor(smoking_status, levels = c("current","previous","never")),
    deprivation_quintile = factor(deprivation_quintile, levels = 1:5, ordered = TRUE)
  )

mortality_vars <- c(
  "all_cause_mortality",
  "cvd_mortality_narrow",
  "cancer_mortality",
  "skin_cancer_mortality",
  "melanoma_skin_cancer_mortality",
  "non_melanoma_mortality",
  "non_skin_cancer_mortality"
)

morta <- morta %>%
  mutate(across(all_of(mortality_vars), ~ ifelse(. == "yes", 1, 0)))

covariates <- c("sex","age","bmi","deprivation_quintile",
                "highest_edu","phys_activity","sleep_score",
                "smoking_status","drinking")

safe_fit <- function(formula, data) {
  tryCatch(coxph(formula, data = data, na.action = na.exclude), error = function(e) NULL)
}

run_one <- function(outcome) {
  f_crude <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ bi_uv_protection"))
  fit_crude <- safe_fit(f_crude, morta)
  
  crude_df <- if (!is.null(fit_crude)) {
    tidy(fit_crude, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "bi_uv_protectionYes") %>%
      mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
             p = sprintf("%.4f", p.value)) %>%
      select(HR_CI, p)
  } else tibble(HR_CI = NA, p = NA)
  
  rhs <- paste(c("bi_uv_protection", covariates), collapse = " + ")
  f_adj <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ ", rhs))
  fit_adj <- safe_fit(f_adj, morta)
  
  adj_df <- if (!is.null(fit_adj)) {
    tidy(fit_adj, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "bi_uv_protectionYes") %>%
      mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
             p = sprintf("%.4f", p.value)) %>%
      select(HR_CI, p)
  } else tibble(HR_CI = NA, p = NA)
  
  tibble(
    outcome         = outcome,
    crude_HR_Yes_vs_No = crude_df$HR_CI,
    crude_p         = crude_df$p,
    adj_HR_Yes_vs_No   = adj_df$HR_CI,
    adj_p           = adj_df$p
  )
}

final_tbl <- map_dfr(mortality_vars, run_one)
final_tbl

write.csv(final_tbl, "../result/5_4_2_cox_bi_uv_protection_crude_adjusted.csv", row.names = FALSE)

morta <- morta_backup
morta <- morta %>%
  mutate(
    bi_swr2008     = factor(bi_swr2008, levels = c(0,1), labels = c("No","Yes")),
    sex            = factor(sex, levels = c("female","male")),
    age            = factor(age, levels = c("37-47","48-54","55-59","60-63","64-73")),
    bmi            = factor(bmi, levels = c("underweight","healthy weight","overweight","obese")),
    sleep_score    = factor(sleep_score, levels = c("poor","intermediate","healthy")),
    drinking       = factor(drinking, levels = c("unhealthy","healthy")),
    phys_activity  = factor(phys_activity, levels = c("low","moderate","high")),
    highest_edu    = factor(highest_edu, levels = c("less than secondary","secondary school","university or above")),
    smoking_status = factor(smoking_status, levels = c("current","previous","never")),
    deprivation_quintile = factor(deprivation_quintile, levels = 1:5, ordered = TRUE)
  )

mortality_vars <- c(
  "all_cause_mortality",
  "cvd_mortality_narrow",
  "cancer_mortality",
  "skin_cancer_mortality",
  "melanoma_skin_cancer_mortality",
  "non_melanoma_mortality",
  "non_skin_cancer_mortality"
)

morta <- morta %>%
  mutate(across(all_of(mortality_vars), ~ ifelse(. == "yes", 1, 0)))

covariates <- c("sex","age","bmi","deprivation_quintile",
                "highest_edu","phys_activity","sleep_score",
                "smoking_status","drinking")

safe_fit <- function(formula, data) {
  tryCatch(coxph(formula, data = data, na.action = na.exclude), error = function(e) NULL)
}

run_one <- function(outcome) {
  f_crude <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ bi_swr2008"))
  fit_crude <- safe_fit(f_crude, morta)
  
  crude_df <- if (!is.null(fit_crude)) {
    tidy(fit_crude, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "bi_swr2008Yes") %>%
      mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
             p = sprintf("%.4f", p.value)) %>%
      select(HR_CI, p)
  } else tibble(HR_CI = NA, p = NA)
  
  rhs <- paste(c("bi_swr2008", covariates), collapse = " + ")
  f_adj <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ ", rhs))
  fit_adj <- safe_fit(f_adj, morta)
  
  adj_df <- if (!is.null(fit_adj)) {
    tidy(fit_adj, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "bi_swr2008Yes") %>%
      mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
             p = sprintf("%.4f", p.value)) %>%
      select(HR_CI, p)
  } else tibble(HR_CI = NA, p = NA)
  
  tibble(
    outcome         = outcome,
    crude_HR_Yes_vs_No = crude_df$HR_CI,
    crude_p         = crude_df$p,
    adj_HR_Yes_vs_No   = adj_df$HR_CI,
    adj_p           = adj_df$p
  )
}

final_tbl <- map_dfr(mortality_vars, run_one)
final_tbl

write.csv(final_tbl, "../result/5_4_3_cox_bi_swr2008_crude_adjusted.csv", row.names = FALSE)

table(morta$solarium_use)

morta <- morta_backup
morta <- morta %>%
  mutate(
    bi_solarium_use = factor(bi_solarium_use, levels = c(0,1), labels = c("No","Yes")),
    sex            = factor(sex, levels = c("female","male")),
    age            = factor(age, levels = c("37-47","48-54","55-59","60-63","64-73")),
    bmi            = factor(bmi, levels = c("underweight","healthy weight","overweight","obese")),
    sleep_score    = factor(sleep_score, levels = c("poor","intermediate","healthy")),
    drinking       = factor(drinking, levels = c("unhealthy","healthy")),
    phys_activity  = factor(phys_activity, levels = c("low","moderate","high")),
    highest_edu    = factor(highest_edu, levels = c("less than secondary","secondary school","university or above")),
    smoking_status = factor(smoking_status, levels = c("current","previous","never")),
    deprivation_quintile = factor(deprivation_quintile, levels = 1:5, ordered = TRUE)
  )

mortality_vars <- c(
  "all_cause_mortality",
  "cvd_mortality_narrow",
  "cancer_mortality",
  "skin_cancer_mortality",
  "melanoma_skin_cancer_mortality",
  "non_melanoma_mortality",
  "non_skin_cancer_mortality"
)

morta <- morta %>%
  mutate(across(all_of(mortality_vars), ~ ifelse(. == "yes", 1, 0)))

library(survival)
library(broom)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

covariates <- c("sex","age","bmi","deprivation_quintile",
                "highest_edu","phys_activity","sleep_score",
                "smoking_status","drinking")

safe_fit <- function(formula, data) {
  tryCatch(coxph(formula, data = data, na.action = na.exclude), error = function(e) NULL)
}

run_one <- function(outcome) {
  f_crude <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ bi_solarium_use"))
  fit_crude <- safe_fit(f_crude, morta)
  
  crude_df <- if (!is.null(fit_crude)) {
    tidy(fit_crude, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "bi_solarium_useYes") %>%
      mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
             p = sprintf("%.4f", p.value)) %>%
      select(HR_CI, p)
  } else tibble(HR_CI = NA, p = NA)
  
  rhs <- paste(c("bi_solarium_use", covariates), collapse = " + ")
  f_adj <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ ", rhs))
  fit_adj <- safe_fit(f_adj, morta)
  
  adj_df <- if (!is.null(fit_adj)) {
    tidy(fit_adj, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "bi_solarium_useYes") %>%
      mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
             p = sprintf("%.4f", p.value)) %>%
      select(HR_CI, p)
  } else tibble(HR_CI = NA, p = NA)
  
  tibble(
    outcome            = outcome,
    crude_HR_Yes_vs_No = crude_df$HR_CI,
    crude_p            = crude_df$p,
    adj_HR_Yes_vs_No   = adj_df$HR_CI,
    adj_p              = adj_df$p
  )
}

final_tbl <- map_dfr(mortality_vars, run_one)
final_tbl

write.csv(final_tbl, "../result/5_4_4_cox_bi_solarium_use_crude_adjusted.csv", row.names = FALSE)

table(morta$solarium_use)
table(morta$bi_solarium_use)
morta <- morta_backup

morta <- morta %>%
  mutate(
    bi_solarium_use2_num = ifelse(solarium_use == 0, 0,
                                  ifelse(!is.na(solarium_use), 1, NA)),
    bi_solarium_use2     = factor(bi_solarium_use2_num, levels = c(0, 1),
                                  labels = c("No", "Yes"))
  )

print(table(morta$solarium_use, morta$bi_solarium_use2, useNA = "ifany"))

morta <- morta %>%
  mutate(
    sex            = factor(sex, levels = c("female","male")),
    age            = factor(age, levels = c("37-47","48-54","55-59","60-63","64-73")),
    bmi            = factor(bmi, levels = c("underweight","healthy weight","overweight","obese")),
    sleep_score    = factor(sleep_score, levels = c("poor","intermediate","healthy")),
    drinking       = factor(drinking, levels = c("unhealthy","healthy")),
    phys_activity  = factor(phys_activity, levels = c("low","moderate","high")),
    highest_edu    = factor(highest_edu, levels = c("less than secondary","secondary school","university or above")),
    smoking_status = factor(smoking_status, levels = c("current","previous","never")),
    deprivation_quintile = factor(deprivation_quintile, levels = 1:5, ordered = TRUE)
  )

mortality_vars <- c(
  "all_cause_mortality",
  "cvd_mortality_narrow",
  "cancer_mortality",
  "skin_cancer_mortality",
  "melanoma_skin_cancer_mortality",
  "non_melanoma_mortality",
  "non_skin_cancer_mortality"
)

morta <- morta %>%
  mutate(across(all_of(mortality_vars), ~ ifelse(. == "yes", 1, 0)))

library(survival)
library(broom)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

covariates <- c("sex","age","bmi","deprivation_quintile",
                "highest_edu","phys_activity","sleep_score",
                "smoking_status","drinking")

safe_fit <- function(formula, data) {
  tryCatch(coxph(formula, data = data, na.action = na.exclude), error = function(e) NULL)
}

run_one <- function(outcome) {
  f_crude <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ bi_solarium_use2"))
  fit_crude <- safe_fit(f_crude, morta)
  
  crude_df <- if (!is.null(fit_crude)) {
    tidy(fit_crude, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "bi_solarium_use2Yes") %>%
      mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
             p = sprintf("%.4f", p.value)) %>%
      select(HR_CI, p)
  } else tibble(HR_CI = NA, p = NA)
  
  rhs <- paste(c("bi_solarium_use2", covariates), collapse = " + ")
  f_adj <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ ", rhs))
  fit_adj <- safe_fit(f_adj, morta)
  
  adj_df <- if (!is.null(fit_adj)) {
    tidy(fit_adj, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "bi_solarium_use2Yes") %>%
      mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
             p = sprintf("%.4f", p.value)) %>%
      select(HR_CI, p)
  } else tibble(HR_CI = NA, p = NA)
  
  tibble(
    outcome            = outcome,
    crude_HR_Yes_vs_No = crude_df$HR_CI,
    crude_p            = crude_df$p,
    adj_HR_Yes_vs_No   = adj_df$HR_CI,
    adj_p              = adj_df$p
  )
}

final_tbl <- map_dfr(mortality_vars, run_one)
print(final_tbl)

dir.create("../result", showWarnings = FALSE, recursive = TRUE)

write.csv(final_tbl, "../result/5_4_5_cox_bi_solarium_use2_crude_adjusted.csv", row.names = FALSE)
