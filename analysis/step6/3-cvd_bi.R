if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())
if (!nzchar(Sys.getenv("SUN_BEEM_CODE_ROOT"))) stop("Set SUN_BEEM_CODE_ROOT or use tools/run_analysis.R.")
source(file.path(Sys.getenv("SUN_BEEM_CODE_ROOT"), "analysis/shared/followup.R"))


library(tidyverse)
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(gtsummary)
library(dplyr)
library(survival)
library(survival)
library(broom)
setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step6/"))

getwd()

cvd <- read.csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "washed_data", "4_2_3_ana_for_cvd_inpatient.csv"))

cvd_backup <- cvd
table(cvd$cvd_mortality_narrow)
table(cvd$diagnosed_after_assessment)

library(dplyr)
cvd <- cvd %>%
  filter(is.na(diagnosed_after_assessment) | diagnosed_after_assessment != "no")

endtime <- registry_end_date("hospital")
endtime <- as.Date(endtime)

cvd <- apply_incidence_followup(cvd, "cvd_mortality_narrow", "diag_date", endtime, baseline_col = "date_assess")
cvd <- cvd %>%
  mutate(
    diag_time = diag_time
  )

library(dplyr)
cvd <- cvd %>%
  mutate(
    cvd_mortality_narrow_bin = coalesce(as.integer(cvd_mortality_narrow == "yes"), 0L)
  )

table(cvd$diagnosed_after_assessment)
table(cvd$cvd_mortality_narrow_bin)

cvd_for_bi <- cvd

library(survival)
library(dplyr)
library(broom)

cvd <- cvd_for_bi %>%
  mutate(
    bi_outdoor = as.integer(bi_outdoor),
    bi_outdoor_cat = factor(bi_outdoor, levels = c(0,1))
  )

covariates <- c("sex","age","bmi","deprivation_quintile",
                "highest_edu","phys_activity","sleep_score",
                "smoking_status","drinking")

safe_fit <- function(formula, data) {
  tryCatch(coxph(formula, data = data, na.action = na.exclude), error = function(e) NULL)
}

fit_crude <- safe_fit(Surv(diag_time, cvd_mortality_narrow_bin) ~ bi_outdoor_cat, cvd)
crude_df <- if (!is.null(fit_crude)) {
  tidy(fit_crude, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "bi_outdoor_cat1") %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(HR_CI, p)
} else tibble(HR_CI = NA, p = NA)

rhs <- paste(c("bi_outdoor_cat", covariates), collapse = " + ")
fit_adj <- safe_fit(as.formula(
  paste0("Surv(diag_time, cvd_mortality_narrow_bin) ~ ", rhs)
), cvd)
adj_df <- if (!is.null(fit_adj)) {
  tidy(fit_adj, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "bi_outdoor_cat1") %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(HR_CI, p)
} else tibble(HR_CI = NA, p = NA)

summary_row <- tibble(
  outcome = "cvd_mortality_narrow",
  crude_HR = crude_df$HR_CI,
  crude_p  = crude_df$p,
  adj_HR   = adj_df$HR_CI,
  adj_p    = adj_df$p
)

print(summary_row)

write.csv(summary_row, "../result/6_3_1_cvd_cox_bi_outdoor_crude_adjusted.csv", row.names = FALSE)

library(survival)
library(dplyr)
library(broom)

cvd <- cvd_for_bi %>%
  mutate(
    bi_uv_protection = as.integer(bi_uv_protection),
    bi_uv_protection_cat = factor(bi_uv_protection, levels = c(0,1))
  )

covariates <- c("sex","age","bmi","deprivation_quintile",
                "highest_edu","phys_activity","sleep_score",
                "smoking_status","drinking")

safe_fit <- function(formula, data) {
  tryCatch(coxph(formula, data = data, na.action = na.exclude), error = function(e) NULL)
}

fit_crude <- safe_fit(Surv(diag_time, cvd_mortality_narrow_bin) ~ bi_uv_protection_cat, cvd)
crude_df <- if (!is.null(fit_crude)) {
  tidy(fit_crude, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "bi_uv_protection_cat1") %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(HR_CI, p)
} else tibble(HR_CI = NA, p = NA)

rhs <- paste(c("bi_uv_protection_cat", covariates), collapse = " + ")
fit_adj <- safe_fit(as.formula(
  paste0("Surv(diag_time, cvd_mortality_narrow_bin) ~ ", rhs)
), cvd)
adj_df <- if (!is.null(fit_adj)) {
  tidy(fit_adj, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "bi_uv_protection_cat1") %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(HR_CI, p)
} else tibble(HR_CI = NA, p = NA)

summary_row <- tibble(
  outcome = "cvd_mortality_narrow",
  crude_HR = crude_df$HR_CI,
  crude_p  = crude_df$p,
  adj_HR   = adj_df$HR_CI,
  adj_p    = adj_df$p
)

print(summary_row)

write.csv(summary_row, "../result/6_3_2_cvd_cox_bi_uv_protection_crude_adjusted.csv", row.names = FALSE)

library(survival)
library(dplyr)
library(broom)

cvd <- cvd_for_bi %>%
  mutate(
    bi_swr2008 = as.integer(bi_swr2008),
    bi_swr2008_cat = factor(bi_swr2008, levels = c(0,1))
  )

covariates <- c("sex","age","bmi","deprivation_quintile",
                "highest_edu","phys_activity","sleep_score",
                "smoking_status","drinking")

safe_fit <- function(formula, data) {
  tryCatch(coxph(formula, data = data, na.action = na.exclude), error = function(e) NULL)
}

fit_crude <- safe_fit(Surv(diag_time, cvd_mortality_narrow_bin) ~ bi_swr2008_cat, cvd)
crude_df <- if (!is.null(fit_crude)) {
  tidy(fit_crude, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "bi_swr2008_cat1") %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(HR_CI, p)
} else tibble(HR_CI = NA, p = NA)

rhs <- paste(c("bi_swr2008_cat", covariates), collapse = " + ")
fit_adj <- safe_fit(as.formula(
  paste0("Surv(diag_time, cvd_mortality_narrow_bin) ~ ", rhs)
), cvd)
adj_df <- if (!is.null(fit_adj)) {
  tidy(fit_adj, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "bi_swr2008_cat1") %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(HR_CI, p)
} else tibble(HR_CI = NA, p = NA)

summary_row <- tibble(
  outcome = "cvd_mortality_narrow",
  crude_HR = crude_df$HR_CI,
  crude_p  = crude_df$p,
  adj_HR   = adj_df$HR_CI,
  adj_p    = adj_df$p
)

print(summary_row)

write.csv(summary_row, "../result/6_3_3_cvd_cox_bi_swr2008_crude_adjusted.csv", row.names = FALSE)

library(survival)
library(dplyr)
library(broom)

cvd <- cvd_for_bi %>%
  mutate(
    bi_solarium_use = as.integer(bi_solarium_use),
    bi_solarium_use_cat = factor(bi_solarium_use, levels = c(0,1))
  )

covariates <- c("sex","age","bmi","deprivation_quintile",
                "highest_edu","phys_activity","sleep_score",
                "smoking_status","drinking")

safe_fit <- function(formula, data) {
  tryCatch(coxph(formula, data = data, na.action = na.exclude), error = function(e) NULL)
}

fit_crude <- safe_fit(Surv(diag_time, cvd_mortality_narrow_bin) ~ bi_solarium_use_cat, cvd)
crude_df <- if (!is.null(fit_crude)) {
  tidy(fit_crude, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "bi_solarium_use_cat1") %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(HR_CI, p)
} else tibble(HR_CI = NA, p = NA)

rhs <- paste(c("bi_solarium_use_cat", covariates), collapse = " + ")
fit_adj <- safe_fit(as.formula(
  paste0("Surv(diag_time, cvd_mortality_narrow_bin) ~ ", rhs)
), cvd)
adj_df <- if (!is.null(fit_adj)) {
  tidy(fit_adj, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "bi_solarium_use_cat1") %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(HR_CI, p)
} else tibble(HR_CI = NA, p = NA)

summary_row <- tibble(
  outcome = "cvd_mortality_narrow",
  crude_HR = crude_df$HR_CI,
  crude_p  = crude_df$p,
  adj_HR   = adj_df$HR_CI,
  adj_p    = adj_df$p
)

print(summary_row)

write.csv(summary_row, "../result/6_3_4_cvd_cox_bi_solarium_use_crude_adjusted.csv", row.names = FALSE)
