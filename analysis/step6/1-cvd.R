if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(here)
library(ggplot2)
library(tidyverse)
library(here)
library(ggplot2)
library(gtsummary)
library(dplyr)
library(survival)
library(survival)
library(broom)
setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step6/"))

getwd()

cvd <- read.csv(here("washed_data", "4_2_3_ana_for_cvd_inpatient.csv"))

cvd_backup <- cvd
table(cvd$cvd_mortality_narrow)
table(cvd$diagnosed_after_assessment)

library(dplyr)
cvd <- cvd %>%
  filter(is.na(diagnosed_after_assessment) | diagnosed_after_assessment != "no")

endtime <- "2025-08-04" # most recent date code has been run
endtime <- as.Date(endtime)

cvd <- cvd %>%
  mutate(
    diag_time = ifelse(
      !is.na(cvd_mortality_narrow) & cvd_mortality_narrow == "yes",
      as.numeric(as.Date(diag_date) - as.Date(date_assess)),
      as.numeric(endtime) - as.numeric(as.Date(date_assess))
    )
  )

library(dplyr)
cvd <- cvd %>%
  mutate(
    cvd_narrow_bin = coalesce(as.integer(cvd_mortality_narrow == "yes"), 0L)
  )

table(cvd$diagnosed_after_assessment)
table(cvd$cvd_narrow_bin)

cvd <- cvd %>%
  dplyr::mutate(
    BEEM = factor(BEEM, levels = c("low", "medium", "high")),
    sex            = factor(sex, levels = c("female","male")),
    age            = factor(age, levels = c("37-47","48-54","55-59","60-63","64-73")),
    bmi            = factor(bmi, levels = c("underweight","healthy weight","overweight","obese")),
    sleep_score    = factor(sleep_score, levels = c("poor","intermediate","healthy")),
    drinking       = factor(drinking, levels = c("unhealthy","healthy")),
    phys_activity  = factor(phys_activity, levels = c("low","moderate","high")),
    highest_edu    = factor(highest_edu, levels = c("less than secondary","secondary school","university or above")),
    smoking_status = factor(smoking_status, levels = c("current","previous","never")),
    deprivation_quintile = factor(deprivation_quintile,
                                  levels = 1:5,
                                  ordered = TRUE) # 有序因子
  )

library(survival)
library(dplyr)
library(broom)

cvd <- cvd %>%
  mutate(
    BEEM_cat     = factor(BEEM, levels = c("low","medium","high"), ordered = FALSE),
    BEEM_numeric = as.numeric(factor(BEEM, levels = c("low","medium","high"))) - 1L
  )

covariates <- c("sex","age","bmi","deprivation_quintile",
                "highest_edu","phys_activity","sleep_score",
                "smoking_status","drinking")

safe_fit <- function(formula, data) {
  tryCatch(coxph(formula, data = data, na.action = na.exclude), error = function(e) NULL)
}

f_crude_cat <- Surv(diag_time, cvd_narrow_bin) ~ BEEM_cat
fit_crude_cat <- safe_fit(f_crude_cat, cvd)
crude_df <- if (!is.null(fit_crude_cat)) {
  tidy(fit_crude_cat, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c("BEEM_catmedium","BEEM_cathigh")) %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(term, HR_CI, p)
} else tibble(term=c("BEEM_catmedium","BEEM_cathigh"), HR_CI=NA, p=NA)

rhs <- paste(c("BEEM_cat", covariates), collapse = " + ")
f_adj_cat <- as.formula(paste0("Surv(diag_time, cvd_narrow_bin) ~ ", rhs))
fit_adj_cat <- safe_fit(f_adj_cat, cvd)
adj_df <- if (!is.null(fit_adj_cat)) {
  tidy(fit_adj_cat, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c("BEEM_catmedium","BEEM_cathigh")) %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(term, HR_CI, p)
} else tibble(term=c("BEEM_catmedium","BEEM_cathigh"), HR_CI=NA, p=NA)

fit_crude_trend <- safe_fit(Surv(diag_time, cvd_narrow_bin) ~ BEEM_numeric, cvd)
fit_adj_trend   <- safe_fit(as.formula(paste0(
  "Surv(diag_time, cvd_narrow_bin) ~ BEEM_numeric + ",
  paste(covariates, collapse = " + ")
)), cvd)

p_trend_crude <- if (!is.null(fit_crude_trend))
  sprintf("%.4f", tidy(fit_crude_trend)$p.value[ tidy(fit_crude_trend)$term == "BEEM_numeric" ]) else NA
p_trend_adj   <- if (!is.null(fit_adj_trend))
  sprintf("%.4f", tidy(fit_adj_trend)$p.value[ tidy(fit_adj_trend)$term == "BEEM_numeric" ]) else NA

summary_row <- tibble(
  outcome = "cvd_incidence_narrow",
  crude_HR_medium = crude_df$HR_CI[crude_df$term=="BEEM_catmedium"],
  crude_p_medium  = crude_df$p[crude_df$term=="BEEM_catmedium"],
  crude_HR_high   = crude_df$HR_CI[crude_df$term=="BEEM_cathigh"],
  crude_p_high    = crude_df$p[crude_df$term=="BEEM_cathigh"],
  p_trend_crude   = p_trend_crude,
  adj_HR_medium   = adj_df$HR_CI[adj_df$term=="BEEM_catmedium"],
  adj_p_medium    = adj_df$p[adj_df$term=="BEEM_catmedium"],
  adj_HR_high     = adj_df$HR_CI[adj_df$term=="BEEM_cathigh"],
  adj_p_high      = adj_df$p[adj_df$term=="BEEM_cathigh"],
  p_trend_adj     = p_trend_adj
)

print(summary_row)

write.csv(summary_row, "../result/6_1_cvd_cox_BEEM_crude_adjusted_ptrend.csv", row.names = FALSE)

N_total <- nrow(cvd)
N_event <- sum(cvd$cvd_narrow_bin, na.rm = TRUE)

event_summary <- tibble(
  outcome = "cvd_incidence_narrow",
  N_total = N_total,
  N_event = N_event
)

print(event_summary)

write.csv(event_summary, "../result/6_2_cvd_event_summary.csv", row.names = FALSE)
