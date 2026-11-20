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
setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step7/"))

getwd()

save_event_summary <- function(data, event_col, outcome_name, out_file) {
  N_total <- nrow(data)
  N_event <- sum(data[[event_col]], na.rm = TRUE)
  event_summary <- tibble(
    outcome = outcome_name,
    N_total = N_total,
    N_event = N_event
  )
  print(event_summary)
  write.csv(event_summary, out_file, row.names = FALSE)
}

cancer <- read.csv(here("washed_data", "4_2_4_ana_for_allcancer_reg.csv"))

cancer_backup <- cancer
table(cancer$cancer_status)
table(cancer$diagnosed_after_assessment)

library(dplyr)
cancer <- cancer %>%
  filter(is.na(diagnosed_after_assessment) | diagnosed_after_assessment != "no")

endtime <- "2025-08-04" # most recent date code has been run
endtime <- as.Date(endtime)

cancer <- cancer %>%
  mutate(
    diag_time = ifelse(
      !is.na(cancer_status) & cancer_status == "yes",
      as.numeric(as.Date(cancer_date) - as.Date(date_assess)),
      as.numeric(endtime) - as.numeric(as.Date(date_assess))
    )
  )

library(dplyr)
cancer <- cancer %>%
  mutate(
    cancer_status_bin = coalesce(as.integer(cancer_status == "yes"), 0L)
  )

table(cancer$diagnosed_after_assessment)
table(cancer$cancer_status_bin)

cancer <- cancer %>%
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

cancer <- cancer %>%
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

f_crude_cat <- Surv(diag_time, cancer_status_bin) ~ BEEM_cat
fit_crude_cat <- safe_fit(f_crude_cat, cancer)
crude_df <- if (!is.null(fit_crude_cat)) {
  tidy(fit_crude_cat, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c("BEEM_catmedium","BEEM_cathigh")) %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(term, HR_CI, p)
} else tibble(term=c("BEEM_catmedium","BEEM_cathigh"), HR_CI=NA, p=NA)

rhs <- paste(c("BEEM_cat", covariates), collapse = " + ")
f_adj_cat <- as.formula(paste0("Surv(diag_time, cancer_status_bin) ~ ", rhs))
fit_adj_cat <- safe_fit(f_adj_cat, cancer)
adj_df <- if (!is.null(fit_adj_cat)) {
  tidy(fit_adj_cat, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c("BEEM_catmedium","BEEM_cathigh")) %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(term, HR_CI, p)
} else tibble(term=c("BEEM_catmedium","BEEM_cathigh"), HR_CI=NA, p=NA)

fit_crude_trend <- safe_fit(Surv(diag_time, cancer_status_bin) ~ BEEM_numeric, cancer)
fit_adj_trend   <- safe_fit(as.formula(paste0(
  "Surv(diag_time, cancer_status_bin) ~ BEEM_numeric + ",
  paste(covariates, collapse = " + ")
)), cancer)

p_trend_crude <- if (!is.null(fit_crude_trend))
  sprintf("%.4f", tidy(fit_crude_trend)$p.value[ tidy(fit_crude_trend)$term == "BEEM_numeric" ]) else NA
p_trend_adj   <- if (!is.null(fit_adj_trend))
  sprintf("%.4f", tidy(fit_adj_trend)$p.value[ tidy(fit_adj_trend)$term == "BEEM_numeric" ]) else NA

summary_row <- tibble(
  outcome = "cancer_status",
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

write.csv(summary_row, "../result/7_1_1_cancerall_cox_BEEM_crude_adjusted_ptrend.csv", row.names = FALSE)
save_event_summary(
  data = cancer,
  event_col = "cancer_status_bin",
  outcome_name = "all_cancer_incidence",
  out_file = "../result/7_1_1_cancerall_event_summary.csv"
)

cancer <- read.csv(here("washed_data", "4_2_5_ana_for_skin_cancer_reg.csv"))

cancer_backup <- cancer
table(cancer$skin_cancer_status)
table(cancer$skin_cancer_after_assessment)

library(dplyr)
cancer <- cancer %>%
  filter(is.na(skin_cancer_after_assessment) | skin_cancer_after_assessment != "no")

endtime <- "2025-08-04" # most recent date code has been run
endtime <- as.Date(endtime)

cancer <- cancer %>%
  mutate(
    diag_time = ifelse(
      !is.na(skin_cancer_status) & skin_cancer_status == "yes",
      as.numeric(as.Date(cancer_date) - as.Date(date_assess)),
      as.numeric(endtime) - as.numeric(as.Date(date_assess))
    )
  )

library(dplyr)
cancer <- cancer %>%
  mutate(
    skin_cancer_status_bin = coalesce(as.integer(skin_cancer_status == "yes"), 0L)
  )

table(cancer$skin_cancer_after_assessment)
table(cancer$skin_cancer_status)

cancer <- cancer %>%
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

cancer <- cancer %>%
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

f_crude_cat <- Surv(diag_time, skin_cancer_status_bin) ~ BEEM_cat
fit_crude_cat <- safe_fit(f_crude_cat, cancer)
crude_df <- if (!is.null(fit_crude_cat)) {
  tidy(fit_crude_cat, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c("BEEM_catmedium","BEEM_cathigh")) %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(term, HR_CI, p)
} else tibble(term=c("BEEM_catmedium","BEEM_cathigh"), HR_CI=NA, p=NA)

rhs <- paste(c("BEEM_cat", covariates), collapse = " + ")
f_adj_cat <- as.formula(paste0("Surv(diag_time, skin_cancer_status_bin) ~ ", rhs))
fit_adj_cat <- safe_fit(f_adj_cat, cancer)
adj_df <- if (!is.null(fit_adj_cat)) {
  tidy(fit_adj_cat, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c("BEEM_catmedium","BEEM_cathigh")) %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(term, HR_CI, p)
} else tibble(term=c("BEEM_catmedium","BEEM_cathigh"), HR_CI=NA, p=NA)

fit_crude_trend <- safe_fit(Surv(diag_time, skin_cancer_status_bin) ~ BEEM_numeric, cancer)
fit_adj_trend   <- safe_fit(as.formula(paste0(
  "Surv(diag_time, skin_cancer_status_bin) ~ BEEM_numeric + ",
  paste(covariates, collapse = " + ")
)), cancer)

p_trend_crude <- if (!is.null(fit_crude_trend))
  sprintf("%.4f", tidy(fit_crude_trend)$p.value[ tidy(fit_crude_trend)$term == "BEEM_numeric" ]) else NA
p_trend_adj   <- if (!is.null(fit_adj_trend))
  sprintf("%.4f", tidy(fit_adj_trend)$p.value[ tidy(fit_adj_trend)$term == "BEEM_numeric" ]) else NA

summary_row <- tibble(
  outcome = "skin_cancer_status",
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

write.csv(summary_row, "../result/7_1_2_skin_cancer_cox_BEEM_crude_adjusted_ptrend.csv", row.names = FALSE)
save_event_summary(
  data = cancer,
  event_col = "skin_cancer_status_bin",
  outcome_name = "skin_cancer_incidence",
  out_file = "../result/7_1_2_skin_cancer_event_summary.csv"
)

cancer <- read.csv(here("washed_data", "4_2_6_ana_for_melanoma_reg.csv"))

cancer_backup <- cancer
table(cancer$melanoma_status)
table(cancer$melanoma_after_assessment)

library(dplyr)
cancer <- cancer %>%
  filter(is.na(melanoma_after_assessment) | melanoma_after_assessment != "no")

endtime <- "2025-08-04" # most recent date code has been run
endtime <- as.Date(endtime)

cancer <- cancer %>%
  mutate(
    diag_time = ifelse(
      !is.na(melanoma_status) & melanoma_status == "yes",
      as.numeric(as.Date(cancer_date) - as.Date(date_assess)),
      as.numeric(endtime) - as.numeric(as.Date(date_assess))
    )
  )

library(dplyr)
cancer <- cancer %>%
  mutate(
    melanoma_status_bin = coalesce(as.integer(melanoma_status == "yes"), 0L)
  )

table(cancer$melanoma_after_assessment)
table(cancer$melanoma_status)

cancer <- cancer %>%
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

cancer <- cancer %>%
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

f_crude_cat <- Surv(diag_time, melanoma_status_bin) ~ BEEM_cat
fit_crude_cat <- safe_fit(f_crude_cat, cancer)
crude_df <- if (!is.null(fit_crude_cat)) {
  tidy(fit_crude_cat, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c("BEEM_catmedium","BEEM_cathigh")) %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(term, HR_CI, p)
} else tibble(term=c("BEEM_catmedium","BEEM_cathigh"), HR_CI=NA, p=NA)

rhs <- paste(c("BEEM_cat", covariates), collapse = " + ")
f_adj_cat <- as.formula(paste0("Surv(diag_time, melanoma_status_bin) ~ ", rhs))
fit_adj_cat <- safe_fit(f_adj_cat, cancer)
adj_df <- if (!is.null(fit_adj_cat)) {
  tidy(fit_adj_cat, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c("BEEM_catmedium","BEEM_cathigh")) %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(term, HR_CI, p)
} else tibble(term=c("BEEM_catmedium","BEEM_cathigh"), HR_CI=NA, p=NA)

fit_crude_trend <- safe_fit(Surv(diag_time, melanoma_status_bin) ~ BEEM_numeric, cancer)
fit_adj_trend   <- safe_fit(as.formula(paste0(
  "Surv(diag_time, melanoma_status_bin) ~ BEEM_numeric + ",
  paste(covariates, collapse = " + ")
)), cancer)

p_trend_crude <- if (!is.null(fit_crude_trend))
  sprintf("%.4f", tidy(fit_crude_trend)$p.value[ tidy(fit_crude_trend)$term == "BEEM_numeric" ]) else NA
p_trend_adj   <- if (!is.null(fit_adj_trend))
  sprintf("%.4f", tidy(fit_adj_trend)$p.value[ tidy(fit_adj_trend)$term == "BEEM_numeric" ]) else NA

summary_row <- tibble(
  outcome = "melanoma_status",
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

write.csv(summary_row, "../result/7_1_3_melanoma_cox_BEEM_crude_adjusted_ptrend.csv", row.names = FALSE)

save_event_summary(
  data = cancer,
  event_col = "melanoma_status_bin",
  outcome_name = "melanoma_incidence",
  out_file = "../result/7_1_3_melanoma_event_summary.csv"
)

cancer <- read.csv(here("washed_data", "4_2_7_ana_for_other_skin_cancer_reg.csv"))

cancer_backup <- cancer
table(cancer$other_skin_cancer_status)
table(cancer$other_skin_cancer_after_assessment)

library(dplyr)
cancer <- cancer %>%
  filter(is.na(other_skin_cancer_after_assessment) | other_skin_cancer_after_assessment != "no")

endtime <- "2025-08-04" # most recent date code has been run
endtime <- as.Date(endtime)

cancer <- cancer %>%
  mutate(
    diag_time = ifelse(
      !is.na(other_skin_cancer_status) & other_skin_cancer_status == "yes",
      as.numeric(as.Date(cancer_date) - as.Date(date_assess)),
      as.numeric(endtime) - as.numeric(as.Date(date_assess))
    )
  )

library(dplyr)
cancer <- cancer %>%
  mutate(
    other_skin_cancer_status_bin = coalesce(as.integer(other_skin_cancer_status == "yes"), 0L)
  )

table(cancer$other_skin_cancer_after_assessment)
table(cancer$other_skin_cancer_status)

cancer <- cancer %>%
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

cancer <- cancer %>%
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

f_crude_cat <- Surv(diag_time, other_skin_cancer_status_bin) ~ BEEM_cat
fit_crude_cat <- safe_fit(f_crude_cat, cancer)
crude_df <- if (!is.null(fit_crude_cat)) {
  tidy(fit_crude_cat, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c("BEEM_catmedium","BEEM_cathigh")) %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(term, HR_CI, p)
} else tibble(term=c("BEEM_catmedium","BEEM_cathigh"), HR_CI=NA, p=NA)

rhs <- paste(c("BEEM_cat", covariates), collapse = " + ")
f_adj_cat <- as.formula(paste0("Surv(diag_time, other_skin_cancer_status_bin) ~ ", rhs))
fit_adj_cat <- safe_fit(f_adj_cat, cancer)
adj_df <- if (!is.null(fit_adj_cat)) {
  tidy(fit_adj_cat, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c("BEEM_catmedium","BEEM_cathigh")) %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(term, HR_CI, p)
} else tibble(term=c("BEEM_catmedium","BEEM_cathigh"), HR_CI=NA, p=NA)

fit_crude_trend <- safe_fit(Surv(diag_time, other_skin_cancer_status_bin) ~ BEEM_numeric, cancer)
fit_adj_trend   <- safe_fit(as.formula(paste0(
  "Surv(diag_time, other_skin_cancer_status_bin) ~ BEEM_numeric + ",
  paste(covariates, collapse = " + ")
)), cancer)

p_trend_crude <- if (!is.null(fit_crude_trend))
  sprintf("%.4f", tidy(fit_crude_trend)$p.value[ tidy(fit_crude_trend)$term == "BEEM_numeric" ]) else NA
p_trend_adj   <- if (!is.null(fit_adj_trend))
  sprintf("%.4f", tidy(fit_adj_trend)$p.value[ tidy(fit_adj_trend)$term == "BEEM_numeric" ]) else NA

summary_row <- tibble(
  outcome = "other_skin_cancer_status",
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

write.csv(summary_row, "../result/7_1_4_other_skin_cancer_cox_BEEM_crude_adjusted_ptrend.csv", row.names = FALSE)

save_event_summary(
  data = cancer,
  event_col = "other_skin_cancer_status_bin",
  outcome_name = "non_melanoma_skin_cancer_incidence",
  out_file = "../result/7_1_4_other_skin_cancer_event_summary.csv"
)

cancer <- read.csv(here("washed_data", "4_2_8_ana_for_non_skin_cancer_reg.csv"))

cancer_backup <- cancer
table(cancer$non_skin_cancer_status)
table(cancer$non_skin_cancer_after_assessment)

library(dplyr)
cancer <- cancer %>%
  filter(is.na(non_skin_cancer_after_assessment) | non_skin_cancer_after_assessment != "no")

endtime <- "2025-08-04" # most recent date code has been run
endtime <- as.Date(endtime)

cancer <- cancer %>%
  mutate(
    diag_time = ifelse(
      !is.na(non_skin_cancer_status) & non_skin_cancer_status == "yes",
      as.numeric(as.Date(cancer_date) - as.Date(date_assess)),
      as.numeric(endtime) - as.numeric(as.Date(date_assess))
    )
  )

library(dplyr)
cancer <- cancer %>%
  mutate(
    non_skin_cancer_status_bin = coalesce(as.integer(non_skin_cancer_status == "yes"), 0L)
  )

table(cancer$non_skin_cancer_after_assessment)
table(cancer$non_skin_cancer_status)

cancer <- cancer %>%
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

cancer <- cancer %>%
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

f_crude_cat <- Surv(diag_time, non_skin_cancer_status_bin) ~ BEEM_cat
fit_crude_cat <- safe_fit(f_crude_cat, cancer)
crude_df <- if (!is.null(fit_crude_cat)) {
  tidy(fit_crude_cat, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c("BEEM_catmedium","BEEM_cathigh")) %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(term, HR_CI, p)
} else tibble(term=c("BEEM_catmedium","BEEM_cathigh"), HR_CI=NA, p=NA)

rhs <- paste(c("BEEM_cat", covariates), collapse = " + ")
f_adj_cat <- as.formula(paste0("Surv(diag_time, non_skin_cancer_status_bin) ~ ", rhs))
fit_adj_cat <- safe_fit(f_adj_cat, cancer)
adj_df <- if (!is.null(fit_adj_cat)) {
  tidy(fit_adj_cat, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c("BEEM_catmedium","BEEM_cathigh")) %>%
    mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value)) %>%
    select(term, HR_CI, p)
} else tibble(term=c("BEEM_catmedium","BEEM_cathigh"), HR_CI=NA, p=NA)

fit_crude_trend <- safe_fit(Surv(diag_time, non_skin_cancer_status_bin) ~ BEEM_numeric, cancer)
fit_adj_trend   <- safe_fit(as.formula(paste0(
  "Surv(diag_time, non_skin_cancer_status_bin) ~ BEEM_numeric + ",
  paste(covariates, collapse = " + ")
)), cancer)

p_trend_crude <- if (!is.null(fit_crude_trend))
  sprintf("%.4f", tidy(fit_crude_trend)$p.value[ tidy(fit_crude_trend)$term == "BEEM_numeric" ]) else NA
p_trend_adj   <- if (!is.null(fit_adj_trend))
  sprintf("%.4f", tidy(fit_adj_trend)$p.value[ tidy(fit_adj_trend)$term == "BEEM_numeric" ]) else NA

summary_row <- tibble(
  outcome = "non_skin_cancer_status",
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

write.csv(summary_row, "../result/7_1_5_non_skin_cancer_cox_BEEM_crude_adjusted_ptrend.csv", row.names = FALSE)
save_event_summary(
  data = cancer,
  event_col = "non_skin_cancer_status_bin",
  outcome_name = "non_skin_cancer_incidence",
  out_file = "../result/7_1_5_non_skin_cancer_event_summary.csv"
)
