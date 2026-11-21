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
setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step7/"))

getwd()
library(tidyverse)
library(survival)
library(broom)
library(rlang)
endtime <- registry_end_date("cancer")

outcomes_cfg <- tribble(
  ~in_csv,                                 ~status_var,               ~after_var,
  "4_2_4_ana_for_allcancer_reg.csv",       "cancer_status",           "diagnosed_after_assessment",
  "4_2_5_ana_for_skin_cancer_reg.csv",     "skin_cancer_status",      "skin_cancer_after_assessment",
  "4_2_6_ana_for_melanoma_reg.csv",        "melanoma_status",         "melanoma_after_assessment",
  "4_2_7_ana_for_other_skin_cancer_reg.csv","other_skin_cancer_status","other_skin_cancer_after_assessment",
  "4_2_8_ana_for_non_skin_cancer_reg.csv", "non_skin_cancer_status",  "non_skin_cancer_after_assessment"
)

fit_one <- function(in_csv, status_var, after_var, exposure_var,
                    date_event = "cancer_date", date_assess = "date_assess") {
  
  status_sym <- sym(status_var)
  after_sym  <- sym(after_var)
  ev_sym     <- sym(date_event)
  as_sym     <- sym(date_assess)
  expo_sym   <- sym(exposure_var)
  
  dat <- read.csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "washed_data", in_csv))
  
  dat <- dat %>% filter(is.na(!!after_sym) | !!after_sym != "no")
  
  dat <- apply_incidence_followup(dat, status_var, date_event, endtime, baseline_col = date_assess)
dat <- dat %>%
    mutate(
      diag_time = diag_time
    )
  
  event_bin <- paste0(status_var, "_bin")
  dat <- dat %>%
    mutate(!!event_bin := coalesce(as.integer(!!status_sym == "yes"), 0L))
  
  dat <- dat %>%
    mutate(
      BEEM = factor(BEEM, levels = c("low","medium","high")),
      sex            = factor(sex, levels = c("female","male")),
      age            = factor(age, levels = c("37-47","48-54","55-59","60-63","64-73")),
      bmi            = factor(bmi, levels = c("underweight","healthy weight","overweight","obese")),
      sleep_score    = factor(sleep_score, levels = c("poor","intermediate","healthy")),
      drinking       = factor(drinking, levels = c("unhealthy","healthy")),
      phys_activity  = factor(phys_activity, levels = c("low","moderate","high")),
      highest_edu    = factor(highest_edu, levels = c("less than secondary","secondary school","university or above")),
      smoking_status = factor(smoking_status, levels = c("current","previous","never")),
      deprivation_quintile = factor(deprivation_quintile, levels = 1:5, ordered = TRUE),
      !!expo_sym := as.integer(!!expo_sym),
      expo_cat  = factor(!!expo_sym, levels = c(0,1))   # 0=ref
    )
  
  covariates <- c("sex","age","bmi","deprivation_quintile",
                  "highest_edu","phys_activity","sleep_score",
                  "smoking_status","drinking")
  
  safe_fit <- function(formula, data) tryCatch(coxph(formula, data = data, na.action = na.exclude), error = function(e) NULL)
  
  f_crude <- as.formula(paste0("Surv(diag_time, ", event_bin, ") ~ expo_cat"))
  fit_crude <- safe_fit(f_crude, dat)
  crude_df <- if (!is.null(fit_crude)) {
    tidy(fit_crude, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "expo_cat1") %>%
      transmute(
        crude_HR = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
        crude_p  = sprintf("%.4f", p.value)
      )
  } else tibble(crude_HR = NA_character_, crude_p = NA_character_)
  
  rhs <- paste(c("expo_cat", covariates), collapse = " + ")
  f_adj <- as.formula(paste0("Surv(diag_time, ", event_bin, ") ~ ", rhs))
  fit_adj <- safe_fit(f_adj, dat)
  adj_df <- if (!is.null(fit_adj)) {
    tidy(fit_adj, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "expo_cat1") %>%
      transmute(
        adj_HR = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
        adj_p  = sprintf("%.4f", p.value)
      )
  } else tibble(adj_HR = NA_character_, adj_p = NA_character_)
  
  tibble(
    exposure = exposure_var,
    outcome  = status_var,
    crude_HR = crude_df$crude_HR,
    crude_p  = crude_df$crude_p,
    adj_HR   = adj_df$adj_HR,
    adj_p    = adj_df$adj_p
  )
}

run_one_exposure <- function(exposure_var, out_csv_path) {
  res <- purrr::pmap_dfr(
    outcomes_cfg,
    ~ fit_one(in_csv = ..1, status_var = ..2, after_var = ..3, exposure_var = exposure_var)
  )
  print(res)
  write.csv(res, out_csv_path, row.names = FALSE)
  invisible(res)
}

res_outdoor <- run_one_exposure(
  "bi_outdoor",
  "../result/7_3_1_bi_outdoor_cox_all_outcomes_crude_adjusted.csv"
)

res_uvprot <- run_one_exposure(
  "bi_uv_protection",
  "../result/7_3_2_bi_uv_protection_cox_all_outcomes_crude_adjusted.csv"
)

res_swr2008 <- run_one_exposure(
  "bi_swr2008",
  "../result/7_3_3_bi_swr2008_cox_all_outcomes_crude_adjusted.csv"
)

res_solarium <- run_one_exposure(
  "bi_solarium_use",
  "../result/7_3_4_bi_solarium_use_cox_all_outcomes_crude_adjusted.csv"
)
