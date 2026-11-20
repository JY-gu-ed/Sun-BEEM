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
setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step5/"))

getwd()

morta <- read.csv(here("washed_data", "4_2_2_ana_for_mortality.csv"))

morta_backup<- morta 

morta <- morta %>%
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

morta <- morta %>%
  mutate(
    BEEM = factor(BEEM, levels = c("low","medium","high")),
    BEEM_numeric = as.numeric(BEEM)  # low=1, medium=2, high=3
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

covariates <- c("sex","age","bmi","deprivation_quintile",
                "highest_edu","phys_activity","sleep_score",
                "smoking_status","drinking")

safe_fit <- function(formula, data) {
  tryCatch(coxph(formula, data = data, na.action = na.exclude), error = function(e) NULL)
}

run_one <- function(outcome) {
  f_crude_cat <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ BEEM"))
  fit_crude_cat <- safe_fit(f_crude_cat, morta)
  
  crude_df <- if (!is.null(fit_crude_cat)) {
    tidy(fit_crude_cat, exponentiate = TRUE, conf.int = TRUE) |>
      filter(term %in% c("BEEMmedium","BEEMhigh")) |>
      mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
             p = sprintf("%.4f", p.value)) |>
      select(term, HR_CI, p)
  } else tibble(term=c("BEEMmedium","BEEMhigh"), HR_CI=NA, p=NA)
  
  rhs <- paste(c("BEEM", covariates), collapse = " + ")
  f_adj_cat <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ ", rhs))
  fit_adj_cat <- safe_fit(f_adj_cat, morta)
  
  adj_df <- if (!is.null(fit_adj_cat)) {
    tidy(fit_adj_cat, exponentiate = TRUE, conf.int = TRUE) |>
      filter(term %in% c("BEEMmedium","BEEMhigh")) |>
      mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
             p = sprintf("%.4f", p.value)) |>
      select(term, HR_CI, p)
  } else tibble(term=c("BEEMmedium","BEEMhigh"), HR_CI=NA, p=NA)
  
  f_crude_trend <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ BEEM_numeric"))
  f_adj_trend   <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ ",
                                     paste(c("BEEM_numeric", covariates), collapse = " + ")))
  fit_crude_trend <- safe_fit(f_crude_trend, morta)
  fit_adj_trend   <- safe_fit(f_adj_trend,   morta)
  
  p_trend_crude <- if (!is.null(fit_crude_trend))
    sprintf("%.4f", tidy(fit_crude_trend)$p.value[ tidy(fit_crude_trend)$term == "BEEM_numeric" ]) else NA
  p_trend_adj   <- if (!is.null(fit_adj_trend))
    sprintf("%.4f", tidy(fit_adj_trend)$p.value[ tidy(fit_adj_trend)$term == "BEEM_numeric" ]) else NA
  
  tibble(
    outcome = outcome,
    crude_HR_medium = crude_df$HR_CI[crude_df$term=="BEEMmedium"],
    crude_p_medium  = crude_df$p[crude_df$term=="BEEMmedium"],
    crude_HR_high   = crude_df$HR_CI[crude_df$term=="BEEMhigh"],
    crude_p_high    = crude_df$p[crude_df$term=="BEEMhigh"],
    p_trend_crude   = p_trend_crude,
    adj_HR_medium   = adj_df$HR_CI[adj_df$term=="BEEMmedium"],
    adj_p_medium    = adj_df$p[adj_df$term=="BEEMmedium"],
    adj_HR_high     = adj_df$HR_CI[adj_df$term=="BEEMhigh"],
    adj_p_high      = adj_df$p[adj_df$term=="BEEMhigh"],
    p_trend_adj     = p_trend_adj
  )
}

final_tbl <- map_dfr(mortality_vars, run_one)

write.csv(final_tbl, "../result/5_2_1_cox_BEEM_crude_adjusted_ptrend.csv", row.names = FALSE)
final_tbl

n_cases_tbl <- morta %>%
  summarise(across(all_of(mortality_vars), ~ sum(. == 1, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(),
               names_to = "outcome",
               values_to = "No_cases")

write.csv(n_cases_tbl,
          "../result/5_2_2_no_of_cases.csv",
          row.names = FALSE)

n_cases_tbl

final_tbl <- read.csv("../result/5_2_1_cox_BEEM_crude_adjusted_ptrend.csv", stringsAsFactors = FALSE)

head(final_tbl)

library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(ggplot2)
library(scales)

final_tbl <- read_csv("../result/5_2_1_cox_BEEM_crude_adjusted_ptrend.csv")

parse_hr <- function(x){
  as.numeric(str_extract(x, "^[0-9]+\\.?[0-9]*"))
}

final_tbl$outcome
hr_tbl <- final_tbl %>%
  transmute(
    outcome,
    RR_med_crude = parse_hr(crude_HR_medium),
    RR_hi_crude  = parse_hr(crude_HR_high),
    RR_med_adj   = parse_hr(adj_HR_medium),
    RR_hi_adj    = parse_hr(adj_HR_high)
  )

p_tbl <- morta %>%
  filter(!is.na(BEEM)) %>%
  count(BEEM) %>%
  mutate(p = n / sum(n)) %>%
  select(BEEM, p) %>%
  pivot_wider(names_from = BEEM, values_from = p, names_prefix = "p_") %>%
  mutate(across(everything(), ~replace_na(., 0)))  # p_low, p_medium, p_high

p_low   <- p_tbl$p_low
p_med   <- p_tbl$p_medium
p_high  <- p_tbl$p_high

calc_paf_A <- function(p_low, p_med, p_high, RR_med, RR_hi){
  denom <- p_low*1 + p_med*RR_med + p_high*RR_hi
  if (is.na(denom) || denom == 0) return(NA_real_)
  (denom - 1) / denom
}

calc_paf_B <- function(p_low, p_med, p_high, RR_med, RR_hi){
  denom <- p_low*1 + p_med*RR_med + p_high*RR_hi
  if (is.na(denom) || denom == 0) return(NA_real_)
  (denom - RR_hi) / denom
}

paf_tbl <- hr_tbl %>%
  rowwise() %>%
  mutate(
    PAF_A_crude = calc_paf_A(p_low, p_med, p_high, RR_med_crude, RR_hi_crude),
    PAF_A_adj   = calc_paf_A(p_low, p_med, p_high, RR_med_adj,   RR_hi_adj),
    
    PAF_B_crude = calc_paf_B(p_low, p_med, p_high, RR_med_crude, RR_hi_crude),
    PAF_B_adj   = calc_paf_B(p_low, p_med, p_high, RR_med_adj,   RR_hi_adj)
  ) %>%
  ungroup() %>%
  mutate(
    PAF_A_crude_pct = sprintf("%.1f%%", 100 * PAF_A_crude),
    PAF_A_adj_pct   = sprintf("%.1f%%", 100 * PAF_A_adj),
    PAF_B_crude_pct = sprintf("%.1f%%", 100 * PAF_B_crude),
    PAF_B_adj_pct   = sprintf("%.1f%%", 100 * PAF_B_adj)
  )

final_with_paf <- final_tbl %>%
  left_join(
    paf_tbl %>% select(outcome, 
                       PAF_A_crude, PAF_A_adj, PAF_B_crude, PAF_B_adj,
                       PAF_A_crude_pct, PAF_A_adj_pct, PAF_B_crude_pct, PAF_B_adj_pct),
    by = "outcome"
  )

write.csv(final_with_paf, "../result/5_2_2_cox_BEEM_with_PAF_AB.csv", row.names = FALSE)

outcomes_in_tbl <- unique(final_with_paf$outcome)
event_counts <- tibble(
  outcome = outcomes_in_tbl,
  deaths_total = sapply(outcomes_in_tbl, function(v) sum(morta[[v]] == 1, na.rm = TRUE))
)

final_with_paf_counts <- final_with_paf %>%
  left_join(event_counts, by = "outcome") %>%
  mutate(
    deaths_delta_A_crude = round(PAF_A_crude * deaths_total),
    deaths_delta_A_adj   = round(PAF_A_adj   * deaths_total),
    deaths_delta_B_crude = round(PAF_B_crude * deaths_total),
    deaths_delta_B_adj   = round(PAF_B_adj   * deaths_total)
  )

write.csv(final_with_paf_counts, "../result/5_2_3_cox_BEEM_with_PAF_AB_and_deaths.csv", row.names = FALSE)
