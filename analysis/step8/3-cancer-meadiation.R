if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step8/"))
getwd()

library(tidyverse)
library(here)
library(readxl)
library(survival)
library(furrr)
library(data.table)
library(readr)
library(future)

mea_olink_set <- readRDS("mea_set_merged.rds")
biomarkers_olink <- readRDS("biomarkers_olink.rds")
options(future.globals.maxSize = 64 * 1024^3)  # 64GB
library(future)
library(furrr)
plan(multisession, workers = 4)

library(dplyr)
library(readr)
library(survival)
library(purrr)
library(data.table) # for fwrite

# Set SUN_BEEM_N_BOOT to the desired number of bootstrap resamples.
# The default preserves the setting in the archived source script.
n_boot <- as.numeric(Sys.getenv("SUN_BEEM_N_BOOT", unset = "50"))
if (length(n_boot) != 1L || !is.finite(n_boot) || n_boot < 2 ||
    n_boot != floor(n_boot) || n_boot > .Machine$integer.max) {
  stop("SUN_BEEM_N_BOOT must be an integer of at least 2.")
}
n_boot <- as.integer(n_boot)
message("Bootstrap resamples requested: ", n_boot)
batch_size <- 10

mea_olink_set <- mea_olink_set %>%
  mutate(
    BEEM = factor(BEEM, levels = c("low","medium","high")),
    BEEM_num = case_when(
      BEEM == "low" ~ 0,
      BEEM == "medium" ~ 1,
      BEEM == "high" ~ 2,
      TRUE ~ NA_real_
    )
  )

bootstrap_mediation <- function(data, biomarker) {
  base_formula <- as.formula(paste(
    "Surv(death_time, cancer_mortality) ~ BEEM_num + age + sex + bmi +",
    "deprivation_quintile + highest_edu + smoking_status + phys_activity + sleep_score + drinking"
  ))
  full_formula <- update(base_formula, paste("~ . +", biomarker))
  
  cox1 <- coxph(base_formula, data = data, ties = "efron")
  cox2 <- coxph(full_formula, data = data, ties = "efron")
  
  coef1 <- coef(cox1)[["BEEM_num"]]
  coef2 <- coef(cox2)[["BEEM_num"]]
  
  indirect <- coef1 - coef2
  prop_med <- ifelse(is.na(coef1) || coef1 == 0, NA_real_, indirect / coef1)
  
  c(indirect_effect = indirect, proportion_mediated = prop_med)
}

run_mediation <- function(biomarker) {
  cat("Running:", biomarker, "\n")
  
  biomarker_mea <- mea_olink_set %>%
    select(death_time, cancer_mortality, BEEM_num, age, sex, bmi,
           deprivation_quintile, highest_edu, smoking_status, phys_activity,
           sleep_score, drinking, !!sym(biomarker)) %>%
    tidyr::drop_na(BEEM_num, !!sym(biomarker), death_time, cancer_mortality)
  
  if (nrow(biomarker_mea) < 100 || dplyr::n_distinct(biomarker_mea$BEEM_num) < 2) return(NULL)
  
  base_formula <- as.formula(paste(
    "Surv(death_time, cancer_mortality) ~ BEEM_num + age + sex + bmi +",
    "deprivation_quintile + highest_edu + smoking_status + phys_activity + sleep_score + drinking"
  ))
  full_formula <- update(base_formula, paste("~ . +", biomarker))
  
  cox1 <- tryCatch(coxph(base_formula, data = biomarker_mea, ties = "efron"), error = function(e) NULL)
  cox2 <- tryCatch(coxph(full_formula, data = biomarker_mea, ties = "efron"), error = function(e) NULL)
  if (is.null(cox1) || is.null(cox2)) return(NULL)
  
  coef1 <- coef(cox1)[["BEEM_num"]]
  coef2 <- coef(cox2)[["BEEM_num"]]
  if (is.na(coef1) || is.na(coef2)) return(NULL)
  
  indirect_obs <- coef1 - coef2
  prop_med_obs <- ifelse(coef1 == 0, NA_real_, indirect_obs / coef1)
  
  boot_res <- replicate(n_boot, {
    boot_idx <- sample(nrow(biomarker_mea), replace = TRUE)
    boot_sample <- biomarker_mea[boot_idx, , drop = FALSE]
    tryCatch({
      bootstrap_mediation(boot_sample, biomarker)
    }, error = function(e) c(indirect_effect = NA_real_, proportion_mediated = NA_real_))
  })
  
  boot_res <- t(boot_res)
  colnames(boot_res) <- c("indirect_effect", "proportion_mediated")
  boot_res <- boot_res[stats::complete.cases(boot_res), , drop = FALSE]
  if (nrow(boot_res) < max(10, ceiling(0.2 * n_boot))) return(NULL)
  
  indirect_ci <- stats::quantile(boot_res[, "indirect_effect"], c(0.025, 0.975), na.rm = TRUE)
  prop_ci     <- stats::quantile(boot_res[, "proportion_mediated"], c(0.025, 0.975), na.rm = TRUE)
  
  data.frame(
    biomarker = biomarker,
    total_effect = coef1,                         # log-HR per one BEEM level increase
    direct_effect = coef2,                        # adjusted for mediator
    indirect_effect = indirect_obs,               # mediated log-HR
    indirect_ci_low = indirect_ci[1],
    indirect_ci_high = indirect_ci[2],
    proportion_mediated = prop_med_obs,           # fraction on log-HR scale
    prop_mediated_ci_low = prop_ci[1],
    prop_mediated_ci_high = prop_ci[2]
  )
}

result_path <- "results/olink_cancer_mediation_BEEM.csv"

done_vars <- if (file.exists(result_path)) {
  readr::read_csv(result_path, show_col_types = FALSE) %>% dplyr::pull(biomarker) %>% unique()
} else {
  character()
}

biomarkers_todo <- setdiff(biomarkers_olink, done_vars)
if (length(biomarkers_todo) == 0) {
  message("All biomarkers already processed in: ", result_path)
} else {
  biomarkers_batches <- split(biomarkers_todo, ceiling(seq_along(biomarkers_todo) / batch_size))
  
  for (i in seq_along(biomarkers_batches)) {
    cat("\n▶▶ Batch", i, "of", length(biomarkers_batches), "\n")
    batch <- biomarkers_batches[[i]]
    
    batch_results <- furrr::future_map(
      batch, run_mediation, .progress = TRUE,
      .options = furrr::furrr_options(seed = TRUE)
    )
    batch_df <- dplyr::bind_rows(batch_results)
    
    if (nrow(batch_df) > 0) {
      data.table::fwrite(batch_df, result_path, append = file.exists(result_path))
    }
  }
}

plan(sequential)
