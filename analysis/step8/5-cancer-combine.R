if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step8/"))
getwd()
getwd()

library(tidyverse)
library(readxl)
library(survival)
library(furrr)
library(data.table)
library(readr)
library(future)

input_dir <- file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step8/results/")
output_dir <- file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step8/results/combine/")

n_boot <- as.numeric(Sys.getenv("SUN_BEEM_N_BOOT", unset = "1000"))
bootstrap_seed <- as.numeric(Sys.getenv("SUN_BEEM_SEED", unset = "123"))
if (!is.finite(n_boot) || n_boot < 10 || n_boot != floor(n_boot) || n_boot > .Machine$integer.max ||
    !is.finite(bootstrap_seed) || bootstrap_seed < 0 || bootstrap_seed != floor(bootstrap_seed) ||
    bootstrap_seed > .Machine$integer.max) stop("Invalid bootstrap execution settings.")
run_name <- paste0("bootstrap_", as.integer(n_boot), "_seed_", as.integer(bootstrap_seed))
mediation_dir <- file.path(input_dir, run_name)
output_dir <- file.path(output_dir, run_name, "")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- "cancer"  # 或 "allcause"、"diabetes" 等

med <- read_csv(file.path(mediation_dir, paste0("olink_", prefix, "_mediation_BEEM.csv")), show_col_types = FALSE)# 读取中介分析结果和线性回归结果
lm <- read_csv(paste0(input_dir, "2_olink_BEEM_lm_trend_results_FDR.csv"), show_col_types = FALSE)

merged <- med %>%
  inner_join(lm, by = "biomarker") %>%
  rename(
    mediation_effect = indirect_effect,
    mediation_ci_low = indirect_ci_low,
    mediation_ci_high = indirect_ci_high,
    prop_mediated = proportion_mediated,
    prop_mediated_ci_low = prop_mediated_ci_low,
    prop_mediated_ci_high = prop_mediated_ci_high,
    uv_beta = estimate,
    uv_fdr = p_adj
  ) %>%
  mutate(
    mediation_direction = case_when(
      mediation_effect > 0 ~ "positive",
      mediation_effect < 0 ~ "negative",
      TRUE ~ "neutral"
    ),
    uv_direction = case_when(
      uv_beta > 0 ~ "increase",
      uv_beta < 0 ~ "decrease",
      TRUE ~ "neutral"
    ),
    mediation_type = case_when(
      uv_beta > 0 & mediation_effect > 0 ~ "UV increases harmful mediator",
      uv_beta > 0 & mediation_effect < 0 ~ "UV increases protective mediator",
      uv_beta < 0 & mediation_effect > 0 ~ "UV reduces protective mediator",
      uv_beta < 0 & mediation_effect < 0 ~ "UV reduces harmful mediator",
      TRUE ~ "unclear"
    )
  )

n_significant_mediation <- merged %>%
  filter(mediation_ci_low > 0 | mediation_ci_high < 0) %>%
  nrow()

print(n_significant_mediation)

n_uv_fdr_significant <- merged %>%
  filter(uv_fdr < 0.01) %>%
  nrow()

print(n_uv_fdr_significant)

significant_hits <- merged %>%
  filter(
    (mediation_ci_low > 0 | mediation_ci_high < 0),  # CI 不跨 0
    uv_fdr < 0.01  # UV → biomarker 显著
  )
write_csv(significant_hits, paste0(output_dir, prefix, "_mt_uv_mediation_significant_hits_annotated.csv"))

nonsignificant_hits <- merged %>%
  filter(
    mediation_ci_low <= 0 & mediation_ci_high >= 0,
    uv_fdr >= 0.01
  )

write_csv(nonsignificant_hits, paste0(output_dir, prefix, "_uvmediation_CIcross0_UVns.csv"))
