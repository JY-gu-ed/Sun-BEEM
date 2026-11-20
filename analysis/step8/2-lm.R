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
setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step8/"))

getwd()

mea_olink_set <- readRDS("mea_set_merged.rds")
message("Protein analysis rows: ", nrow(mea_olink_set))
if (!dir.exists("results")) dir.create("results", recursive = TRUE)
biomarkers_olink <- readRDS("biomarkers_olink.rds")

olink_lm_results <- list()

mea_olink_set <- mea_olink_set %>%
  mutate(
    BEEM_num = case_when(
      BEEM == "low" ~ 1,
      BEEM == "medium" ~ 2,
      BEEM == "high" ~ 3,
      TRUE ~ NA_real_
    )
  )

covariate_formula <- paste(
  "BEEM_num",
  "age + sex + bmi + deprivation_quintile + highest_edu + smoking_status + phys_activity + sleep_score + drinking",
  sep = " + "
)

for (biomarker in biomarkers_olink) {
  df <- mea_olink_set %>%
    select(all_of(biomarker), BEEM_num, age, sex, bmi,
           deprivation_quintile, highest_edu, smoking_status,
           phys_activity, sleep_score, drinking) %>%
    drop_na()
  
  if (n_distinct(df[[biomarker]]) < 2 || n_distinct(df$BEEM_num) < 2) next
  
  formula <- as.formula(paste(biomarker, "~", covariate_formula))
  
  fit <- tryCatch(lm(formula, data = df), error = function(e) NULL)
  if (is.null(fit)) next
  
  tt <- tidy(fit, conf.int = TRUE) %>%
    filter(term == "BEEM_num") %>%
    transmute(
      biomarker = biomarker,
      estimate  = estimate,
      conf.low  = conf.low,
      conf.high = conf.high,
      p_trend   = p.value
    )
  
  olink_lm_results[[biomarker]] <- tt
}

olink_lm_trend_df <- bind_rows(olink_lm_results) %>%
  arrange(p_trend)

olink_lm_trend_df <- olink_lm_trend_df %>%
  mutate(
    p_adj = p.adjust(p_trend, method = "fdr"),
    direction = case_when(
      p_adj >= 0.05 ~ "no_change",
      p_adj < 0.05 & estimate > 0 ~ "up",
      p_adj < 0.05 & estimate < 0 ~ "down"
    )
  )

table(olink_lm_trend_df$direction)

write_csv(olink_lm_trend_df, "results/2_olink_BEEM_lm_trend_results_FDR.csv")

library(dplyr)
library(broom)
library(readr)
library(tidyr)

mea_olink_set <- readRDS("mea_set_merged.rds")
biomarkers_olink <- readRDS("biomarkers_olink.rds")

exposure_vars <- c("bi_outdoor", "bi_solarium_use", "bi_uv_protection", "bi_swr2008")

for (exposure in exposure_vars) {
  
  olink_lm_results <- list()
  
  mea_olink_set <- mea_olink_set %>%
    mutate(!!exposure := factor(.data[[exposure]], levels = c(0, 1)))
  
  covariate_formula <- paste(
    exposure,
    "age + sex + bmi + deprivation_quintile + highest_edu + smoking_status + phys_activity + sleep_score + drinking",
    sep = " + "
  )
  
  for (biomarker in biomarkers_olink) {
    df <- mea_olink_set %>%
      select(all_of(biomarker), all_of(exposure), age, sex, bmi,
             deprivation_quintile, highest_edu, smoking_status,
             phys_activity, sleep_score, drinking) %>%
      drop_na()
    
    if (n_distinct(df[[biomarker]]) < 2 || n_distinct(df[[exposure]]) < 2) next
    
    formula <- as.formula(paste(biomarker, "~", covariate_formula))
    fit <- tryCatch(lm(formula, data = df), error = function(e) NULL)
    if (is.null(fit)) next
    
    tt <- tidy(fit, conf.int = TRUE)
    term_name <- if (paste0(exposure, "1") %in% tt$term) paste0(exposure, "1") else exposure
    
    tt <- tt %>%
      filter(term == term_name) %>%
      transmute(
        biomarker = biomarker,
        estimate  = estimate,   # 均值(1) − 均值(0)
        conf.low  = conf.low,
        conf.high = conf.high,
        p_value   = p.value
      )
    
    olink_lm_results[[biomarker]] <- tt
  }
  
  olink_lm_df <- bind_rows(olink_lm_results) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      direction = case_when(
        p_value >= 0.05 ~ "no_change",
        p_value < 0.05 & estimate > 0 ~ "up (1 > 0)",
        p_value < 0.05 & estimate < 0 ~ "down (1 < 0)"
      )
    ) %>%
    arrange(p_adj)
  
  cat("\n", exposure, ":\n")
  print(table(olink_lm_df$direction))
  
  write_csv(olink_lm_df, paste0("results/olink_", exposure, "_results_FDR.csv"))
}

library(dplyr)
library(broom)
library(readr)
library(purrr)
library(tidyr)

mea_olink_set <- readRDS("mea_set_merged.rds")
biomarkers_olink <- readRDS("biomarkers_olink.rds")  # 字符向量：biomarker 列名

covars <- c("age","sex","bmi","deprivation_quintile","highest_edu",
            "smoking_status","phys_activity","sleep_score","drinking")
covar_str <- paste(covars, collapse = " + ")

mea_olink_set <- mea_olink_set %>%
  mutate(uva_2dayavg = as.numeric(uva_2dayavg))

results <- map_dfr(biomarkers_olink, function(bm) {
  df <- mea_olink_set %>%
    select(all_of(bm), uva_2dayavg, any_of(covars)) %>%
    drop_na()
  
  if (n_distinct(df[[bm]]) < 2 || n_distinct(df$uva_2dayavg) < 2) return(NULL)
  if (!is.numeric(df[[bm]])) df[[bm]] <- as.numeric(df[[bm]])
  
  form <- as.formula(paste(bm, "~ uva_2dayavg +", covar_str))
  fit  <- tryCatch(lm(form, data = df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  tidy(fit, conf.int = TRUE) %>%
    filter(term == "uva_2dayavg") %>%
    transmute(
      biomarker = bm,
      beta      = estimate,      # 每 ↑1 单位 UVA 对 biomarker 的效应
      std.error = std.error,
      conf.low  = conf.low,
      conf.high = conf.high,
      p_value   = p.value,
      n         = nrow(df)
    )
})

olink_uva_df <- results %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    direction = case_when(
      p_value >= 0.05 ~ "no_change",
      beta  >  0   ~ "positive",
      TRUE         ~ "negative"
    )
  ) %>%
  arrange(p_value)

print(table(olink_uva_df$direction))
write_csv(olink_uva_df, "results/olink_biomarker_vs_uva2dayavg_FDR.csv")

library(dplyr)
library(broom)
library(readr)
library(purrr)
library(tidyr)

mea_olink_set <- readRDS("mea_set_merged.rds")
biomarkers_olink <- readRDS("biomarkers_olink.rds")  # biomarker 名称向量

covars <- c("age","sex","bmi","deprivation_quintile","highest_edu",
            "smoking_status","phys_activity","sleep_score","drinking")
covar_str <- paste(covars, collapse = " + ")

uva_vars <- c("uva_7_dayavg", "uva_14_dayavg", "uva_30_dayavg")

mea_olink_set <- mea_olink_set %>%
  mutate(across(all_of(uva_vars), as.numeric))

for (uva in uva_vars) {
  results <- map_dfr(biomarkers_olink, function(bm) {
    df <- mea_olink_set %>%
      select(all_of(bm), all_of(uva), any_of(covars)) %>%
      drop_na()
    
    if (n_distinct(df[[bm]]) < 2 || n_distinct(df[[uva]]) < 2) return(NULL)
    if (!is.numeric(df[[bm]])) df[[bm]] <- as.numeric(df[[bm]])
    
    form <- as.formula(paste(bm, "~", uva, "+", covar_str))
    fit  <- tryCatch(lm(form, data = df), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    
    tidy(fit, conf.int = TRUE) %>%
      filter(term == uva) %>%
      transmute(
        biomarker = bm,
        beta      = estimate,
        std.error = std.error,
        conf.low  = conf.low,
        conf.high = conf.high,
        p_value   = p.value,
        n         = nrow(df)
      )
  })
  
  results <- results %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      direction = case_when(
        p_value >= 0.05 ~ "no_change",
        beta  >  0      ~ "positive",
        TRUE            ~ "negative"
      )
    ) %>%
    arrange(p_value)
  
  fout <- paste0("results/olink_biomarker_vs_", uva, "_FDR.csv")
  write_csv(results, fout)
  
  message("✅ 保存结果: ", fout)
}
