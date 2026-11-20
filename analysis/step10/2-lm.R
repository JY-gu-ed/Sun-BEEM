if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(broom)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step10/"))
getwd()

mea_nmr_set   <- readRDS("mea_set_merged.rds")
biomarkers_nmr <- readRDS("biomarkers_nmr.rds")

message("N individuals: ", nrow(mea_nmr_set))
message("N metabolites: ", length(biomarkers_nmr))

mea_nmr_set <- mea_nmr_set %>%
  mutate(
    UV_num = case_when(
      BEEM == "low" ~ 1,
      BEEM == "medium" ~ 2,
      BEEM == "high" ~ 3,
      TRUE ~ NA_real_
    )
  )

covariate_formula <- paste(
  "UV_num",
  "age + sex + bmi + deprivation_quintile + highest_edu + smoking_status + phys_activity + sleep_score + drinking",
  sep = " + "
)

nmr_lm_results <- list()

for (biomarker in biomarkers_nmr) {
  
  df <- mea_nmr_set %>%
    dplyr::select(all_of(biomarker), UV_num, age, sex, bmi,
           deprivation_quintile, highest_edu, smoking_status,
           phys_activity, sleep_score, drinking) %>%
    drop_na()
  
  if (n_distinct(df[[biomarker]]) < 2 || n_distinct(df$UV_num) < 2) next
  
  if (!is.numeric(df[[biomarker]])) df[[biomarker]] <- as.numeric(df[[biomarker]])
  
  formula <- as.formula(paste(biomarker, "~", covariate_formula))
  
  fit <- tryCatch(lm(formula, data = df), error = function(e) NULL)
  if (is.null(fit)) next
  
  tt <- tidy(fit, conf.int = TRUE) %>%
    filter(term == "UV_num") %>%
    transmute(
      biomarker = biomarker,
      estimate  = estimate,
      conf.low  = conf.low,
      conf.high = conf.high,
      p_trend   = p.value,
      n         = nrow(df)
    )
  
  nmr_lm_results[[biomarker]] <- tt
}

nmr_lm_trend_df <- bind_rows(nmr_lm_results) %>%
  arrange(p_trend) %>%
  mutate(
    p_adj = p.adjust(p_trend, method = "fdr"),
    direction = case_when(
      p_adj >= 0.05 ~ "no_change",
      p_adj < 0.05 & estimate > 0 ~ "up",
      p_adj < 0.05 & estimate < 0 ~ "down"
    )
  )

print(table(nmr_lm_trend_df$direction, useNA = "ifany"))

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

write_csv(nmr_lm_trend_df, "results/2_nmr_UV_lm_trend_results_FDR.csv")
message("✅ Saved: results/2_nmr_UV_lm_trend_results_FDR.csv")
