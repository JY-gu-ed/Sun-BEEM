if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(survival)
library(broom)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step5/"))

morta <- read.csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "washed_data", "4_2_2_ana_for_mortality.csv"))

morta_backup <- morta

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
                                  ordered = TRUE),  # 有序因子
    across(c(all_cause_mortality, cvd_mortality_narrow, cancer_mortality,
             skin_cancer_mortality, melanoma_skin_cancer_mortality,
             non_melanoma_mortality, non_skin_cancer_mortality),
           ~ ifelse(. == "yes", 1, 0)),
    BEEM_numeric = as.numeric(BEEM)
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

run_one <- function(outcome, data) {
  f_crude_cat <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ BEEM"))
  fit_crude_cat <- safe_fit(f_crude_cat, data)
  
  crude_df <- if (!is.null(fit_crude_cat)) {
    tidy(fit_crude_cat, exponentiate = TRUE, conf.int = TRUE) |>
      filter(term %in% c("BEEMmedium","BEEMhigh")) |>
      mutate(HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
             p = sprintf("%.4f", p.value)) |>
      select(term, HR_CI, p)
  } else tibble(term=c("BEEMmedium","BEEMhigh"), HR_CI=NA, p=NA)
  
  rhs <- paste(c("BEEM", covariates), collapse = " + ")
  f_adj_cat <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ ", rhs))
  fit_adj_cat <- safe_fit(f_adj_cat, data)
  
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
  fit_crude_trend <- safe_fit(f_crude_trend, data)
  fit_adj_trend   <- safe_fit(f_adj_trend,   data)
  
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

datasets <- list(
  main  = morta,
  sens2 = morta %>% filter(death_time > 730),   # 排除前2年死亡
  sens5 = morta %>% filter(death_time > 1825)   # 排除前5年死亡
)

outdir <- "../result/sen_ana/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

results <- lapply(names(datasets), function(name){
  res <- map_dfr(mortality_vars, run_one, data = datasets[[name]]) %>%
    mutate(analysis_type = name)
  out_file <- paste0(outdir, "cox_BEEM_", name, ".csv")
  write.csv(res, out_file, row.names = FALSE)
  message("已保存: ", out_file)
  res
})

final_tbl_all <- bind_rows(results)
write.csv(final_tbl_all, paste0(outdir, "cox_BEEM_all_combined.csv"), row.names = FALSE)
message("合并结果表已保存: ", paste0(outdir, "cox_BEEM_all_combined.csv"))

final_tbl_all

library(ggplot2)
library(stringr)
library(tidyr)
library(purrr)

outcome_labels_en <- c(
  all_cause_mortality            = "All-cause mortality",
  cvd_mortality_narrow           = "CVD mortality",
  cancer_mortality               = "Cancer mortality",
  melanoma_skin_cancer_mortality = "Melanoma mortality",
  non_melanoma_mortality         = "Other skin cancer mortality",
  non_skin_cancer_mortality      = "Non-skin cancer mortality"
)

parse_hrci <- function(x){
  ifelse(is.na(x), NA, x) %>%
    str_match("([0-9.]+) \\(([0-9.]+)-([0-9.]+)\\)") %>%
    as.data.frame() %>%
    transmute(
      estimate = as.numeric(V2),
      conf.low = as.numeric(V3),
      conf.high= as.numeric(V4)
    )
}

hr_long <- final_tbl_all %>%
  pivot_longer(
    cols = c(crude_HR_medium, crude_HR_high, adj_HR_medium, adj_HR_high),
    names_to = c("model","metric","group"),
    names_sep = "_",
    values_to = "HR_CI"
  ) %>%
  filter(metric == "HR") %>%
  mutate(tmp = map(HR_CI, parse_hrci)) %>%
  unnest(tmp) %>%
  mutate(
    model = recode(model, crude = "Crude", adj = "Adjusted"),
    BEEM_group = ifelse(group == "medium", "Medium vs Low", "High vs Low"),
    BEEM_group = factor(BEEM_group, levels = c("Medium vs Low","High vs Low")),
    outcome_label = factor(outcome,
                           levels = names(outcome_labels_en),
                           labels = unname(outcome_labels_en)),
    analysis_type = factor(analysis_type, levels = c("main","sens2","sens5"))
  )

make_forest <- function(df, outcome_order, file_png, file_pdf = NULL){
  df2 <- df %>%
    filter(outcome %in% outcome_order) %>%
    mutate(
      outcome_label = factor(outcome,
                             levels = outcome_order,
                             labels = unname(outcome_labels_en[outcome_order]))
    )
  
  pd <- position_dodge(width = 0.6)
  
  p <- ggplot(df2, aes(x = estimate, y = outcome_label,
                       xmin = conf.low, xmax = conf.high,
                       color = analysis_type, shape = model)) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_errorbarh(height = 0.2, position = pd) +
    geom_point(size = 2.8, position = pd) +
    geom_text(aes(x = conf.high,
                  label = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high)),
              hjust = -0.15, vjust = 0.5, position = pd, size = 3.2, show.legend = FALSE) +
    facet_grid(. ~ BEEM_group) +
    scale_x_log10() +
    labs(x = "Hazard Ratio (log scale)",
         y = NULL,
         color = "Analysis",
         shape = "Model") +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 10)
    )
  
  print(p)
  ggsave(file_png, p, width = 25, height = 6, dpi = 300, bg = "white")
  if(!is.null(file_pdf)) ggsave(file_pdf, p, width = 25, height = 6, bg = "white")
}

set1 <- c("non_skin_cancer_mortality",
          "cancer_mortality",
          "cvd_mortality_narrow",
          "all_cause_mortality")

set2 <- c("non_melanoma_mortality",
          "melanoma_skin_cancer_mortality",
          "skin_cancer_mortality")

make_forest(hr_long, set1,
            "../result/sen_ana/cox_BEEM_all_set1.png",
            "../result/sen_ana/cox_BEEM_all_set1.pdf")

make_forest(hr_long, set2,
            "../result/sen_ana/cox_BEEM_all_set2.png",
            "../result/sen_ana/cox_BEEM_all_set2.pdf")

all_outcomes <- c(

  "non_melanoma_mortality",
  "melanoma_skin_cancer_mortality",  
  "non_skin_cancer_mortality",
  "cancer_mortality",
"cvd_mortality_narrow",
"all_cause_mortality")

make_forest_sens <- function(df, ana, outcome_order, file_png, file_pdf = NULL) {
  df2 <- df %>%
    dplyr::filter(
      analysis_type == ana,          # sens2 或 sens5
      model == "Adjusted",           # 只看矫正模型
      outcome %in% outcome_order
    ) %>%
    mutate(
      outcome_label = factor(
        outcome,
        levels = outcome_order,
        labels = unname(outcome_labels_en[outcome_order])
      )
    )
  
  pd <- position_dodge(width = 0.6)
  
  p <- ggplot(
    df2,
    aes(
      x    = estimate,
      y    = outcome_label,
      xmin = conf.low,
      xmax = conf.high,
      shape = BEEM_group,            # 用形状区分 Medium vs High
      group = BEEM_group
    )
  ) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_errorbarh(height = 0.2, position = pd) +
    geom_point(size = 2.8, position = pd) +
    geom_text(
      aes(
        x = conf.high,
        label = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high)
      ),
      hjust = -0.15,
      vjust = 0.5,
      position = pd,
      size = 3.2,
      show.legend = FALSE
    ) +
    scale_x_log10() +
    scale_shape_discrete(name = "UV exposure group") +
    labs(
      x = "Hazard Ratio (log scale)",
      y = NULL,
      title = ifelse(
        ana == "sens2",
        "Adjusted HRs for mortality outcomes\n(excluding deaths in first 2 years)",
        "Adjusted HRs for mortality outcomes\n(excluding deaths in first 5 years)"
      )
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      axis.text.y      = element_text(size = 10),
      plot.title       = element_text(hjust = 0.5)
    )
  
  print(p)
  ggsave(file_png, p, width = 28, height = 8, dpi = 300, bg = "white")
  if (!is.null(file_pdf)) {
    ggsave(file_pdf, p, width = 28, height = 8, bg = "white")
  }
}

make_forest_sens(
  hr_long,
  ana = "sens2",
  outcome_order = all_outcomes,
  file_png = "../result/sen_ana/cox_BEEM_adj_sens2_all.png",
  file_pdf = "../result/sen_ana/cox_BEEM_adj_sens2_all.pdf"
)

make_forest_sens(
  hr_long,
  ana = "sens5",
  outcome_order = all_outcomes,
  file_png = "../result/sen_ana/cox_BEEM_adj_sens5_all.png",
  file_pdf = "../result/sen_ana/cox_BEEM_adj_sens5_all.pdf"
)
