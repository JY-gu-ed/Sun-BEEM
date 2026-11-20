if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(here)
library(survival)
library(broom)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step5/"))

morta <- read.csv(here("washed_data", "4_2_2_ana_for_mortality.csv"))

morta <- morta %>%
  mutate(
    BEEM = factor(BEEM, levels = c("low","medium","high")),
    BEEM_numeric = as.numeric(BEEM),
    sex          = factor(sex, levels = c("female","male")),
    age          = factor(age, levels = c("37-47","48-54","55-59","60-63","64-73")),
    age_group2   = ifelse(age %in% c("37-47","48-54","55-59"), "<60", "≥60"),
    bmi            = factor(bmi, levels = c("underweight","healthy weight","overweight","obese")),
    sleep_score    = factor(sleep_score, levels = c("poor","intermediate","healthy")),
    drinking       = factor(drinking, levels = c("unhealthy","healthy")),
    phys_activity  = factor(phys_activity, levels = c("low","moderate","high")),
    highest_edu    = factor(highest_edu, levels = c("less than secondary","secondary school","university or above")),
    smoking_status = factor(smoking_status, levels = c("current","previous","never")),
    deprivation_quintile = factor(deprivation_quintile, levels = 1:5, ordered = TRUE),
    across(c(all_cause_mortality, cvd_mortality_narrow, cancer_mortality,
             skin_cancer_mortality, melanoma_skin_cancer_mortality,
             non_melanoma_mortality, non_skin_cancer_mortality),
           ~ ifelse(. == "yes", 1, 0))
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

run_one_group <- function(outcome, data, group_var, group_level, covariates){
  data_sub <- data %>% filter(!!sym(group_var) == group_level)
  covs <- covariates[!covariates %in% group_var]
  
  rhs <- paste(c("BEEM", covs), collapse = " + ")
  f_cat <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ ", rhs))
  fit_cat <- safe_fit(f_cat, data_sub)
  
  hr_df <- if (!is.null(fit_cat)) {
    tidy(fit_cat, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term %in% c("BEEMmedium","BEEMhigh")) %>%
      mutate(HR_CI = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high)) %>%
      select(term, estimate, conf.low, conf.high, HR_CI, p.value)
  } else {
    tibble(term=c("BEEMmedium","BEEMhigh"), estimate=NA, conf.low=NA, conf.high=NA,
           HR_CI=NA, p.value=NA)
  }
  
  f_trend <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ BEEM_numeric + ", paste(covs, collapse=" + ")))
  fit_trend <- safe_fit(f_trend, data_sub)
  p_trend <- if (!is.null(fit_trend)) {
    tidy(fit_trend)$p.value[ tidy(fit_trend)$term == "BEEM_numeric" ]
  } else NA
  
  hr_df %>%
    mutate(
      outcome = outcome,
      analysis_type = paste0(group_var,"_",group_level),
      p_trend = p_trend
    )
}

sex_levels <- levels(morta$sex)
sex_results <- map_dfr(sex_levels, ~map_dfr(mortality_vars, run_one_group,
                                            data=morta, group_var="sex", group_level=.x, covariates=covariates))

age_levels2 <- unique(morta$age_group2)
age_results <- map_dfr(age_levels2, ~map_dfr(mortality_vars, run_one_group,
                                             data=morta, group_var="age_group2", group_level=.x, covariates=covariates))

all_res <- bind_rows(sex_results, age_results)

outdir <- "../result/strata/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(all_res, paste0(outdir, "cox_BEEM_mortality_groupAnalysis_age60_withPtrend.csv"), row.names = FALSE)

print(all_res)
all_res <- all_res %>% dplyr::filter(outcome != "skin_cancer_mortality")  # ← 新增这一行

library(ggplot2)

outcome_labels_en <- c(
  non_melanoma_mortality         = "other skin cancer mortality",
  melanoma_skin_cancer_mortality = "Melanoma mortality",
  non_skin_cancer_mortality      = "non-skin cancer  mortality" , 
  cancer_mortality               = "Cancer mortality",
  cvd_mortality_narrow           = "CVD mortality",
  all_cause_mortality            = "All-cause mortality"
)

plot_df <- all_res %>%
  mutate(
    BEEM_group = recode(term,
                        "BEEMmedium" = "Medium vs Low",
                        "BEEMhigh"   = "High vs Low"),
    BEEM_group = factor(BEEM_group, levels = c("Medium vs Low","High vs Low")),
    outcome_label = factor(outcome,
                           levels = names(outcome_labels_en),
                           labels = unname(outcome_labels_en)),
    group = analysis_type
  )

plot_sex <- plot_df %>% filter(grepl("^sex_", group))
plot_age <- plot_df %>% filter(grepl("^age_", group))

make_forest_noPtrend <- function(df, file_png, file_pdf=NULL, group_title="Group"){
  pd <- position_dodge(width = 0.6)
  
  p <- ggplot(df, aes(x = estimate, y = outcome_label,
                      xmin = conf.low, xmax = conf.high,
                      color = group)) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_errorbarh(height = 0.2, position = pd) +
    geom_point(size = 2.8, position = pd) +
    geom_text(aes(x = conf.high, label = HR_CI),
              hjust = -0.1, vjust = 0.5, position = pd,
              size = 3.0, show.legend = FALSE) +
    facet_grid(. ~ BEEM_group) +
    scale_x_log10() +
    labs(x = "Hazard Ratio (log scale)", y = NULL, color = group_title) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 11))
  
  print(p)
  ggsave(file_png, p, width = 30, height = 6, dpi = 300, bg = "white")
  if(!is.null(file_pdf)) ggsave(file_pdf, p, width = 30, height = 6, bg = "white")
}

make_forest_noPtrend(plot_sex,
            paste0(outdir, "forest_BEEM_mortality_bySex_withPtrend.png"),
            paste0(outdir, "forest_BEEM_mortality_bySex_withPtrend.pdf"),
            group_title = "Sex")

make_forest_noPtrend(plot_age,
            paste0(outdir, "forest_BEEM_mortality_byAge60_withPtrend.png"),
            paste0(outdir, "forest_BEEM_mortality_byAge60_withPtrend.pdf"),
            group_title = "Age group")

plot_sex_shape <- plot_sex %>%
  mutate(
    group_color = recode(group,
                         "sex_female" = "Female",
                         "sex_male"   = "Male")
  )

plot_age_shape <- plot_age %>%
  mutate(
    group_color = recode(group,
                         "age_group2_<60" = "<60",
                         "age_group2_≥60" = "≥60")
  )

make_forest_shape_color <- function(df, file_png, file_pdf = NULL, group_title = "Group") {
  pd <- position_dodge(width = 0.6)
  
  p <- ggplot(
    df,
    aes(
      x    = estimate,
      y    = outcome_label,
      xmin = conf.low,
      xmax = conf.high,
      color = group_color,      # 颜色表示性别/年龄
      shape = BEEM_group        # 形状表示 Medium vs High
    )
  ) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_errorbarh(height = 0.2, position = pd) +
    geom_point(size = 2.8, position = pd) +
    geom_text(
      aes(
        x = conf.high,
        label = HR_CI
      ),
      hjust = -0.1,
      vjust = 0.5,
      position = pd,
      size = 3.0,
      show.legend = FALSE
    ) +
    scale_x_log10() +
    labs(
      x = "Hazard Ratio (log scale)",
      y = NULL,
      color = group_title,
      shape = "UV exposure group"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 11)
    )
  
  print(p)
  ggsave(file_png, p, width = 20, height = 10, dpi = 300, bg = "white")
  if (!is.null(file_pdf)) {
    ggsave(file_pdf, p, width = 20, height = 10, bg = "white")
  }
}

make_forest_shape_color(
  plot_sex_shape,
  paste0(outdir, "forest_BEEM_mortality_bySex_shapeColor.png"),
  paste0(outdir, "forest_BEEM_mortality_bySex_shapeColor.pdf"),
  group_title = "Sex"
)

make_forest_shape_color(
  plot_age_shape,
  paste0(outdir, "forest_BEEM_mortality_byAge60_shapeColor.png"),
  paste0(outdir, "forest_BEEM_mortality_byAge60_shapeColor.pdf"),
  group_title = "Age group"
)
