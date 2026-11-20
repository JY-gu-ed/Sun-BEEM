if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(survival)
library(broom)
library(purrr)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step5/"))

morta <- read.csv("../washed_data/4_2_2_ana_for_mortality.csv")

morta <- morta %>%
  dplyr::mutate(
    BEEM = factor(BEEM, levels = c("low", "medium", "high")),
    BEEM_numeric = as.numeric(BEEM),  # low=1, medium=2, high=3
    sex            = factor(sex, levels = c("female","male")),
    age            = factor(age, levels = c("37-47","48-54","55-59","60-63","64-73")),
    bmi            = factor(bmi, levels = c("underweight","healthy weight","overweight","obese")),
    sleep_score    = factor(sleep_score, levels = c("poor","intermediate","healthy")),
    drinking       = factor(drinking, levels = c("unhealthy","healthy")),
    phys_activity  = factor(phys_activity, levels = c("low","moderate","high")),
    highest_edu    = factor(highest_edu, levels = c("less than secondary","secondary school","university or above")),
    smoking_status = factor(smoking_status, levels = c("current","previous","never")),
    deprivation_quintile = factor(deprivation_quintile, levels = 1:5, ordered = TRUE)
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

covariates <- c("sex","age","bmi","deprivation_quintile",
                "highest_edu","phys_activity","sleep_score",
                "smoking_status","drinking")

safe_fit <- function(formula, data) {
  tryCatch(coxph(formula, data = data, na.action = na.exclude), error = function(e) NULL)
}

run_one <- function(outcome) {
  results <- list()
  
  f_crude <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ BEEM"))
  fit_crude <- safe_fit(f_crude, morta)
  
  crude_df <- if (!is.null(fit_crude)) {
    tidy(fit_crude, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term %in% c("BEEMmedium","BEEMhigh")) %>%
      mutate(model = "Model 1: crude") %>%
      select(model, term, estimate, conf.low, conf.high, p.value)
  } else tibble(model="Model 1: crude", term=c("BEEMmedium","BEEMhigh"),
                estimate=NA, conf.low=NA, conf.high=NA, p.value=NA)
  
  f_crude_trend <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ BEEM_numeric"))
  fit_crude_trend <- safe_fit(f_crude_trend, morta)
  p_trend_crude <- if (!is.null(fit_crude_trend))
    tidy(fit_crude_trend)$p.value[ tidy(fit_crude_trend)$term == "BEEM_numeric" ] else NA
  
  crude_df <- crude_df %>%
    mutate(outcome = outcome,
           HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value),
           p_trend = sprintf("%.4f", p_trend_crude))
  
  results[[1]] <- crude_df
  
  for (i in seq_along(covariates)) {
    rhs <- paste(c("BEEM", covariates[1:i]), collapse = " + ")
    f <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ ", rhs))
    fit <- safe_fit(f, morta)
    
    adj_df <- if (!is.null(fit)) {
      tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
        filter(term %in% c("BEEMmedium","BEEMhigh")) %>%
        mutate(model = paste0("Model ", i+1, ": crude + ", 
                              paste(covariates[1:i], collapse = " + "))) %>%
        select(model, term, estimate, conf.low, conf.high, p.value)
    } else tibble(model=paste0("Model ", i+1, ": crude + ",
                               paste(covariates[1:i], collapse = " + ")),
                  term=c("BEEMmedium","BEEMhigh"),
                  estimate=NA, conf.low=NA, conf.high=NA, p.value=NA)
    
    f_trend <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ BEEM_numeric + ",
                                 paste(covariates[1:i], collapse = " + ")))
    fit_trend <- safe_fit(f_trend, morta)
    p_trend <- if (!is.null(fit_trend))
      tidy(fit_trend)$p.value[ tidy(fit_trend)$term == "BEEM_numeric" ] else NA
    
    adj_df <- adj_df %>%
      mutate(outcome = outcome,
             HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
             p = sprintf("%.4f", p.value),
             p_trend = sprintf("%.4f", p_trend))
    
    results[[i+1]] <- adj_df
  }
  
  bind_rows(results)
}

final_tbl <- map_dfr(mortality_vars, run_one)

write.csv(final_tbl, "cox_BEEM_stepwise_models_with_ptrend.csv", row.names=FALSE)

head(final_tbl)

library(ggplot2)
library(stringr)
library(grid)

parse_hrci <- function(estimate, low, high){
  sprintf("%.2f (%.2f–%.2f)", estimate, low, high)
}

plot_df <- final_tbl %>%
  mutate(
    model = factor(model, levels = unique(model)),  # 按顺序排列模型
    HR_text = parse_hrci(estimate, conf.low, conf.high),
    outcome_label = outcome
  )

outcome_labels_en <- c(
  all_cause_mortality            = "All-cause mortality",
  cvd_mortality_narrow           = "CVD mortality (narrow)",
  cancer_mortality               = "Cancer mortality",
  skin_cancer_mortality          = "Skin cancer mortality",
  melanoma_skin_cancer_mortality = "Melanoma mortality",
  non_melanoma_mortality         = "Non-melanoma skin cancer mortality",
  non_skin_cancer_mortality      = "Non-skin cancer mortality"
)

plot_df$outcome_label <- factor(plot_df$outcome,
                                levels = mortality_vars,
                                labels = unname(outcome_labels_en[mortality_vars]))

make_forest <- function(df, outcome_order, file_png, file_pdf=NULL){
  df2 <- df %>%
    filter(outcome %in% outcome_order) %>%
    mutate(
      outcome_label = factor(outcome,
                             levels = outcome_order,
                             labels = unname(outcome_labels_en[outcome_order]))
    )
  
  pd <- position_dodge(width = 0.8)
  
  x_min <- min(df2$conf.low, na.rm=TRUE) / 1.5
  x_max <- max(df2$conf.high, na.rm=TRUE) * 1.5
  
  p <- ggplot(df2, aes(y=model, x=estimate, xmin=conf.low, xmax=conf.high,
                       color=term)) +
    geom_vline(xintercept = 1, linetype=2) +
    geom_errorbarh(height=0.2, position=pd) +
    geom_point(size=2.5, position=pd) +
    geom_text(aes(x=conf.high, label=HR_text),
              hjust=-0.1, vjust=0.5, size=3, position=pd, show.legend=FALSE) +
    scale_x_log10(limits=c(x_min, x_max)) +
    scale_color_manual(values=c("BEEMmedium"="#1f78b4","BEEMhigh"="#e31a1c"),
                       labels=c("BEEMmedium"="Medium vs Low","BEEMhigh"="High vs Low")) +
    labs(x="Hazard Ratio (log scale)", y=NULL, color=NULL) +
    facet_grid(outcome_label ~ ., scales="free_y", space="free_y") +
    theme_minimal() +
    theme(
      legend.position="bottom",
      axis.text.y=element_text(size=9),
      strip.text.y=element_text(size=11, face="bold"),
      plot.margin=margin(5.5, 12, 5.5, 5.5)
    )
  
  print(p)
  ggsave(file_png, p, width=12, height=12, dpi=300, bg="white")
  if(!is.null(file_pdf)) ggsave(file_pdf, p, width=12, height=12, bg="white")
}
make_forest_shape <- function(df, outcome_order, file_png, file_pdf = NULL){
  df2 <- df %>%
    dplyr::filter(outcome %in% outcome_order) %>%
    dplyr::mutate(
      outcome_label = factor(
        outcome,
        levels = outcome_order,
        labels = unname(outcome_labels_en[outcome_order])
      )
    )
  
  pd <- position_dodge(width = 0.8)
  
  x_min <- min(df2$conf.low, na.rm = TRUE) / 1.5
  x_max <- max(df2$conf.high, na.rm = TRUE) * 1.5
  
  p <- ggplot(
    df2,
    aes(
      y     = model,
      x     = estimate,
      xmin  = conf.low,
      xmax  = conf.high,
      shape = term          # ★ 用形状区分 BEEM 组别
    )
  ) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_errorbarh(height = 0.2, position = pd) +
    geom_point(size = 2.5, position = pd) +
    geom_text(
      aes(x = conf.high, label = HR_text),
      hjust = -0.1, vjust = 0.5, size = 3,
      position = pd, show.legend = FALSE
    ) +
    scale_x_log10(limits = c(x_min, x_max)) +
    scale_shape_manual(
      values = c(
        BEEMmedium = 16,  # 圆点 = medium vs low
        BEEMhigh   = 15   # 方块 = high vs low
      ),
      labels = c(
        BEEMmedium = "Medium vs Low",
        BEEMhigh   = "High vs Low"
      )
    ) +
    labs(
      x     = "Hazard Ratio (log scale)",
      y     = NULL,
      shape = NULL
    ) +
    facet_wrap(
      ~ outcome_label,
      ncol = 1,
      scales = "free_y",
      strip.position = "top"
    ) +
    theme_minimal() +
    theme(
      legend.position  = "bottom",
      axis.text.y      = element_text(size = 9),
      strip.placement  = "outside",
      strip.background = element_blank(),
      strip.text.x     = element_text(size = 11, face = "bold"),
      plot.margin      = margin(5.5, 12, 5.5, 5.5)
    )
  
  print(p)
  ggsave(file_png, p, width = 18, height = 18, dpi = 300, bg = "white")
  if (!is.null(file_pdf)) ggsave(file_pdf, p, width = 18, height = 18, bg = "white")
}

set1 <- c("non_skin_cancer_mortality",
          "cancer_mortality",
          "cvd_mortality_narrow",
          "all_cause_mortality")

set2 <- c("non_melanoma_mortality",
          "melanoma_skin_cancer_mortality")

make_forest(plot_df, set1,
            "5_4_mortality_forest_BEEM_stepwise_set1.png",
            "5_4_mortality_forest_BEEM_stepwise_set1.pdf")

make_forest(plot_df, set2,
            "5_4_mortality_forest_BEEM_stepwise_set2.png",
            "5_4_mortality_forest_BEEM_stepwise_set2.pdf")

make_forest_shape(
  plot_df, set1,
  "5_4_mortality_forest_BEEM_stepwise_set1_shape.png",
  "5_4_mortality_forest_BEEM_stepwise_set1_shape.pdf"
)

make_forest_shape(
  plot_df, set2,
  "5_4_mortality_forest_BEEM_stepwise_set2_shape.png",
  "5_4_mortality_forest_BEEM_stepwise_set2_shape.pdf"
)

combined_set <- c(
  "all_cause_mortality",
  "cvd_mortality_narrow",
  "cancer_mortality",
  "non_skin_cancer_mortality",
  "melanoma_skin_cancer_mortality",
  "non_melanoma_mortality"
)

make_forest_shape(
  plot_df,
  combined_set,
  "5_4_mortality_forest_BEEM_stepwise_combined_shape.png",
  "5_4_mortality_forest_BEEM_stepwise_combined_shape.pdf"
)

make_forest(
  plot_df,
  combined_set,
  "5_4_mortality_forest_BEEM_stepwise_combined.png",
  "5_4_mortality_forest_BEEM_stepwise_combined.pdf"
)
