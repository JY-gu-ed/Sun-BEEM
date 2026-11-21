if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())
if (!nzchar(Sys.getenv("SUN_BEEM_CODE_ROOT"))) stop("Set SUN_BEEM_CODE_ROOT or use tools/run_analysis.R.")
source(file.path(Sys.getenv("SUN_BEEM_CODE_ROOT"), "analysis/shared/followup.R"))


library(tidyverse)
library(survival)
library(broom)
library(purrr)
library(ggplot2)
setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step7/"))

outdir <- "../result/stepwise/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

prep_cancer_data <- function(df, status_col, after_col, date_col="cancer_date") {
  endtime <- registry_end_date("cancer")
  
  df <- apply_incidence_followup(df, status_col, date_col, endtime, baseline_col = "date_assess")
df <- df %>%
    filter(is.na(.data[[after_col]]) | .data[[after_col]] != "no") %>%
    mutate(
      diag_time = diag_time,
      status_bin = coalesce(as.integer(.data[[status_col]] == "yes"), 0L),
      BEEM = factor(BEEM, levels=c("low","medium","high")),
      BEEM_numeric = as.numeric(BEEM),  # low=1, medium=2, high=3
      sex  = factor(sex, levels=c("female","male")),
      age  = factor(age, levels=c("37-47","48-54","55-59","60-63","64-73")),
      bmi  = factor(bmi, levels=c("underweight","healthy weight","overweight","obese")),
      sleep_score = factor(sleep_score, levels=c("poor","intermediate","healthy")),
      drinking    = factor(drinking, levels=c("unhealthy","healthy")),
      phys_activity = factor(phys_activity, levels=c("low","moderate","high")),
      highest_edu   = factor(highest_edu, levels=c("less than secondary","secondary school","university or above")),
      smoking_status = factor(smoking_status, levels=c("current","previous","never")),
      deprivation_quintile = factor(deprivation_quintile, levels=1:5, ordered=TRUE)
    )
  
  return(df)
}

covariates <- c("sex","age","bmi","deprivation_quintile",
                "highest_edu","phys_activity","sleep_score",
                "smoking_status","drinking")

safe_fit <- function(formula, data) {
  tryCatch(coxph(formula, data = data, na.action = na.exclude), error = function(e) NULL)
}

run_one_incidence <- function(df, time, status, outcome_name) {
  results <- list()
  
  f_crude <- as.formula(paste0("Surv(", time, ", ", status, ") ~ BEEM"))
  fit_crude <- safe_fit(f_crude, df)
  
  crude_df <- if (!is.null(fit_crude)) {
    tidy(fit_crude, exponentiate=TRUE, conf.int=TRUE) %>%
      filter(term %in% c("BEEMmedium","BEEMhigh")) %>%
      mutate(model="Model 1: crude") %>%
      select(model, term, estimate, conf.low, conf.high, p.value)
  } else tibble(model="Model 1: crude", term=c("BEEMmedium","BEEMhigh"),
                estimate=NA, conf.low=NA, conf.high=NA, p.value=NA)
  
  f_crude_trend <- as.formula(paste0("Surv(", time, ", ", status, ") ~ BEEM_numeric"))
  fit_crude_trend <- safe_fit(f_crude_trend, df)
  p_trend_crude <- if (!is.null(fit_crude_trend))
    tidy(fit_crude_trend)$p.value[ tidy(fit_crude_trend)$term=="BEEM_numeric" ] else NA
  
  crude_df <- crude_df %>%
    mutate(outcome=outcome_name,
           HR_CI=sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
           p=sprintf("%.4f", p.value),
           p_trend=sprintf("%.4f", p_trend_crude))
  
  results[[1]] <- crude_df
  
  for (i in seq_along(covariates)) {
    rhs <- paste(c("BEEM", covariates[1:i]), collapse=" + ")
    f <- as.formula(paste0("Surv(", time, ", ", status, ") ~ ", rhs))
    fit <- safe_fit(f, df)
    
    adj_df <- if (!is.null(fit)) {
      tidy(fit, exponentiate=TRUE, conf.int=TRUE) %>%
        filter(term %in% c("BEEMmedium","BEEMhigh")) %>%
        mutate(model=paste0("Model ", i+1, ": crude + ", paste(covariates[1:i], collapse=" + "))) %>%
        select(model, term, estimate, conf.low, conf.high, p.value)
    } else tibble(model=paste0("Model ", i+1, ": crude + ", paste(covariates[1:i], collapse=" + ")),
                  term=c("BEEMmedium","BEEMhigh"),
                  estimate=NA, conf.low=NA, conf.high=NA, p.value=NA)
    
    f_trend <- as.formula(paste0("Surv(", time, ", ", status, ") ~ BEEM_numeric + ",
                                 paste(covariates[1:i], collapse=" + ")))
    fit_trend <- safe_fit(f_trend, df)
    p_trend <- if (!is.null(fit_trend))
      tidy(fit_trend)$p.value[ tidy(fit_trend)$term=="BEEM_numeric" ] else NA
    
    adj_df <- adj_df %>%
      mutate(outcome=outcome_name,
             HR_CI=sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
             p=sprintf("%.4f", p.value),
             p_trend=sprintf("%.4f", p_trend))
    
    results[[i+1]] <- adj_df
  }
  
  bind_rows(results)
}

outcomes <- list(
  all_cancer       = list(file="4_2_4_ana_for_allcancer_reg.csv", status="cancer_status", after="diagnosed_after_assessment"),
  skin_cancer      = list(file="4_2_5_ana_for_skin_cancer_reg.csv", status="skin_cancer_status", after="skin_cancer_after_assessment"),
  melanoma         = list(file="4_2_6_ana_for_melanoma_reg.csv", status="melanoma_status", after="melanoma_after_assessment"),
  other_skin_cancer= list(file="4_2_7_ana_for_other_skin_cancer_reg.csv", status="other_skin_cancer_status", after="other_skin_cancer_after_assessment"),
  non_skin_cancer  = list(file="4_2_8_ana_for_non_skin_cancer_reg.csv", status="non_skin_cancer_status", after="non_skin_cancer_after_assessment")
)

final_tbl_inc <- map_dfr(names(outcomes), function(oc){
  info <- outcomes[[oc]]
  
  cancer <- read.csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "washed_data", info$file)) %>%
    prep_cancer_data(status_col=info$status, after_col=info$after)
  
  run_one_incidence(cancer, "diag_time", "status_bin", oc)
})

write.csv(final_tbl_inc, file.path(outdir, "cox_BEEM_stepwise_incidence.csv"), row.names=FALSE)

outcome_labels_inc <- c(
  all_cancer        = "All cancer incidence",
  skin_cancer       = "Skin cancer incidence",
  melanoma          = "Melanoma incidence",
  other_skin_cancer = "Non-melanoma skin cancer incidence",
  non_skin_cancer   = "Non-skin cancer incidence"
)

plot_df_inc <- final_tbl_inc %>%
  mutate(
    model = factor(model, levels = unique(model)),
    HR_text = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
    outcome_label = factor(outcome,
                           levels = names(outcome_labels_inc),
                           labels = unname(outcome_labels_inc[names(outcome_labels_inc)]))
  )

make_forest_inc <- function(df, outcome_order, file_png, file_pdf=NULL){
  df2 <- df %>%
    filter(outcome %in% outcome_order) %>%
    mutate(
      outcome_label = factor(outcome,
                             levels = outcome_order,
                             labels = unname(outcome_labels_inc[outcome_order]))
    )
  
  pd <- position_dodge(width=0.8)
  
  x_min <- min(df2$conf.low, na.rm=TRUE) / 1.5
  x_max <- max(df2$conf.high, na.rm=TRUE) * 1.5
  
  p <- ggplot(df2, aes(y=model, x=estimate, xmin=conf.low, xmax=conf.high, color=term)) +
    geom_vline(xintercept=1, linetype=2) +
    geom_errorbarh(height=0.2, position=pd) +
    geom_point(size=2.5, position=pd) +
    geom_text(aes(x=conf.high, label=HR_text),
              hjust=-0.1, vjust=0.5, size=3, position=pd, show.legend=FALSE) +
    scale_x_log10(limits=c(x_min, x_max)) +
    scale_color_manual(values=c("BEEMmedium"="#1f78b4","BEEMhigh"="#e31a1c"),
                       labels=c("BEEMmedium"="Medium vs Low", "BEEMhigh"="High vs Low")) +
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
  ggsave(file_png, p, width=12, height=10, dpi=300, bg="white")
  if(!is.null(file_pdf)) ggsave(file_pdf, p, width=12, height=10, bg="white")
}

make_forest_inc_shape <- function(df, outcome_order, file_png, file_pdf = NULL){
  df2 <- df %>%
    filter(outcome %in% outcome_order) %>%
    mutate(
      outcome_label = factor(
        outcome,
        levels = outcome_order,
        labels = unname(outcome_labels_inc[outcome_order])
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
        BEEMmedium = 16,   # 圆点 = medium vs low
        BEEMhigh   = 15    # 方块 = high vs low
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
    facet_grid(outcome_label ~ ., scales = "free_y", space = "free_y") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.y     = element_text(size = 9),
      strip.text.y    = element_text(size = 11, face = "bold"),
      plot.margin     = margin(5.5, 12, 5.5, 5.5)
    )
  
  print(p)
  ggsave(file_png, p, width = 12, height = 10, dpi = 300, bg = "white")
  if (!is.null(file_pdf)) ggsave(file_pdf, p, width = 12, height = 10, bg = "white")
}

set1_inc <- c("all_cancer","non_skin_cancer","skin_cancer")
set2_inc <- c("melanoma","other_skin_cancer")

make_forest_inc(plot_df_inc, set1_inc,
                file.path(outdir, "5_4_incidence_forest_BEEM_stepwise_set1.png"),
                file.path(outdir, "5_4_incidence_forest_BEEM_stepwise_set1.pdf"))

make_forest_inc(plot_df_inc, set2_inc,
                file.path(outdir, "5_4_incidence_forest_BEEM_stepwise_set2.png"),
                file.path(outdir, "5_4_incidence_forest_BEEM_stepwise_set2.pdf"))

make_forest_inc_shape(
  plot_df_inc, set1_inc,
  file.path(outdir, "5_4_incidence_forest_BEEM_stepwise_set1_shape.png"),
  file.path(outdir, "5_4_incidence_forest_BEEM_stepwise_set1_shape.pdf")
)

make_forest_inc_shape(
  plot_df_inc, set2_inc,
  file.path(outdir, "5_4_incidence_forest_BEEM_stepwise_set2_shape.png"),
  file.path(outdir, "5_4_incidence_forest_BEEM_stepwise_set2_shape.pdf")
)

cvd_step <- read.csv("../result/6_1_cvd_cox_BEEM_stepwise.csv")

plot_df_cvd <- cvd_step %>%
  dplyr::select(-p_trend) %>%
  mutate(
    model   = factor(model, levels = unique(model)),
    HR_text = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
    outcome = "cvd_incidence"   # 新的 outcome 名，用来和癌症区分
  )

outcome_labels_all <- c(
  cvd_incidence     = "Cardiovascular disease incidence",
  all_cancer        = "All cancer incidence",
  non_skin_cancer   = "Non-skin cancer incidence",
  melanoma          = "Melanoma incidence",
  other_skin_cancer = "Other skin cancer incidence"
)

plot_df_all <- bind_rows(
  plot_df_cvd,
  plot_df_inc %>%
    dplyr::select(-p, -p_trend) %>%
    dplyr::filter(outcome != "skin_cancer")   # 不合并 skin cancer 这一栏
) %>%
  mutate(
    outcome_label = factor(
      outcome,
      levels = names(outcome_labels_all),
      labels = unname(outcome_labels_all[names(outcome_labels_all)])
    ),
    model = factor(model, levels = unique(model))
  )

sapply(plot_df_cvd[, c("term","estimate","conf.low","conf.high","HR_text","outcome")], class)

sapply(plot_df_inc[, c("term","estimate","conf.low","conf.high","HR_text","outcome")], class)

plot_df_all <- bind_rows(
  plot_df_cvd,
  plot_df_inc %>% dplyr::select(-p, -p_trend)
) %>%
  mutate(
    outcome_label = factor(
      outcome,
      levels = names(outcome_labels_all),
      labels = unname(outcome_labels_all[names(outcome_labels_all)])
    ),
    model = factor(model, levels = unique(model))
  )

str(plot_df_all)

dplyr::count(plot_df_all, outcome, outcome_label)

head(plot_df_all)

make_forest_inc_cvd_shape <- function(df, file_png, file_pdf = NULL){
  pd <- position_dodge(width = 0.8)
  
  x_min <- min(df$conf.low, na.rm = TRUE) / 1.5
  x_max <- max(df$conf.high, na.rm = TRUE) * 1.5
  
  p <- ggplot(
    df,
    aes(
      y     = model,
      x     = estimate,
      xmin  = conf.low,
      xmax  = conf.high,
      shape = term       # 用形状区分 BEEM 组别
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
        BEEMmedium = 16,  # 圆点
        BEEMhigh   = 15   # 方块
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
  ggsave(file_png, p, width = 12, height = 18, dpi = 300, bg = "white")
  if (!is.null(file_pdf)) ggsave(file_pdf, p, width = 12, height = 18, bg = "white")
}

make_forest_inc_cvd_shape(
  plot_df_all,
  file.path(outdir, "5_4_incidence_cvd_forest_BEEM_stepwise_combined_shape.png"),
  file.path(outdir, "5_4_incidence_cvd_forest_BEEM_stepwise_combined_shape.pdf")
)
