if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(survival)
library(broom)
library(purrr)
library(grid)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step6/"))

cvd <- read.csv("../washed_data/4_2_3_ana_for_cvd_inpatient.csv")

cvd <- cvd %>%
  filter(is.na(diagnosed_after_assessment) | diagnosed_after_assessment != "no")

endtime <- as.Date("2025-08-04")
cvd <- cvd %>%
  mutate(
    diag_time = ifelse(
      !is.na(cvd_mortality_narrow) & cvd_mortality_narrow == "yes",
      as.numeric(as.Date(diag_date) - as.Date(date_assess)),
      as.numeric(endtime) - as.numeric(as.Date(date_assess))
    ),
    cvd_mortality_narrow_bin = coalesce(as.integer(cvd_mortality_narrow == "yes"), 0L)
  ) %>%
  mutate(
    BEEM = factor(BEEM, levels = c("low","medium","high")),
    BEEM_numeric = as.numeric(BEEM)
  )

covariates <- c("sex","age","bmi","deprivation_quintile",
                "highest_edu","phys_activity","sleep_score",
                "smoking_status","drinking")

safe_fit <- function(formula, data) {
  tryCatch(coxph(formula, data = data, na.action = na.exclude),
           error = function(e) NULL)
}

run_stepwise <- function(outcome_bin, time_var, data) {
  results <- list()
  
  f_crude <- as.formula(paste0("Surv(", time_var, ", ", outcome_bin, ") ~ BEEM"))
  fit_crude <- safe_fit(f_crude, data)
  
  crude_df <- if (!is.null(fit_crude)) {
    tidy(fit_crude, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term %in% c("BEEMmedium","BEEMhigh")) %>%
      mutate(model = "Model 1: crude")
  } else tibble(term=c("BEEMmedium","BEEMhigh"),
                estimate=NA, conf.low=NA, conf.high=NA, p.value=NA,
                model="Model 1: crude")
  
  f_trend <- as.formula(paste0("Surv(", time_var, ", ", outcome_bin, ") ~ BEEM_numeric"))
  fit_trend <- safe_fit(f_trend, data)
  p_trend <- if (!is.null(fit_trend))
    tidy(fit_trend)$p.value[ tidy(fit_trend)$term=="BEEM_numeric" ] else NA
  
  crude_df <- crude_df %>%
    mutate(outcome=outcome_bin, p_trend=sprintf("%.4f", p_trend))
  results[[1]] <- crude_df
  
  for (i in seq_along(covariates)) {
    rhs <- paste(c("BEEM", covariates[1:i]), collapse=" + ")
    f <- as.formula(paste0("Surv(", time_var, ", ", outcome_bin, ") ~ ", rhs))
    fit <- safe_fit(f, data)
    
    adj_df <- if (!is.null(fit)) {
      tidy(fit, exponentiate=TRUE, conf.int=TRUE) %>%
        filter(term %in% c("BEEMmedium","BEEMhigh")) %>%
        mutate(model=paste0("Model ", i+1, ": crude + ", paste(covariates[1:i], collapse=" + ")))
    } else tibble(term=c("BEEMmedium","BEEMhigh"),
                  estimate=NA, conf.low=NA, conf.high=NA, p.value=NA,
                  model=paste0("Model ", i+1, ": crude + ", paste(covariates[1:i], collapse=" + ")))
    
    f_trend <- as.formula(paste0("Surv(", time_var, ", ", outcome_bin, ") ~ BEEM_numeric + ",
                                 paste(covariates[1:i], collapse=" + ")))
    fit_trend <- safe_fit(f_trend, data)
    p_trend <- if (!is.null(fit_trend))
      tidy(fit_trend)$p.value[ tidy(fit_trend)$term=="BEEM_numeric" ] else NA
    
    adj_df <- adj_df %>%
      mutate(outcome=outcome_bin, p_trend=sprintf("%.4f", p_trend))
    results[[i+1]] <- adj_df
  }
  
  bind_rows(results) %>%
    mutate(HR_CI = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
           p = sprintf("%.4f", p.value))
}

final_tbl <- run_stepwise("cvd_mortality_narrow_bin", "diag_time", cvd)
write.csv(final_tbl, "../result/6_1_cvd_cox_BEEM_stepwise.csv", row.names=FALSE)
head(final_tbl)

parse_hrci <- function(estimate, low, high){
  sprintf("%.2f (%.2f–%.2f)", estimate, low, high)
}

plot_df <- final_tbl %>%
  mutate(
    model = factor(model, levels=unique(model)),
    HR_text = parse_hrci(estimate, conf.low, conf.high),
    outcome_label = "CVD incidence"
  )

make_forest <- function(df, file_png, file_pdf=NULL){
  pd <- position_dodge(width=0.8)
  x_min <- min(df$conf.low, na.rm=TRUE) / 1.5
  x_max <- max(df$conf.high, na.rm=TRUE) * 1.5
  
  p <- ggplot(df, aes(y=model, x=estimate, xmin=conf.low, xmax=conf.high, color=term)) +
    geom_vline(xintercept=1, linetype=2) +
    geom_errorbarh(height=0.2, position=pd) +
    geom_point(size=2.5, position=pd) +
    geom_text(aes(x=conf.high, label=HR_text),
              hjust=-0.1, vjust=0.5, size=3, position=pd, show.legend=FALSE) +
    scale_x_log10(limits=c(x_min, x_max)) +
    scale_color_manual(values=c("BEEMmedium"="#1f78b4","BEEMhigh"="#e31a1c"),
                       labels=c("BEEMmedium"="Medium vs Low","BEEMhigh"="High vs Low")) +
    labs(x="Hazard Ratio (log scale)", y=NULL, color=NULL, title="CVD incidence") +
    theme_minimal() +
    theme(
      legend.position="bottom",
      axis.text.y=element_text(size=8),
      plot.margin=margin(5.5, 20, 5.5, 5.5)
    )
  
  print(p)
  ggsave(file_png, p, width=10, height=6, dpi=300, bg="white")
  if(!is.null(file_pdf)) ggsave(file_pdf, p, width=10, height=6, bg="white")
}

make_forest(plot_df, "../result/6_2_cvd_forest_stepwise.png",
            "../result/6_2_cvd_forest_stepwise.pdf")
