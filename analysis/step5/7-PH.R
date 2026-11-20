if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(survival)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step5/"))
getwd()

morta <- read.csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "washed_data", "4_2_2_ana_for_mortality.csv"))

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
  mutate(across(all_of(mortality_vars), ~ as.numeric(ifelse(. == "yes", 1, 0))))

covariates <- c("BEEM","sex","age","bmi","deprivation_quintile",
                "highest_edu","phys_activity","sleep_score",
                "smoking_status","drinking")

ph_out_dir <- "result/ph_tests/"
if (!dir.exists(ph_out_dir)) dir.create(ph_out_dir, recursive = TRUE)

run_ph <- function(outcome, data = morta) {
  message("▶ 正在检验 outcome: ", outcome)
  rhs <- paste(covariates, collapse = " + ")
  f <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ ", rhs))
  
  fit <- tryCatch(coxph(f, data = data, na.action = na.exclude),
                  error = function(e) { 
                    message("❌ 出错: ", outcome, " -- ", e$message)
                    return(NULL)
                  })
  
  if (is.null(fit)) return(NULL)
  
  ph <- tryCatch(cox.zph(fit),
                 error = function(e) {
                   message("❌ PH 检验失败: ", outcome, " -- ", e$message)
                   return(NULL)
                 })
  
  if (is.null(ph)) return(NULL)
  
  ph_df <- as.data.frame(ph$table)
  ph_df$term <- rownames(ph$table)
  ph_df$outcome <- outcome
  return(ph_df)
}

ph_results <- lapply(mortality_vars, run_ph)

ph_summary <- bind_rows(ph_results)
write.csv(ph_summary, file.path(ph_out_dir, "ph_test_summary.csv"), row.names = FALSE)

ph_summary
