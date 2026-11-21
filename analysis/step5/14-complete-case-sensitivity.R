if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT.")
suppressPackageStartupMessages({library(data.table); library(survival); library(broom)})
root <- Sys.getenv("SUN_BEEM_PROJECT_ROOT")
d <- fread(file.path(root, "washed_data/4_2_2_ana_for_mortality.csv"))
membership <- fread(file.path(root, "washed_data/2_1_0_covariate_complete_cases.csv"))
stopifnot(!anyDuplicated(membership$eid))
d <- d[eid %in% membership[core_complete_case == TRUE, eid]]
d$BEEM <- factor(d$BEEM, levels = c("low", "medium", "high"))
covariates <- c("age", "sex", "bmi", "deprivation_quintile", "highest_edu",
                "smoking_status", "drinking", "phys_activity", "sleep_score")
for (v in covariates) d[[v]] <- factor(d[[v]])
d$deprivation_quintile <- ordered(d$deprivation_quintile, levels = 1:5)
outcomes <- c("all_cause_mortality", "cvd_mortality_narrow", "cancer_mortality",
              "skin_cancer_mortality", "melanoma_skin_cancer_mortality",
              "non_melanoma_mortality", "non_skin_cancer_mortality")
results <- list()
for (outcome in outcomes) {
  d[[outcome]] <- as.integer(as.character(d[[outcome]]) %in% c("yes", "1"))
  if (sum(d[[outcome]]) == 0) next
  formula <- as.formula(paste0("Surv(death_time, ", outcome, ") ~ BEEM + ", paste(covariates, collapse = " + ")))
  fit <- coxph(formula, data = d, ties = "efron")
  row <- as.data.table(broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE))[grepl("^BEEM", term)]
  row[, `:=`(outcome = outcome, N = fit$n, events = fit$nevent)]
  results[[outcome]] <- row
}
out <- file.path(root, "result/complete_case")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(rbindlist(results), file.path(out, "adjusted_mortality.csv"))
