if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(here)
library(dplyr)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), ""))  # 先到父目录
if (!dir.exists("step10")) dir.create("step10", recursive = TRUE)
setwd("step10")
getwd()

morta <- readr::read_csv(here("washed_data", "4_2_2_ana_for_mortality.csv"), show_col_types = FALSE)
morta_backup <- morta

morta <- morta %>%
  mutate(
    BEEM = factor(BEEM, levels = c("low", "medium", "high")),
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

acute <- readr::read_csv(here("step8", "daily_uva_temp.csv"), show_col_types = FALSE)
sum(duplicated(acute$n_eid))

acute_first <- acute %>%
  arrange(n_eid, year) %>%                 # 确保“第一”的定义
  distinct(n_eid, .keep_all = TRUE)

acute2 <- readr::read_csv(here("step8", "UVA_UKB_allethn.csv"), show_col_types = FALSE)

sum(duplicated(acute2$n_eid))
acute2 %>%
  filter(duplicated(n_eid) | duplicated(n_eid, fromLast = TRUE)) %>%
  arrange(n_eid)

acute_unique <- acute2 %>%
  distinct(n_eid, .keep_all = TRUE)

morta <- morta %>%
  left_join(
    acute_first %>% select(n_eid, uva_2dayavg),
    by = c("eid" = "n_eid")
  ) %>%
  left_join(
    acute_unique %>% dplyr::select(n_eid, uva_7_dayavg, uva_14_dayavg, uva_30_dayavg),
    by = c("eid" = "n_eid")
  )

metabolomics_nmr <- readr::read_csv(
  here("washed_data", "1_7_metabolomics_nmr.csv"),
  show_col_types = FALSE
)

dim(metabolomics_nmr)
head(metabolomics_nmr[, 1:6])

metabolomics_nmr <- metabolomics_nmr %>%
  filter(
    rowSums(!is.na(dplyr::select(., -eid))) > 0
  )

dim(metabolomics_nmr)

biomarkers_nmr <- setdiff(names(metabolomics_nmr), "eid")
mea_nmr_set <- inner_join(morta, metabolomics_nmr, by = "eid")

saveRDS(mea_nmr_set, file = "mea_set_merged.rds")

biomarkers_nmr <- setdiff(names(metabolomics_nmr), "eid")

saveRDS(biomarkers_nmr, file = "biomarkers_nmr.rds")
