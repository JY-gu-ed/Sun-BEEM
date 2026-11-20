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

morta <- read.csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "washed_data", "4_2_2_ana_for_mortality.csv"))

morta_backup<- morta 

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
  mutate(across(all_of(mortality_vars), ~ ifelse(. == "yes", 1, 0)))

acute <- read.csv("daily_uva_temp.csv")
sum(duplicated(acute$n_eid))

acute_first <- acute %>%
  arrange(n_eid, year) %>%   # 确保“第一”的定义
  distinct(n_eid, .keep_all = TRUE)

acute2 <- read.csv("UVA_UKB_allethn.csv")

sum(duplicated(acute2$n_eid))
acute2 %>%
  filter(duplicated(n_eid) | duplicated(n_eid, fromLast = TRUE)) %>%
  arrange(n_eid)

acute_unique <- acute2 %>%
  distinct(n_eid, .keep_all = TRUE)

library(dplyr)

morta <- morta %>%
  left_join(
    acute_first %>% select(n_eid, uva_2dayavg),
    by = c("eid" = "n_eid")
  ) %>%
  left_join(
    acute_unique %>% select(n_eid, uva_7_dayavg, uva_14_dayavg, uva_30_dayavg),
    by = c("eid" = "n_eid")
  )

olink <- read_csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "raw_data_3_8", "olink.csv"))

olink <- olink %>% 
  rename(eid = olink_instance_0.eid)

colnames(olink) <- gsub("olink_instance_0\\.", "", colnames(olink))

biomarkers <- setdiff(names(olink), "eid")
saveRDS(biomarkers, file = "biomarkers_olink.rds")

morta <- morta %>% 
  inner_join(olink, "eid")

saveRDS(morta, file = "mea_set_merged.rds")
