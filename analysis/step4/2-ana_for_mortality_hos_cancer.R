if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(tidyverse)
library(ggplot2)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step4/"))

getwd()

beem <- read.csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "washed_data", "4_1_beem_follow_cauca.csv"))
social <- read.csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "washed_data", "2_1_2_social_mice_imputed.csv"))
sunburn <-  read.csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "washed_data", "1_4_sunburn_child_data.csv"))
vitd  <-  read.csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "washed_data", "1_5_vitd_child_data.csv"))
morta <- read.csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "washed_data", "3_2_recoded_death.csv"))
hos <- read.csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "washed_data", "3_3_recoded_inpatient_cvd.csv"))

followup_dates <- morta %>% select(eid, date_assess, death_date)
loss_dates <- readr::read_csv(
  file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "raw_data_3_8", "request_an_columns.csv"),
  col_select = c(eid, p191), show_col_types = FALSE
) %>% rename(loss_to_followup_date = p191)
stopifnot(!anyDuplicated(followup_dates$eid), !anyDuplicated(loss_dates$eid))
followup_dates <- left_join(followup_dates, loss_dates, by = "eid")

cancer <- read.csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "washed_data", "3_4_recoded_cancer.csv"))

vitd <- vitd %>%
  rename(vitd = p30890_i0)
beem_vitd <- beem %>%
  left_join(vitd, by = "eid")
sum(is.na(beem_vitd $vitd))

beem_vitd_sunburn <- beem_vitd %>%
  left_join(sunburn, by = "eid")
sum(is.na(beem_vitd_sunburn$sunburn_child))

beem_vitd_sunburn_social <- beem_vitd_sunburn %>%
  left_join(social, by = "eid")

write_csv(beem_vitd_sunburn_social,
          "../washed_data/4_2_1_beem_vitd_sunburn_social.csv")

beem_morta <- beem_vitd_sunburn_social %>%
  inner_join(morta, by = "eid")

readr::write_csv(beem_morta, "../washed_data/4_2_2_ana_for_mortality.csv")

hos <- hos %>%
  mutate(diag_date = as.Date(diag_date)) %>%  # 转换为日期
  group_by(eid) %>%
  slice_min(order_by = diag_date, n = 1, with_ties = FALSE) %>% # 保留最早
  ungroup()

beem_hos <- beem_vitd_sunburn_social %>%
  left_join(hos %>% select(-p53_i0, -assessment_date), by = "eid")

beem_hos <- beem_hos %>%
  left_join(followup_dates, by = "eid")

readr::write_csv(beem_hos, "../washed_data/4_2_3_ana_for_cvd_inpatient.csv")

cancer_sub <- cancer %>%
  select(eid, cancer_date, cancer_icd10, histology,
         cancer_status, diagnosed_after_assessment) %>%
  mutate(cancer_date = as.Date(cancer_date)) %>%
  group_by(eid) %>%
  slice_min(order_by = cancer_date, n = 1, with_ties = FALSE) %>%
  ungroup()

beem_cancer <- beem_vitd_sunburn_social %>%
  left_join(cancer_sub, by = "eid")

beem_cancer <- beem_cancer %>%
  left_join(followup_dates, by = "eid")

readr::write_csv(beem_cancer, "../washed_data/4_2_4_ana_for_allcancer_reg.csv")

library(dplyr)

skin_cancer_sub <- cancer %>%
  filter(skin_cancer_status == "yes") %>%
  select(eid, cancer_date, cancer_icd10, histology,
         skin_cancer_status, skin_cancer_after_assessment) %>%
  mutate(cancer_date = as.Date(cancer_date)) %>%
  group_by(eid) %>%
  slice_min(order_by = cancer_date, n = 1, with_ties = FALSE) %>%
  ungroup()

beem_skin_cancer_sub <- beem_vitd_sunburn_social %>%
  left_join(skin_cancer_sub, by = "eid")

beem_skin_cancer_sub <- beem_skin_cancer_sub %>%
  left_join(followup_dates, by = "eid")

readr::write_csv(beem_skin_cancer_sub, "../washed_data/4_2_5_ana_for_skin_cancer_reg.csv")

melanoma_sub <- cancer %>%
  filter(melanoma_status == "yes") %>%
  select(eid, cancer_date, cancer_icd10, histology,
         melanoma_status, melanoma_after_assessment) %>%
  mutate(cancer_date = as.Date(cancer_date)) %>%
  group_by(eid) %>%
  slice_min(order_by = cancer_date, n = 1, with_ties = FALSE) %>%
  ungroup()

beem_melanoma_sub <- beem_vitd_sunburn_social %>%
  left_join(melanoma_sub, by = "eid") %>%
  left_join(followup_dates, by = "eid")

readr::write_csv(beem_melanoma_sub,
                 "../washed_data/4_2_6_ana_for_melanoma_reg.csv")

library(dplyr)

other_skin_cancer_sub <- cancer %>%
  filter(other_skin_cancer_status == "yes") %>%
  select(eid, cancer_date, cancer_icd10, histology,
         other_skin_cancer_status, other_skin_cancer_after_assessment) %>%
  mutate(cancer_date = as.Date(cancer_date)) %>%
  group_by(eid) %>%
  slice_min(order_by = cancer_date, n = 1, with_ties = FALSE) %>%
  ungroup()

beem_other_skin_cancer_sub <- beem_vitd_sunburn_social %>%
  left_join(other_skin_cancer_sub, by = "eid") %>%
  left_join(followup_dates, by = "eid")

readr::write_csv(beem_other_skin_cancer_sub,
                 "../washed_data/4_2_7_ana_for_other_skin_cancer_reg.csv")

non_skin_cancer_sub <- cancer %>%
  filter(non_skin_cancer_status == "yes") %>%
  select(eid, cancer_date, cancer_icd10, histology,
         non_skin_cancer_status, non_skin_cancer_after_assessment) %>%
  mutate(cancer_date = as.Date(cancer_date)) %>%
  group_by(eid) %>%
  slice_min(order_by = cancer_date, n = 1, with_ties = FALSE) %>%
  ungroup()

beem_non_skin_cancer_sub <- beem_vitd_sunburn_social %>%
  left_join(non_skin_cancer_sub, by = "eid") %>%
  left_join(followup_dates, by = "eid")

readr::write_csv(beem_non_skin_cancer_sub,
                 "../washed_data/4_2_8_ana_for_non_skin_cancer_reg.csv")
