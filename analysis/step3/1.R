if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lubridate)
library(mice)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step3/"))

getwd()

library(tidyverse)
library(lubridate)

hos_de_ca <- read_csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "raw_data_3_8", "hos_de_ca_3_8.csv"))
date_assess <- read_csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "raw_data_3_8", "socia_field_columns.csv")) %>% 
  select(eid, p53_i0)

date_an <- read_csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "raw_data_3_8", "request_an_columns.csv"))

lost_eids <- date_an %>%
  filter(!is.na(p191)) %>%
  pull(eid)

head(lost_eids)

table(hos_de_ca$p41280_a0)

table(hos_de_ca$p41280_a1)
write.table(lost_eids, "../washed_data/3_1_0_eid_lost_to_followup.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

p41280_cols <- names(hos_de_ca)[str_detect(names(hos_de_ca), "^p41280_a")]

p40006_cols <- names(hos_de_ca)[str_detect(names(hos_de_ca), "^p40006")]
p40005_cols <- names(hos_de_ca)[str_detect(names(hos_de_ca), "^p40005")]
p40011_cols <- names(hos_de_ca)[str_detect(names(hos_de_ca), "^p40011")]

outcomes <- hos_de_ca %>%
  full_join(date_assess, by = "eid") %>%
  select(
    eid,
    date_assess = p53_i0,
    diag_icd10 = p41270,
    primary_cause = p40001_i0,
    death_date = p40000_i0,
    all_of(p40006_cols),
    all_of(p40005_cols),
    all_of(p40011_cols),
    all_of(p41280_cols)  # 加入所有诊断日期列
  )

library(dplyr)
library(stringr)
library(tidyr)

hos_icd <- hos_de_ca %>%
  select(eid, icd_raw = p41270) %>%
  mutate(
    icd_raw = as.character(icd_raw),
    icd_raw = str_remove_all(icd_raw, "\\[|\\]|'")  # 清除中括号和引号
  )

max_icd <- hos_icd$icd_raw %>%
  replace_na("") %>%
  str_count(",") %>%
  max(na.rm = TRUE) + 1

hos_icd_split <- hos_icd %>%
  separate(
    icd_raw,
    into = paste0("icd10_", 0:(max_icd - 1)),
    sep = ",\\s*",
    fill = "right"
  )

hos_icd_long <- hos_icd_split %>%
  pivot_longer(
    cols = starts_with("icd10_"),
    names_to = "diag_index",
    names_pattern = "icd10_(\\d+)",
    values_to = "icd10"
  ) %>%
  mutate(
    diag_index = as.integer(diag_index),
    icd10 = str_trim(icd10)
  )

p41280_cols <- paste0("p41280_a", 0:(max_icd - 1))
outcomes_long_date <- hos_de_ca %>%
  select(eid, all_of(p41280_cols)) %>%
  pivot_longer(
    cols = all_of(p41280_cols),
    names_to = "date_index",
    names_pattern = "p41280_a(\\d+)",
    values_to = "diag_date"
  ) %>%
  mutate(date_index = as.integer(date_index))

matched_outcomes <- hos_icd_long %>%
  left_join(outcomes_long_date,
            by = c("eid" = "eid", "diag_index" = "date_index")) %>%
  filter(!is.na(icd10))  # 去除空 ICD10

cancer_date_long <- hos_de_ca %>%
  select(eid, all_of(p40005_cols)) %>%
  pivot_longer(
    cols = starts_with("p40005"),
    names_to = "index",
    names_pattern = "p40005_i(\\d+)",
    values_to = "cancer_date"
  )

table(cancer_date_long$index)

cancer_icd10_long <- hos_de_ca %>%
  select(eid, all_of(p40006_cols)) %>%
  pivot_longer(
    cols = starts_with("p40006"),
    names_to = "index",
    names_pattern = "p40006_i(\\d+)",
    values_to = "cancer_icd10"
  )

histology_long <- hos_de_ca %>%
  select(eid, all_of(p40011_cols)) %>%
  pivot_longer(
    cols = starts_with("p40011"),
    names_to = "index",
    names_pattern = "p40011_i(\\d+)",
    values_to = "histology"
  )

cancer_long <- cancer_date_long %>%
  full_join(cancer_icd10_long, by = c("eid", "index")) %>%
  full_join(histology_long, by = c("eid", "index")) %>%
  arrange(eid, index) %>%
  filter(!is.na(cancer_icd10))

outcomes <- hos_de_ca %>%
  inner_join(date_assess, by = "eid") %>%
  select(
    eid,
    date_assess = p53_i0,
    primary_cause = p40001_i0,
    death_date = p40000_i0 # 加入所有诊断日期列
  )

icd10_date_long <- matched_outcomes %>%
  select(eid, icd10, diag_date)

icd10_date_long_ass <- icd10_date_long %>%
  inner_join(date_assess, by = "eid")

cancer_long <- cancer_long %>%
  select(eid, cancer_date, cancer_icd10, histology)

cancer_long_ass <- cancer_long %>%
  inner_join(date_assess, by = "eid")

write_csv(outcomes, "../washed_data/3_1_1_death_data.csv")
write_csv(icd10_date_long_ass, "../washed_data/3_1_2_inpatient_data.csv")
write_csv(cancer_long_ass, "../washed_data/3_1_3_cancer_data.csv")
