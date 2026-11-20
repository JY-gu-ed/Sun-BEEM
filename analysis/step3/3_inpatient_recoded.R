if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(here)
library(ggplot2)
library(tidyverse)
library(here)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(here)
library(lubridate)
library(mice)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step3/"))

getwd()

outcomes <- read.csv("../washed_data/3_1_2_inpatient_data.csv")

cvd_cases <- outcomes %>%
  mutate(
    icd10 = replace_na(icd10, "no"),
    diag_date = lubridate::ymd(diag_date),
    assessment_date = lubridate::ymd(p53_i0),
    cvd_mortality_narrow = ifelse(str_detect(icd10, "I21|I22|I23|I24|I25|I60|I61|I62|I63|I64"), "yes", "no")
  ) %>%
  filter(cvd_mortality_narrow == "yes") %>%
  mutate(
    diagnosed_after_assessment = ifelse(diag_date > assessment_date, "yes", "no")
  )

n_records_after_assessment <- cvd_cases %>%
  filter(diagnosed_after_assessment == "yes") %>%
  nrow()

n_patients_after_assessment <- cvd_cases %>%
  filter(diagnosed_after_assessment == "yes") %>%
  distinct(eid) %>%
  nrow()

cat("诊断在体检后的记录数:", n_records_after_assessment, "\n")
cat("诊断在体检后的唯一患者数:", n_patients_after_assessment, "\n")

write.csv(cvd_cases, "../washed_data/3_3_recoded_inpatient_cvd.csv", row.names = FALSE)
