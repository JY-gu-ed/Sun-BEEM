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

outcomes <- read.csv("../washed_data/3_1_3_cancer_data.csv")

library(dplyr)
library(lubridate)

cancer_cases <- outcomes %>%
  mutate(
    cancer_status = ifelse(!is.na(cancer_icd10) & cancer_icd10 != "", "yes", "no"),
    cancer_date = ymd(cancer_date),
    assessment_date = ymd(p53_i0),
    diagnosed_after_assessment = ifelse(cancer_status == "yes" & cancer_date > assessment_date, "yes", "no")
  )

n_all_cancer <- cancer_cases %>% filter(cancer_status == "yes") %>% nrow()
n_post_assessment_cancer <- cancer_cases %>% filter(diagnosed_after_assessment == "yes") %>% nrow()
n_unique_post_assessment_cancer <- cancer_cases %>% filter(diagnosed_after_assessment == "yes") %>% distinct(eid) %>% nrow()

cat("有癌症记录的总条数:", n_all_cancer, "\n")
cat("癌症诊断在体检之后的条数:", n_post_assessment_cancer, "\n")
cat("癌症诊断在体检之后的唯一人数:", n_unique_post_assessment_cancer, "\n")

cancer_cases <- cancer_cases %>%
  mutate(
    skin_cancer_status = ifelse(str_detect(cancer_icd10, "C43|C44"), "yes", "no"),
    skin_cancer_after_assessment = ifelse(skin_cancer_status == "yes" & cancer_date > assessment_date, "yes", "no")
  )

n_skin_cancer <- cancer_cases %>%
  filter(skin_cancer_status == "yes") %>%
  nrow()

n_skin_cancer_after <- cancer_cases %>%
  filter(skin_cancer_after_assessment == "yes") %>%
  nrow()

n_skin_cancer_after_patients <- cancer_cases %>%
  filter(skin_cancer_after_assessment == "yes") %>%
  distinct(eid) %>%
  nrow()

cat("皮肤癌记录总数:", n_skin_cancer, "\n")
cat("皮肤癌诊断在体检后的记录数:", n_skin_cancer_after, "\n")
cat("皮肤癌诊断在体检后的唯一患者数:", n_skin_cancer_after_patients, "\n")

cancer_cases <- cancer_cases %>%
  mutate(
    melanoma_status = ifelse(str_detect(cancer_icd10, "C43"), "yes", "no"),
    melanoma_after_assessment = ifelse(melanoma_status == "yes" & cancer_date > assessment_date, "yes", "no")
  )

n_melanoma <- cancer_cases %>%
  filter(melanoma_status == "yes") %>%
  nrow()

n_melanoma_after <- cancer_cases %>%
  filter(melanoma_after_assessment == "yes") %>%
  nrow()

n_melanoma_after_patients <- cancer_cases %>%
  filter(melanoma_after_assessment == "yes") %>%
  distinct(eid) %>%
  nrow()

cat("Melanoma 记录总数:", n_melanoma, "\n")
cat("Melanoma 诊断在体检后的记录数:", n_melanoma_after, "\n")
cat("Melanoma 诊断在体检后的唯一患者数:", n_melanoma_after_patients, "\n")

cancer_cases <- cancer_cases %>%
  mutate(
    other_skin_cancer_status = ifelse(str_detect(cancer_icd10, "C44"), "yes", "no"),
    other_skin_cancer_after_assessment = ifelse(other_skin_cancer_status == "yes" & cancer_date > assessment_date, "yes", "no")
  )

n_other_skin_cancer <- cancer_cases %>%
  filter(other_skin_cancer_status == "yes") %>%
  nrow()

n_other_skin_cancer_after <- cancer_cases %>%
  filter(other_skin_cancer_after_assessment == "yes") %>%
  nrow()

n_other_skin_cancer_after_patients <- cancer_cases %>%
  filter(other_skin_cancer_after_assessment == "yes") %>%
  distinct(eid) %>%
  nrow()

cat("非黑色素瘤皮肤癌记录总数:", n_other_skin_cancer, "\n")
cat("非黑色素瘤体检后诊断记录数:", n_other_skin_cancer_after, "\n")
cat("非黑色素瘤体检后诊断的唯一患者数:", n_other_skin_cancer_after_patients, "\n")

cancer_cases <- cancer_cases %>%
  filter(!is.na(cancer_icd10) & cancer_icd10 != "") %>%
  mutate(
    cancer_date = ymd(cancer_date),
    assessment_date = ymd(p53_i0),
    
    non_skin_cancer_status = ifelse(!str_detect(cancer_icd10, "C43|C44"), "yes", "no"),
    
    non_skin_cancer_after_assessment = ifelse(non_skin_cancer_status == "yes" & cancer_date > assessment_date, "yes", "no")
  )

cat("非皮肤癌记录总数:", sum(cancer_cases$non_skin_cancer_status == "yes"), "\n")
cat("非皮肤癌体检后诊断的记录数:", sum(cancer_cases$non_skin_cancer_after_assessment == "yes"), "\n")
cat("非皮肤癌体检后诊断的唯一患者数:", cancer_cases %>%
      filter(non_skin_cancer_after_assessment == "yes") %>%
      distinct(eid) %>% nrow(), "\n")

write.csv(cancer_cases, "../washed_data/3_4_recoded_cancer.csv", row.names = FALSE)
