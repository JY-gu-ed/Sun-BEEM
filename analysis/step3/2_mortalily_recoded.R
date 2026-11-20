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

outcomes <- read.csv("../washed_data/3_1_1_death_data.csv")

outcomes <- outcomes %>% 
  mutate(all_cause_mortality = ifelse(is.na(primary_cause), "no", "yes"))

endtime <- "2025-08-04" # most recent date code has been run
endtime <- as.Date(endtime)

outcomes <- outcomes %>%
  mutate(death_time = ifelse(all_cause_mortality == "yes", as.Date(death_date) - as.Date(date_assess),
                             ifelse(all_cause_mortality == "no", as.numeric(endtime) - as.numeric(as.Date(date_assess)), NA)))

outcomes <- outcomes %>%
  mutate(primary_cause = replace_na(primary_cause, "no")) %>% 
  mutate(cvd_mortality_narrow = 
           ifelse(str_detect(primary_cause, "I21|I22|I23|I24|I25|I60|I61|I62|I63|I64"), "yes", "no"))

table(outcomes$cvd_mortality_narrow)

outcomes <- outcomes %>% 
  mutate(cancer_mortality = 
           ifelse(str_detect(primary_cause, "C|D37|D38|D39|D40|D41|D42|D43|D44|D45|D46|D47|D48"), "yes", "no"))

table(outcomes$cancer_mortality)

outcomes <- outcomes %>% 
  mutate(skin_cancer_mortality = 
           ifelse(str_detect(primary_cause, "C43|C44"), "yes", "no"))

table(outcomes$skin_cancer_mortality)

outcomes <- outcomes %>% 
  mutate(melanoma_skin_cancer_mortality = 
           ifelse(str_detect(primary_cause, "C43"), "yes", "no"))

table(outcomes$melanoma_skin_cancer_mortality)

outcomes <- outcomes %>% 
  mutate(non_melanoma_mortality = 
           ifelse(str_detect(primary_cause, "C44"), "yes", "no"))

outcomes <- outcomes %>% 
  mutate(non_skin_cancer_mortality = 
           ifelse(cancer_mortality == "yes" & skin_cancer_mortality == "no", "yes", "no"))

top_causes <- outcomes %>%
  filter(all_cause_mortality == "yes") %>%
  count(primary_cause, sort = TRUE) %>%
  slice_max(n, n = 20)

print(top_causes)

write.csv(outcomes, file = "../washed_data/3_2_recoded_death.csv", row.names = FALSE)
