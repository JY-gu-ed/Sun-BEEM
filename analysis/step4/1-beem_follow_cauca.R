if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(here)
library(ggplot2)
library(tidyverse)
library(here)
library(ggplot2)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step4/"))

getwd()

beem <- read.csv(here("washed_data", "1_3_bi_merged_uv.csv"))
nrow(beem)

lost_to_followup <- read.table(here("washed_data", "3_1_0_eid_lost_to_followup.txt"),
                               header = FALSE,  # 没有列名
                               col.names = "eid") # 自定义列名为 eid

beem <- beem %>%
  anti_join(lost_to_followup, by = "eid")
nrow(beem)

cau <- read.csv(here("washed_data", "1_6_ances_data.csv"))
beem <- beem %>%
  semi_join(cau %>% select(eid), by = "eid")

nrow(beem)

readr::write_csv(beem, "../washed_data/4_1_beem_follow_cauca.csv")

summary(beem$swr2008_kJ)
quantile(beem$swr2008_kJ, probs = 0.5, na.rm = TRUE)
quantile(beem$swr2008_kJ, probs = c(1/3, 2/3), na.rm = TRUE)
quantile(beem$swr2008_kJ, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

summary(beem$outdoor_avg)
