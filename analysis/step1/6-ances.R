if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(here)
library(ggplot2)
library(tidyverse)
library(here)
library(ggplot2)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step1/"))

getwd()

data_ances <- read_csv(here("raw_data_3_8", "ancestry_id_UKB.csv"))

table(data_ances$ancestry)

library(dplyr)

data_ances <- data_ances %>%
  rename(eid = n_eid) %>%       # 重命名列
  mutate(Caucasian = "yes")     # 新增列

write.csv(data_ances, "../washed_data/1_6_ances_data.csv", row.names = FALSE)
