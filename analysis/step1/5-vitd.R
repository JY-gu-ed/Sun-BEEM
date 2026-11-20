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

date_an <- read_csv(here("raw_data_3_8", "request_an_columns.csv"))

vitd_data <- date_an %>%
  select(eid, p30890_i0)

write.csv(vitd_data, "../washed_data/1_5_vitd_child_data.csv", row.names = FALSE)
