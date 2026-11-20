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

filtered_uv <- read.csv("../washed_data/1_1_recoded_all_uv_data.csv")
sunburn_data <- filtered_uv %>%
  select(eid, sunburn_child)

write.csv(sunburn_data, "../washed_data/1_4_sunburn_child_data.csv", row.names = FALSE)
