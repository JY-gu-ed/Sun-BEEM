if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(dplyr)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step1/"))

getwd()

filtered_uv <- read.csv("../washed_data/1_recoded_all_uv_data.csv")
colnames(filtered_uv)

selected_uv <- filtered_uv %>%
  select(eid, outdoor_summer, outdoor_winter, uv_protection, solarium_use)

nrow(selected_uv)
selected_uv_clean <- na.omit(selected_uv)
nrow(selected_uv_clean)

swr_uv <- read_csv("../raw_data_3_8/Dataset_neid_SWR.csv")
names(swr_uv)[names(swr_uv) == "n_eid"] <- "eid"
head(names(swr_uv))

swr_uv <- swr_uv[, c("eid", "swr2008_kJ")]

swr_uv_clean <- na.omit(swr_uv)
nrow(swr_uv)        # 原始行数
nrow(swr_uv_clean)  # 删除 NA 后的行数

merged_data <- merge(selected_uv_clean, swr_uv_clean, by = "eid")
nrow(merged_data)

if (!dir.exists("../washed_data")) {
  dir.create("../washed_data")
}

write.csv(merged_data, file = "../washed_data/1_2_merged_uv_data.csv", row.names = FALSE)
