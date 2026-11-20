if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(tidyverse)
library(ggplot2)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step1/"))

getwd()

library(tidyverse)
library(ggplot2)

raw_uv <- read_csv("../raw_data_3_8/uv_3_8.csv")

filtered_uv <- raw_uv[, grepl("_i0$", names(raw_uv)) | names(raw_uv) == "eid"]

head(filtered_uv)

column_mapping <- list(
  '1050' = 'outdoor_summer',
  '1060' = 'outdoor_winter',
  '1717' = 'skin_color',
  '1727' = 'skin_tanning',
  '1737' = 'sunburn_child',
  '1747' = 'hair_color',
  '1757' = 'facial_ageing',
  '2267' = 'uv_protection',
  '2277' = 'solarium_use'
)

names(filtered_uv) <- sapply(names(filtered_uv), function(col) {
  if (grepl("^p\\d+_i0$", col)) {
    field_id <- sub("^p(\\d+)_i0$", "\\1", col)
    if (field_id %in% names(column_mapping)) {
      return(column_mapping[[field_id]])
    }
  }
  return(col)  # eid 或未匹配的列保持不变
})

print(names(filtered_uv))

nrow(filtered_uv)

columns_to_replace <- c("outdoor_summer", "outdoor_winter")

filtered_uv[columns_to_replace] <- lapply(filtered_uv[columns_to_replace], function(x) {
  x <- as.character(x)                     # 转字符以便替换
  x[x == "-10"] <- "0.5"                   # 替换 -10 为 0.5
  x[x %in% c("-1", "-3")] <- NA            # 替换 -1 和 -3 为 NA
  return(as.numeric(x))                    # 转为数值型
})

head(filtered_uv[columns_to_replace], 10)

str(filtered_uv[columns_to_replace])

color_mapping <- c(
  "-1" = NA,
  "-3" = NA
)

columns_to_replace <- c("skin_color")

filtered_uv[columns_to_replace] <- lapply(filtered_uv[columns_to_replace], function(x) {
  x <- as.character(x)
  x <- ifelse(x %in% names(color_mapping), color_mapping[x], as.numeric(x))
  return(as.numeric(x))
})

head(filtered_uv[columns_to_replace])

columns_to_replace <- c("hair_color")

filtered_uv <- filtered_uv %>%
  mutate(across(all_of(columns_to_replace), ~ ifelse(. %in% c(-1, -3), NA, .)))

head(filtered_uv[columns_to_replace])

columns_to_replace <- c("facial_ageing")

filtered_uv <- filtered_uv %>%
  mutate(across(all_of(columns_to_replace), ~ ifelse(. %in% c(-1, -3), NA, .)))

head(filtered_uv[columns_to_replace])

columns_to_replace <- c("skin_tanning")

filtered_uv <- filtered_uv %>%
  mutate(across(all_of(columns_to_replace), ~ ifelse(. %in% c(-1, -3), NA, .)))

head(filtered_uv[columns_to_replace])

library(dplyr)

columns_to_replace <- c("uv_protection")

filtered_uv <- filtered_uv %>%
  mutate(across(all_of(columns_to_replace), ~ ifelse(. %in% c(-1, -3), NA, .)))

head(filtered_uv[columns_to_replace])

columns_to_replace <- c("sunburn_child")

filtered_uv <- filtered_uv %>%
  mutate(across(all_of(columns_to_replace), ~ ifelse(. %in% c(-1, -3), NA, .)))

head(filtered_uv[columns_to_replace])

columns_to_replace <- c("solarium_use")

filtered_uv <- filtered_uv %>%
  mutate(across(all_of(columns_to_replace), ~ ifelse(. %in% c(-1, -3), NA,
                                                     ifelse(. == -10, 0.5, .))))

head(filtered_uv[columns_to_replace])

write.csv(filtered_uv, file = "../washed_data/1_1_recoded_all_uv_data.csv", row.names = FALSE)
