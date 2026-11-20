if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(here)
library(ggplot2)
library(tidyverse)
library(here)
library(ggplot2)
library(dplyr)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step1/"))

getwd()

merged_uv <- read.csv("../washed_data/1_2_merged_uv_data.csv")
table(merged_uv$solarium_use)
anyNA(merged_uv)
library(dplyr)

merged_uv <- merged_uv %>%
  mutate(outdoor_avg = rowMeans(across(c(outdoor_summer, outdoor_winter)), na.rm = TRUE))

library(ggplot2)

if (!dir.exists("pic")) dir.create("pic")

plot_and_save <- function(data, varname, binwidth = NULL) {
  p <- ggplot(data, aes_string(x = varname)) +
    geom_bar(fill = "steelblue") +
    theme_minimal() +
    labs(title = paste("Distribution of", varname), x = varname, y = "Count")
  
  if (!is.null(binwidth)) {
    p <- ggplot(data, aes_string(x = varname)) +
      geom_histogram(binwidth = binwidth, fill = "steelblue", color = "white") +
      theme_minimal() +
      labs(title = paste("Histogram of", varname), x = varname, y = "Frequency")
  }
  
  ggsave(filename = paste0("pic/", varname, ".png"), plot = p, width = 6, height = 4)
}

plot_and_save(merged_uv, "outdoor_summer")
plot_and_save(merged_uv, "outdoor_winter")
plot_and_save(merged_uv, "uv_protection")
plot_and_save(merged_uv, "solarium_use")
plot_and_save(merged_uv, "swr2008_kJ")
plot_and_save(merged_uv,"outdoor_avg")

merged_uv <- merged_uv %>%
  mutate(bi_outdoor = ifelse(outdoor_avg >= 4, 1, 0))

merged_uv %>% 
  ggplot(aes(x = uv_protection)) +
  geom_histogram(bins = 5) # u-shaped

merged_uv %>%
  group_by(uv_protection) %>% 
  summarise(mean_outdoor = mean(outdoor_avg, na.rm = TRUE)) %>%
  ggplot(aes(x = uv_protection, y = mean_outdoor)) +
  geom_col() # u-shaped

merged_uv %>%
  group_by(uv_protection) %>% 
  summarise(mean_swr = mean(swr2008_kJ, na.rm = TRUE)) %>%
  ggplot(aes(x = uv_protection, y = mean_swr)) +
  geom_col() # u-shaped

merged_uv %>%
  group_by(solarium_use, uv_protection) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = uv_protection, y = n, fill = factor(solarium_use))) +
  geom_col(position = "dodge") # u-shaped

table(merged_uv$uv_protection)
merged_uv <- merged_uv %>%
  mutate(bi_uv_protection = ifelse(uv_protection %in% c(1, 5), 0, 1))

summary(merged_uv$swr2008_kJ)
quantile(merged_uv$swr2008_kJ, 
         probs = c(1/3, 2/3), 
         na.rm = TRUE)
merged_uv <- merged_uv %>%
  mutate(bi_swr2008 = ifelse(swr2008_kJ < 9516, 0, 1))

table(merged_uv$solarium_use)
merged_uv <- merged_uv %>%
  mutate(bi_solarium_use = ifelse(solarium_use %in% c(0, 0.5), 0, 1))

table(merged_uv$solarium_use, merged_uv$bi_solarium_use)

library(ggplot2)

plot_and_save_binary <- function(data, varname, folder = "pic") {
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  
  p <- ggplot(data, aes_string(x = varname)) +
    geom_bar(fill = "steelblue") +
    labs(title = paste("Distribution of", varname),
         x = varname,
         y = "Count") +
    theme_minimal()
  
  ggsave(filename = paste0(folder, "/", varname, "_barplot.png"),
         plot = p, width = 5, height = 4, dpi = 300)
}

binary_vars <- c("bi_outdoor", "bi_uv_protection", "bi_swr2008", "bi_solarium_use")

for (var in binary_vars) {
  plot_and_save_binary(merged_uv, var)
}

merged_uv <- merged_uv %>%
  mutate(bi_uv_score = bi_outdoor + bi_uv_protection + bi_swr2008 + bi_solarium_use)
table(merged_uv$bi_uv_score, useNA = "ifany")

p <- ggplot(merged_uv, aes(x = bi_uv_score)) +
  geom_bar(fill = "darkgreen") +
  labs(title = "Distribution of Combined UV Risk Score",
       x = "UV Risk Score (0–4)", y = "Count") +
  theme_minimal()

if (!dir.exists("pic")) dir.create("pic", recursive = TRUE)

ggsave(filename = "pic/bi_uv_score_barplot.png", plot = p,
       width = 5, height = 4, dpi = 300)

library(dplyr)

merged_uv <- merged_uv %>%
  mutate(
    BEEM = case_when(
      bi_uv_score %in% c(0, 1) ~ "low",
      bi_uv_score == 2         ~ "medium",
      bi_uv_score %in% c(3, 4) ~ "high"
    )
  )

head(merged_uv)
table(merged_uv$BEEM)

p <- ggplot(merged_uv, aes(x = BEEM)) +
  geom_bar(fill = "tomato") +
  labs(title = "Distribution of BEEM Risk Categories",
       x = "BEEM Category", y = "Count") +
  theme_minimal()

if (!dir.exists("pic")) dir.create("pic", recursive = TRUE)

ggsave(filename = "pic/BEEM_barplot.png", plot = p,
       width = 5, height = 4, dpi = 300)

old_uv <- read.csv(file.path(Sys.getenv("SUN_BEEM_EXTERNAL_ROOT"), "uv_variable.csv"))

library(dplyr)

old_group <- old_uv %>%
  select(eid, uv_category)

new_group <- merged_uv %>%
  select(eid, BEEM)

compare_df <- inner_join(old_group, new_group, by = "eid")

compare_df <- compare_df %>%
  mutate(match = ifelse(uv_category == BEEM, "match", "mismatch"))

table(compare_df$match)

mismatch_cases <- compare_df %>%
  filter(match == "mismatch")

head(mismatch_cases)

write.csv(merged_uv, file = "../washed_data/1_3_bi_merged_uv.csv", row.names = FALSE)
