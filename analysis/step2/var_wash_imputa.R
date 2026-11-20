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

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step2/"))

getwd()

social_field_columns <- read_csv(here("raw_data_3_8", "socia_field_columns.csv"), show_col_types = FALSE) 

alcohol_columns <- read_csv(here("raw_data_3_8", "alcohol_columns.csv"), show_col_types = FALSE) 

merged_data <- inner_join(social_field_columns, alcohol_columns, by = "eid")

selected_vars <- c(
  "eid",
  "p31", "p21003_i0", "p21001_i0", "p22189", "p6138_i0",
  "p1239_i0", "p1249_i0", "p2644_i0", "p1160_i0", "p1180_i0",
  "p1200_i0", "p1210_i0", "p1220_i0", "p22040_i0",
  "p1568_i0", "p1578_i0", "p1588_i0", "p1598_i0", "p1608_i0", "p5364_i0"
)

final_data <- merged_data %>%
  select(all_of(selected_vars))

head(final_data)

proportional_imputation <- function(data, variable) {
  props <- data %>%
    filter(!is.na({{ variable }})) %>%  # Use curly-curly for tidy evaluation
    count({{ variable }}) %>%
    mutate(prop = n / sum(n))
  data %>%
    mutate({{ variable }} := ifelse(is.na({{ variable }}),
                                    sample(props[[1]], size = 1, prob = props$prop),
                                    {{ variable }})) # Use curly-curly here as well
}

mode_impute <- function(data, col) {
  if (!any(is.na(data[[col]]))) return(data)
  mode_val <- names(which.max(table(data[[col]])))
  data[[col]][is.na(data[[col]])] <- mode_val
  return(data)
}

mean_imputation <- function(data, variable) {
  mean_val <- mean(data[[variable]], na.rm = TRUE)
  data <- data %>%
    mutate(!!variable := ifelse(is.na(!!sym(variable)), mean_val, !!sym(variable)))
  return(data)
}

final_data <- mean_imputation(final_data, "p21001_i0") # bmi
final_data <- mean_imputation(final_data, "p22189") # townsend

quantile(final_data$p21003_i0, probs = seq(0, 1, 0.2), na.rm = TRUE)

extract_levels <- function(x) {
  str_extract_all(x, "[0-9-]+") %>% # Extract all numbers and hyphens
    lapply(as.integer) %>% # Convert to integers
    lapply(function(y) if (length(y) == 0) NA else y) # Handle empty strings
}

clean_social_field_columns <- final_data %>% 
  mutate(sex = case_when(p31 == 0 ~ "female",
                         p31 == 1 ~ "male")) %>%
  mutate(age = case_when(
    p21003_i0 > 36 & p21003_i0 <= 47 ~ "37-47",
    p21003_i0 > 47 & p21003_i0 <= 54 ~ "48-54",
    p21003_i0 > 54 & p21003_i0 <= 59 ~ "55-59",
    p21003_i0 > 59 & p21003_i0 <= 63 ~ "60-63",
    p21003_i0 > 63 & p21003_i0 <= 73 ~ "64-73" )) %>% 
  mutate(bmi = case_when(p21001_i0 < 18.5 ~ "underweight",
                         p21001_i0 >= 18.5 & p21001_i0 < 25 ~ "healthy weight",
                         p21001_i0 >= 25 & p21001_i0 < 30 ~ "overweight",
                         p21001_i0 >= 30 ~ "obese",
                         is.na(p21001_i0) ~ NA)) %>% 
  mutate(deprivation_quintile = case_when(
    p22189 > -6.3627 & p22189 <= -2.4167 ~ "1", # least deprived
    p22189 > -2.4167 & p22189 <= -0.4373 ~ "2",
    p22189 > -0.4373 & p22189 <= 1.7861 ~ "3",
    p22189 > 1.7861 & p22189 <= 4.7426 ~ "4",
    p22189 > 4.7426 & p22189 <= 13.5881 ~ "5" # most deprived
  )) %>%
  mutate(highest_edu = extract_levels(p6138_i0)) %>%
  mutate(highest_edu = case_when(
    sapply(highest_edu, function(x) any(x == 1)) ~ "college or university degree",
    sapply(highest_edu, function(x) any(x == 6)) ~ "other professional qualifications",
    sapply(highest_edu, function(x) any(x == 5)) ~ "NVQ, HND, HNC or equivalent",
    sapply(highest_edu, function(x) any(x == 2)) ~ "A levels or equivalent",
    sapply(highest_edu, function(x) any(x == 3)) ~ "GCSE or equivalent",
    sapply(highest_edu, function(x) any(x == 4)) ~ "CSE or equivalent",
    sapply(highest_edu, function(x) identical(x, -7L)) ~ "no qualification",
    sapply(highest_edu, function(x) identical(x, -3L)) ~ NA,
    TRUE ~ NA)) %>% 
  mutate(highest_edu = case_when(
    highest_edu == "college or university degree" |
      highest_edu == "other professional qualifications" ~ "university or above",
    highest_edu == "NVQ, HND, HNC or equivalent" |
      highest_edu == "A levels or equivalent" |
      highest_edu == "GCSE or equivalent" |
      highest_edu == "CSE or equivalent" ~ "secondary school",
    highest_edu == "no qualification" ~ "less than secondary",
    TRUE ~ NA
  )) %>%
  mutate(smoking_status = case_when(
    p1239_i0 == 1 ~ "current",
    p1249_i0 == 1 | p1249_i0 %in% c(2, 3) & p2644_i0 == 1 ~ "previous",
    p1249_i0 %in% c(2, 3) & p2644_i0 == 0 | p1249_i0 == 4 ~ "never",
    TRUE ~ NA )) %>% 
  mutate(sleep_duration = case_when(
    is.na(p1160_i0) ~ NA,
    p1160_i0 %in% c(7, 8) ~ "1",
    p1160_i0 == -3 | p1160_i0 == -1 ~ NA,
    !(p1160_i0 %in% c(7,8)) ~ "0")) %>% 
  mutate(chronotype = case_when(
    is.na(p1180_i0) ~ NA,
    p1180_i0 %in% c(1, 2) ~ "1",
    p1180_i0 == -3 | p1180_i0 == -1 ~ NA,
    !(p1180_i0 %in% c(1, 2)) ~ "0" )) %>%
  mutate(insomnia = case_when(
    is.na(p1200_i0) ~ NA,
    p1200_i0 == 1 ~ "1",
    p1200_i0 == -3 ~ NA,
    p1200_i0 != 1 ~ "0")) %>% 
  mutate(snoring = case_when(
    is.na(p1210_i0) ~ NA,
    p1210_i0 == 2 ~ "1",
    p1210_i0 == -3 ~ NA,
    p1210_i0 != 2 ~ "0")) %>% 
  mutate(day_sleepy = case_when(
    is.na(p1220_i0) ~ NA,
    p1220_i0 == 0 ~ "1",
    p1220_i0 == -3 ~ NA,
    p1220_i0 != 0 ~ "0")) %>% 
  mutate(sleep_duration = as.numeric(sleep_duration),
         chronotype = as.numeric(chronotype),
         insomnia = as.numeric(insomnia),
         snoring = as.numeric(snoring),
         day_sleepy = as.numeric(day_sleepy)) %>% 
  mutate(phys_activity = case_when(
    p22040_i0 >= 3000 ~ "high",
    p22040_i0 >= 600 & p22040_i0 < 3000 ~ "moderate",
    p22040_i0 < 600 ~ "low",
    TRUE ~ NA)) %>%
  mutate(
    red_wine = ifelse(is.na(p1568_i0) | p1568_i0 %in% c(-3, -7), 0, p1568_i0),
    white_wine = ifelse(is.na(p1578_i0) | p1578_i0 %in% c(-3, -7), 0, p1578_i0),
    beer = ifelse(is.na(p1588_i0) | p1588_i0 %in% c(-3, -7), 0, p1588_i0),
    spirits = ifelse(is.na(p1598_i0) | p1598_i0 %in% c(-3, -7), 0, p1598_i0),
    fortified_wine = ifelse(is.na(p1608_i0) | p1608_i0 %in% c(-3, -7), 0, p1608_i0),
    other_alcohol = ifelse(is.na(p5364_i0) | p5364_i0 %in% c(-3, -7), 0, p5364_i0),
    
    alcohol_perday = (
      red_wine * 2 +
        white_wine * 2 +
        beer * 2 +
        spirits +
        fortified_wine * 2 +
        other_alcohol * 2
    ) / 7,
    
    drinking = case_when(
      alcohol_perday <= 2 ~ "healthy",
      alcohol_perday > 2 ~ "unhealthy",
      TRUE ~ NA_character_
    )
  )

clean_social_field_columns <- proportional_imputation(clean_social_field_columns, sleep_duration)
clean_social_field_columns <- proportional_imputation(clean_social_field_columns, chronotype)
clean_social_field_columns <- proportional_imputation(clean_social_field_columns, insomnia)
clean_social_field_columns <- proportional_imputation(clean_social_field_columns, snoring)
clean_social_field_columns <- proportional_imputation(clean_social_field_columns, day_sleepy)

clean_social_field_columns <- clean_social_field_columns %>% 
  group_by(eid) %>% 
  mutate(sleep_score = sum(sleep_duration, chronotype, insomnia, snoring, day_sleepy)) %>% 
  ungroup() %>% 
  mutate(sleep_score = case_when(
    sleep_score >= 4 ~ "healthy",
    sleep_score >= 2 & sleep_score <= 3 ~ "intermediate",
    sleep_score <= 1 ~ "poor"
  ))

missing_summary <- data.frame(
  variable = names(clean_social_field_columns),
  missing_count = colSums(is.na(clean_social_field_columns)),
  missing_percent = colSums(is.na(clean_social_field_columns)) / nrow(clean_social_field_columns) * 100
)

missing_summary %>%
  filter(missing_count > 0) %>%
  arrange(desc(missing_percent))

write.csv(clean_social_field_columns, file = "../washed_data/2_1_1_social_all.csv", row.names = FALSE)

selected_vars <- clean_social_field_columns %>%
  select(
    eid,
    sex, age, bmi, deprivation_quintile, highest_edu,
    smoking_status, sleep_score, phys_activity, drinking
  )

missing_summary <- data.frame(
  variable = names(selected_vars),
  missing_count = colSums(is.na(selected_vars)),
  missing_percent = colSums(is.na(selected_vars)) / nrow(selected_vars) * 100
)

missing_summary %>%
  filter(missing_count > 0) %>%
  arrange(desc(missing_percent))

library(mice)

mice_data <- selected_vars %>%
  select(sex, age, bmi, deprivation_quintile, sleep_score, drinking,
         phys_activity, highest_edu, smoking_status)

mice_data <- mice_data %>% mutate(across(everything(), as.factor))

imputed <- mice(mice_data, m = 5, method = 'rf', seed = 123)

completed_data <- complete(imputed, 2)

completed_data <- cbind(eid = selected_vars$eid, completed_data)

write.csv(completed_data, "../washed_data/2_1_2_social_mice_imputed.csv", row.names = FALSE)
