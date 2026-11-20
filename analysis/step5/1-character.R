if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(here)
library(ggplot2)
library(tidyverse)
library(here)
library(ggplot2)
library(gtsummary)
library(dplyr)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step5/"))

getwd()

beem <- read.csv(here("washed_data", "4_2_1_beem_vitd_sunburn_social.csv"))

beem_backup <- beem 

beem <- beem %>%
  mutate(
    deprivation_quintile = as.numeric(as.character(deprivation_quintile)),
    BEEM = factor(BEEM, levels = c("low", "medium", "high")),
    bmi = factor(bmi, levels = c("underweight", "healthy weight", "overweight", "obese")),
    smoking_status = factor(smoking_status, levels = c("current", "previous", "never")),
    phys_activity = factor(phys_activity, levels = c("low", "moderate", "high")),
    sleep_score = factor(sleep_score, levels = c("poor", "intermediate", "healthy")),
    drinking = factor(drinking, levels = c("unhealthy", "healthy")),
    across(where(is.character), as.factor)
  )

continuous_vars <- c("vitd", "sunburn_child")
categorical_vars <- c("sex", "age", "bmi", "highest_edu",
                      "smoking_status", "phys_activity",
                      "sleep_score", "drinking","deprivation_quintile")

var_order <- c("sex", "age", "bmi", "deprivation_quintile",
               "highest_edu", "phys_activity",
               "sleep_score","smoking_status",  "drinking", "sunburn_child", "vitd")

tbl <- beem %>%
  select(BEEM, all_of(var_order)) %>%
  tbl_summary(
    by = BEEM,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    missing = "no"
  ) %>%
  modify_table_body(~ .x %>% filter(!(variable == "sex" & label != "male")))

p_values <- sapply(var_order, function(var) {
  if (var %in% categorical_vars) {
    chisq.test(table(beem[[var]], beem$BEEM))$p.value
  } else {
    summary(aov(beem[[var]] ~ beem$BEEM))[[1]]["Pr(>F)"][1, 1]
  }
})

tbl <- tbl %>%
  modify_table_body(
    ~ .x %>%
      left_join(
        tibble(
          variable = var_order,
          p.value = formatC(p_values, format = "f", digits = 4)
        ),
        by = "variable"
      )
  ) %>%
  modify_header(p.value ~ "**p-value**") %>%
  bold_labels()

tbl

tbl %>%
  as_flex_table() %>% 
  flextable::save_as_docx(path = "../result/5_1_baseline_table.docx")

tbl_df <- tbl %>%
  as_tibble(col_names = TRUE)  # 保留列名

write.csv(tbl_df, "../result/5_1_baseline_table.csv", row.names = FALSE)

beem <- beem_backup %>%
  mutate(
    deprivation_quintile = as.numeric(as.character(deprivation_quintile)),
    bi_solarium_use = factor(bi_solarium_use, levels = c(0, 1),
                             labels = c("No", "Yes")),
    bmi = factor(bmi, levels = c("underweight", "healthy weight", "overweight", "obese")),
    smoking_status = factor(smoking_status, levels = c("current", "previous", "never")),
    phys_activity = factor(phys_activity, levels = c("low", "moderate", "high")),
    sleep_score = factor(sleep_score, levels = c("poor", "intermediate", "healthy")),
    drinking = factor(drinking, levels = c("unhealthy", "healthy")),
    across(where(is.character), as.factor)
  )

continuous_vars <- c("vitd", "sunburn_child")
categorical_vars <- c("sex", "age", "bmi", "highest_edu",
                      "smoking_status", "phys_activity",
                      "sleep_score", "drinking", "deprivation_quintile")

var_order <- c("sex", "age", "bmi", "deprivation_quintile",
               "highest_edu", "phys_activity", "sleep_score",
               "smoking_status", "drinking", "sunburn_child", "vitd")

tbl <- beem %>%
  select(bi_solarium_use, all_of(var_order)) %>%
  tbl_summary(
    by = bi_solarium_use,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    missing = "no"
  ) %>%
  modify_table_body(~ .x %>% filter(!(variable == "sex" & label != "male")))

p_values <- sapply(var_order, function(var) {
  if (var %in% categorical_vars) {
    chisq.test(table(beem[[var]], beem$bi_solarium_use))$p.value
  } else {
    summary(aov(beem[[var]] ~ beem$bi_solarium_use))[[1]]["Pr(>F)"][1, 1]
  }
})

tbl <- tbl %>%
  modify_table_body(
    ~ .x %>%
      left_join(
        tibble(
          variable = var_order,
          p.value = formatC(p_values, format = "f", digits = 4)
        ),
        by = "variable"
      )
  ) %>%
  modify_header(p.value ~ "**p-value**") %>%
  bold_labels()

tbl

tbl %>%
  as_flex_table() %>% 
  flextable::save_as_docx(path = "../result/5_2_bi_solarium_baseline_table.docx")

tbl_df <- tbl %>%
  as_tibble(col_names = TRUE)  # 保留列名

write.csv(tbl_df, "../result/5_2_bi_solarium_baseline_table.csv", row.names = FALSE)

beem <- beem_backup %>%
  mutate(
    deprivation_quintile = as.numeric(as.character(deprivation_quintile)),
    bi_swr2008 = factor(bi_swr2008, levels = c(0, 1),
                        labels = c("No", "Yes")),
    bmi = factor(bmi, levels = c("underweight", "healthy weight", "overweight", "obese")),
    smoking_status = factor(smoking_status, levels = c("current", "previous", "never")),
    phys_activity = factor(phys_activity, levels = c("low", "moderate", "high")),
    sleep_score = factor(sleep_score, levels = c("poor", "intermediate", "healthy")),
    drinking = factor(drinking, levels = c("unhealthy", "healthy")),
    across(where(is.character), as.factor)
  )

continuous_vars <- c("vitd", "sunburn_child")
categorical_vars <- c("sex", "age", "bmi", "highest_edu",
                      "smoking_status", "phys_activity",
                      "sleep_score", "drinking", "deprivation_quintile")

var_order <- c("sex", "age", "bmi", "deprivation_quintile",
               "highest_edu", "phys_activity", "sleep_score",
               "smoking_status", "drinking", "sunburn_child", "vitd")

tbl <- beem %>%
  select(bi_swr2008, all_of(var_order)) %>%
  tbl_summary(
    by = bi_swr2008,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    missing = "no"
  ) %>%
  modify_table_body(~ .x %>% filter(!(variable == "sex" & label != "male")))

p_values <- sapply(var_order, function(var) {
  if (var %in% categorical_vars) {
    chisq.test(table(beem[[var]], beem$bi_swr2008))$p.value
  } else {
    summary(aov(beem[[var]] ~ beem$bi_swr2008))[[1]]["Pr(>F)"][1, 1]
  }
})

tbl <- tbl %>%
  modify_table_body(
    ~ .x %>%
      left_join(
        tibble(
          variable = var_order,
          p.value = formatC(p_values, format = "f", digits = 4)
        ),
        by = "variable"
      )
  ) %>%
  modify_header(p.value ~ "**p-value**") %>%
  bold_labels()

tbl %>%
  as_flex_table() %>% 
  flextable::save_as_docx(path = "../result/5_3_bi_swr2008_baseline_table.docx")

tbl

tbl_df <- tbl %>%
  as_tibble(col_names = TRUE)

write.csv(tbl_df, "../result/5_3_bi_swr2008_baseline_table.csv", row.names = FALSE)

beem <- beem_backup %>%
  mutate(
    deprivation_quintile = as.numeric(as.character(deprivation_quintile)),
    bi_uv_protection = factor(bi_uv_protection, levels = c(0, 1),
                              labels = c("No", "Yes")),
    bmi = factor(bmi, levels = c("underweight", "healthy weight", "overweight", "obese")),
    smoking_status = factor(smoking_status, levels = c("current", "previous", "never")),
    phys_activity = factor(phys_activity, levels = c("low", "moderate", "high")),
    sleep_score = factor(sleep_score, levels = c("poor", "intermediate", "healthy")),
    drinking = factor(drinking, levels = c("unhealthy", "healthy")),
    across(where(is.character), as.factor)
  )

continuous_vars <- c("vitd", "sunburn_child")
categorical_vars <- c("sex", "age", "bmi", "highest_edu",
                      "smoking_status", "phys_activity",
                      "sleep_score", "drinking", "deprivation_quintile")

var_order <- c("sex", "age", "bmi", "deprivation_quintile",
               "highest_edu", "phys_activity", "sleep_score",
               "smoking_status", "drinking", "sunburn_child", "vitd")

tbl <- beem %>%
  select(bi_uv_protection, all_of(var_order)) %>%
  tbl_summary(
    by = bi_uv_protection,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    missing = "no"
  ) %>%
  modify_table_body(~ .x %>% filter(!(variable == "sex" & label != "male")))

p_values <- sapply(var_order, function(var) {
  if (var %in% categorical_vars) {
    chisq.test(table(beem[[var]], beem$bi_uv_protection))$p.value
  } else {
    summary(aov(beem[[var]] ~ beem$bi_uv_protection))[[1]]["Pr(>F)"][1, 1]
  }
})

tbl <- tbl %>%
  modify_table_body(
    ~ .x %>%
      left_join(
        tibble(
          variable = var_order,
          p.value = formatC(p_values, format = "f", digits = 4)
        ),
        by = "variable"
      )
  ) %>%
  modify_header(p.value ~ "**p-value**") %>%
  bold_labels()

tbl
tbl %>%
  as_flex_table() %>% 
  flextable::save_as_docx(path = "../result/5_4_bi_uv_protection_baseline_table.docx")

tbl_df <- tbl %>%
  as_tibble(col_names = TRUE)

write.csv(tbl_df, "../result/5_4_bi_uv_protection_baseline_table.csv", row.names = FALSE)

beem <- beem_backup %>%
  mutate(
    deprivation_quintile = as.numeric(as.character(deprivation_quintile)),
    bi_outdoor = factor(bi_outdoor, levels = c(0, 1),
                        labels = c("No", "Yes")),
    bmi = factor(bmi, levels = c("underweight", "healthy weight", "overweight", "obese")),
    smoking_status = factor(smoking_status, levels = c("current", "previous", "never")),
    phys_activity = factor(phys_activity, levels = c("low", "moderate", "high")),
    sleep_score = factor(sleep_score, levels = c("poor", "intermediate", "healthy")),
    drinking = factor(drinking, levels = c("unhealthy", "healthy")),
    across(where(is.character), as.factor)
  )

continuous_vars <- c("vitd", "sunburn_child")
categorical_vars <- c("sex", "age", "bmi", "highest_edu",
                      "smoking_status", "phys_activity",
                      "sleep_score", "drinking", "deprivation_quintile")

var_order <- c("sex", "age", "bmi", "deprivation_quintile",
               "highest_edu", "phys_activity", "sleep_score",
               "smoking_status", "drinking", "sunburn_child", "vitd")

tbl <- beem %>%
  select(bi_outdoor, all_of(var_order)) %>%
  tbl_summary(
    by = bi_outdoor,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    missing = "no"
  ) %>%
  modify_table_body(~ .x %>% filter(!(variable == "sex" & label != "male")))

p_values <- sapply(var_order, function(var) {
  if (var %in% categorical_vars) {
    chisq.test(table(beem[[var]], beem$bi_outdoor))$p.value
  } else {
    summary(aov(beem[[var]] ~ beem$bi_outdoor))[[1]]["Pr(>F)"][1, 1]
  }
})

tbl <- tbl %>%
  modify_table_body(
    ~ .x %>%
      left_join(
        tibble(
          variable = var_order,
          p.value = formatC(p_values, format = "f", digits = 4)
        ),
        by = "variable"
      )
  ) %>%
  modify_header(p.value ~ "**p-value**") %>%
  bold_labels()

tbl %>%
  as_flex_table() %>% 
  flextable::save_as_docx(path = "../result/5_5_bi_outdoor_baseline_table.docx")

tbl_df <- tbl %>%
  as_tibble(col_names = TRUE)

write.csv(tbl_df, "../result/5_5_bi_outdoor_baseline_table.csv", row.names = FALSE)
tbl

library(patchwork) 

dat <- beem_backup
num_vars <- c("solarium_use", "swr2008_kJ", "outdoor_avg")
sunburn_ymax <- quantile(dat$sunburn_child, 0.99, na.rm = TRUE)

base_theme <- theme_minimal(base_size = 12)

p_vitd_uv <- ggplot(dat, aes(x = factor(uv_protection), y = vitd, fill = factor(uv_protection))) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  labs(x = "UV protection", y = "Vitamin D") +
  base_theme + theme(legend.position = "none")

p_burn_uv <- ggplot(dat, aes(x = factor(uv_protection), y = sunburn_child, fill = factor(uv_protection))) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  coord_cartesian(ylim = c(0, sunburn_ymax)) +
  labs(x = "UV protection", y = "Childhood sunburn count") +
  base_theme + theme(legend.position = "none")

mk_scatter <- function(xvar, yvar, ylab){
  ggplot(dat, aes_string(x = xvar, y = yvar)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(x = xvar, y = ylab) +
    base_theme
}

p_vitd_sol  <- mk_scatter("solarium_use", "vitd", "Vitamin D")
p_vitd_swr  <- mk_scatter("swr2008_kJ",  "vitd", "Vitamin D")
p_vitd_out  <- mk_scatter("outdoor_avg", "vitd", "Vitamin D")

p_burn_sol  <- mk_scatter("solarium_use", "sunburn_child", "Childhood sunburn count") + 
  coord_cartesian(ylim = c(0, sunburn_ymax))
p_burn_swr  <- mk_scatter("swr2008_kJ",  "sunburn_child", "Childhood sunburn count") + 
  coord_cartesian(ylim = c(0, sunburn_ymax))
p_burn_out  <- mk_scatter("outdoor_avg", "sunburn_child", "Childhood sunburn count") + 
  coord_cartesian(ylim = c(0, sunburn_ymax))

fig <- (p_vitd_uv | p_burn_uv) /
  (p_vitd_sol | p_burn_sol) /
  (p_vitd_swr | p_burn_swr) /
  (p_vitd_out | p_burn_out) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig

ggsave("../result/6_fig_uv_solarium_swr_outdoor_vs_vitd_sunburn.png",
       fig, width = 14, height = 16, dpi = 300)
ggsave("../result/6_fig_uv_solarium_swr_outdoor_vs_vitd_sunburn.pdf",
       fig, width = 14, height = 16)
