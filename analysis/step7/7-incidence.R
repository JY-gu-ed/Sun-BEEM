if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(stringr)
library(purrr)
library(grid)
library(gridExtra)
library(tidyverse)
library(ggplot2)
library(stringr)
library(purrr)
library(grid)
library(gridExtra)   # ✅ 这里加上
library(cowplot)     # 拼接用

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step6/"))

df_cvd <- read.csv("../result/6_1_cvd_cox_BEEM_crude_adjusted_ptrend.csv")
df_cancerall <- read.csv("../result/7_1_1_cancerall_cox_BEEM_crude_adjusted_ptrend.csv")
df_skin_cancer <- read.csv("../result/7_1_2_skin_cancer_cox_BEEM_crude_adjusted_ptrend.csv")
df_melanoma <- read.csv("../result/7_1_3_melanoma_cox_BEEM_crude_adjusted_ptrend.csv")
df_other_skin_cancer <- read.csv("../result/7_1_4_other_skin_cancer_cox_BEEM_crude_adjusted_ptrend.csv")
df_non_skin_cancer <- read.csv("../result/7_1_5_non_skin_cancer_cox_BEEM_crude_adjusted_ptrend.csv")

library(dplyr)

final_tbl <- bind_rows(
  df_cvd              %>% mutate(outcome = "cvd_incidence"),
  df_cancerall        %>% mutate(outcome = "cancer_incidence"),
  df_skin_cancer      %>% mutate(outcome = "skin_cancer_incidence"),
  df_melanoma         %>% mutate(outcome = "melanoma_incidence"),
  df_other_skin_cancer%>% mutate(outcome = "non_melanoma_skin_cancer_incidence"),
  df_non_skin_cancer  %>% mutate(outcome = "non_skin_cancer_incidence")
)

dplyr::count(final_tbl, outcome)

outcome_labels_en <- c(
  cvd_incidence                     = "CVD incidence (narrow)",
  cancer_incidence                  = "Cancer incidence",
  skin_cancer_incidence             = "Skin cancer incidence",
  melanoma_incidence                = "Melanoma incidence",
  non_melanoma_skin_cancer_incidence= "Non-melanoma skin cancer incidence",
  non_skin_cancer_incidence         = "Non-skin cancer incidence"
)

cvd_event <- read.csv("../result/6_2_cvd_event_summary.csv")
cancerall_event <- read.csv("../result/7_1_1_cancerall_event_summary.csv")
skin_cancer_event <- read.csv("../result/7_1_2_skin_cancer_event_summary.csv")
melanoma_event <- read.csv("../result/7_1_3_melanoma_event_summary.csv")
other_skin_cancer_event <- read.csv("../result/7_1_4_other_skin_cancer_event_summary.csv")
non_skin_cancer_event <- read.csv("../result/7_1_5_non_skin_cancer_event_summary.csv")

library(dplyr)

event_tbl <- bind_rows(
  cvd_event,
  cancerall_event,
  skin_cancer_event,
  melanoma_event,
  other_skin_cancer_event,
  non_skin_cancer_event
)

library(dplyr)

n_cases_tbl <- event_tbl %>%
  select(outcome, No_cases = N_event)

print(n_cases_tbl)

parse_hrci <- function(x){
  tibble(
    estimate = as.numeric(str_extract(x, "^[0-9.]+")),
    conf.low = as.numeric(str_match(x, "\\(([0-9.]+)-")[,2]),
    conf.high= as.numeric(str_match(x, "-([0-9.]+)\\)")[,2])
  )
}

hr_long <- final_tbl %>%
  pivot_longer(
    cols = c(crude_HR_medium, crude_HR_high, adj_HR_medium, adj_HR_high),
    names_to = c("model","metric","group"),
    names_sep = "_",
    values_to = "HR_CI"
  ) %>%
  select(outcome, model, group, HR_CI)

p_long <- final_tbl %>%
  pivot_longer(
    cols = c(crude_p_medium, crude_p_high, adj_p_medium, adj_p_high),
    names_to = c("model","metric","group"),
    names_sep = "_",
    values_to = "p.value"
  ) %>%
  select(outcome, model, group, p.value)

plot_df <- hr_long %>%
  left_join(p_long, by = c("outcome","model","group")) %>%
  mutate(tmp = map(HR_CI, parse_hrci)) %>%
  unnest(tmp) %>%
  mutate(
    model = recode(model, crude = "Crude", adj = "Adjusted"),
    model = factor(model, levels = c("Crude","Adjusted")),
    BEEM_group = ifelse(group == "medium", "Medium vs Low", "High vs Low"),
    BEEM_group = factor(BEEM_group, levels = c("Medium vs Low","High vs Low"))
  )

make_bmj_forest <- function(df, n_cases_tbl, file_png, file_pdf=NULL,
                            xlim.low=0.6, xlim.up=1.2){
  outcome_order <- c(
    "non_melanoma_skin_cancer_incidence",  # ✅ 改这个
    "melanoma_incidence",
    "skin_cancer_incidence",
    "non_skin_cancer_incidence",
    "cancer_incidence",
    "cvd_incidence"
  )
  
  
  if(!"No_cases" %in% names(n_cases_tbl)){
    possible_names <- c("n_cases", "cases", "N_cases", "No.of.cases")
    found <- possible_names[possible_names %in% names(n_cases_tbl)]
    if(length(found) == 1){
      n_cases_tbl <- n_cases_tbl %>%
        dplyr::rename(No_cases = !!found)
    } else {
      stop("❌ n_cases_tbl 没有找到 No_cases 列，请检查输入文件列名。")
    }
  }
  
  df2 <- df %>%
    dplyr::filter(outcome %in% outcome_order, model=="Adjusted") %>%
    dplyr::left_join(n_cases_tbl, by="outcome") %>%
    dplyr::mutate(
      estimate  = as.numeric(estimate),
      conf.low  = as.numeric(conf.low),
      conf.high = as.numeric(conf.high),
      outcome_label = factor(outcome,
                             levels = outcome_order,
                             labels = unname(outcome_labels_en[outcome_order])),
      HR_text = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
      term = ifelse(BEEM_group=="Medium vs Low", "Medium vs Low","High vs Low")
    )
  
  p_init <- ggplot(df2, aes(y=outcome_label, x=estimate, color=term)) +
    geom_point(position=position_dodge(width=0.7))
  pb <- ggplot_build(p_init)
  df2$y_dodged <- pb$data[[1]]$y
  
  df_outcome <- df2 %>%
    dplyr::group_by(outcome_label) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(show_label = as.character(outcome_label))
  
  df_cases <- df_outcome %>%
    dplyr::select(y_dodged, No_cases)
  
  p_outcome <- ggplot(df_outcome, aes(y=y_dodged, x=1, label=show_label)) +
    geom_text(hjust=0.5, size=3.5) +
    scale_y_continuous(limits=c(0, max(df2$y_dodged)+1), expand=expansion()) +
    theme_void()
  
  p_cases <- ggplot(df_cases, aes(y=y_dodged, x=1, label=No_cases)) +
    geom_text(size=3.5) +
    scale_y_continuous(limits=c(0, max(df2$y_dodged)+1), expand=expansion()) +
    theme_void()
  
  p_forest <- ggplot(df2, aes(x=estimate, y=y_dodged, color=term)) +
    geom_vline(xintercept=1, linetype=2) +
    geom_segment(aes(x=pmax(conf.low, xlim.low),
                     xend=pmin(conf.high, xlim.up),
                     y=y_dodged, yend=y_dodged)) +
    geom_point(size=2.6) +
    geom_segment(data=dplyr::filter(df2, conf.low < xlim.low),
                 aes(x=xlim.low*1.01, xend=xlim.low,
                     y=y_dodged, yend=y_dodged, color=term),
                 arrow=arrow(angle=45, length=unit(6,"pt"))) +
    geom_segment(data=dplyr::filter(df2, conf.high > xlim.up),
                 aes(x=xlim.up*0.99, xend=xlim.up,
                     y=y_dodged, yend=y_dodged, color=term),
                 arrow=arrow(angle=45, length=unit(6,"pt"))) +
    scale_x_log10(limits=c(xlim.low, xlim.up)) +   # ⬅️ 严格截断
    scale_y_continuous(limits=c(0, max(df2$y_dodged)+1), expand=expansion()) +
    scale_color_manual(values=c("Medium vs Low"="#1f78b4","High vs Low"="#e31a1c")) +
    labs(x="Hazard Ratio", y=NULL, color=NULL) +
    theme_minimal() +
    theme(legend.position="bottom",
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  p_hr <- ggplot(df2, aes(y=y_dodged, x=1, label=HR_text)) +
    geom_text(hjust=0.5, size=3.5) +   # ⬅️ 向中间偏移
    scale_y_continuous(limits=c(0, max(df2$y_dodged)+1), expand=expansion()) +
    theme_void()
  
  g_all <- plot_grid(
    p_outcome, p_cases, p_forest, p_hr,
    ncol=4, rel_widths=c(0.25,0.15,0.35,0.25), align="h"
  )
  
  header <- cowplot::plot_grid(
    ggplot() + annotate("text", x=1, y=1, label="Outcome", size=4.2) + theme_void(),
    ggplot() + annotate("text", x=1, y=1, label="No. of cases", size=4.2) + theme_void(),
    ggplot() + annotate("text", x=1, y=1, label="Forest plot", size=4.2) + theme_void(),
    ggplot() + annotate("text", x=1, y=1, label="HR (95% CI)", size=4.2) + theme_void(),
    ncol=4, rel_widths=c(0.25,0.15,0.35,0.25)
  )
  
  g_final <- plot_grid(header, g_all, ncol=1, rel_heights=c(0.05,1))
  
  print(g_final)
  ggsave(file_png, g_final, width=12, height=8, dpi=300, bg="white")
  if(!is.null(file_pdf)) ggsave(file_pdf, g_final, width=12, height=8, bg="white")
}

n_cases_tbl <- n_cases_tbl
n_cases_tbl <- n_cases_tbl %>%
  mutate(outcome = recode(outcome,
                          cvd_incidence_narrow = "cvd_incidence",
                          all_cancer_incidence = "cancer_incidence"))

make_bmj_forest(
  plot_df,
  n_cases_tbl,
  file_png = "../result/incidence_forest_bmj_style_split.png",
  file_pdf = "../result/incidence_forest_bmj_style_split.pdf",
  xlim.low = 0.9,
  xlim.up  = 1.2
)

df_cvd <- read.csv("../result/6_1_cvd_cox_BEEM_crude_adjusted_ptrend.csv")
df_cancerall <- read.csv("../result/7_1_1_cancerall_cox_BEEM_crude_adjusted_ptrend.csv")
df_melanoma <- read.csv("../result/7_1_3_melanoma_cox_BEEM_crude_adjusted_ptrend.csv")
df_other_skin_cancer <- read.csv("../result/7_1_4_other_skin_cancer_cox_BEEM_crude_adjusted_ptrend.csv")
df_non_skin_cancer <- read.csv("../result/7_1_5_non_skin_cancer_cox_BEEM_crude_adjusted_ptrend.csv")

library(dplyr)

final_tbl <- bind_rows(
  df_cvd              %>% mutate(outcome = "cvd_incidence"),
  df_cancerall        %>% mutate(outcome = "cancer_incidence"),
  df_melanoma         %>% mutate(outcome = "melanoma_incidence"),
  df_other_skin_cancer%>% mutate(outcome = "non_melanoma_skin_cancer_incidence"),
  df_non_skin_cancer  %>% mutate(outcome = "non_skin_cancer_incidence")
)

dplyr::count(final_tbl, outcome)

outcome_labels_en <- c(
  cvd_incidence                     = "Cardiovascular disease",
  cancer_incidence                  = "Cancer",
  melanoma_incidence                = "Melanoma",
  non_melanoma_skin_cancer_incidence= "Other skin cancer",
  non_skin_cancer_incidence         = "Cancer (excluding skin cancer)"
)

cvd_event <- read.csv("../result/6_2_cvd_event_summary.csv")
cancerall_event <- read.csv("../result/7_1_1_cancerall_event_summary.csv")
skin_cancer_event <- read.csv("../result/7_1_2_skin_cancer_event_summary.csv")
melanoma_event <- read.csv("../result/7_1_3_melanoma_event_summary.csv")
other_skin_cancer_event <- read.csv("../result/7_1_4_other_skin_cancer_event_summary.csv")
non_skin_cancer_event <- read.csv("../result/7_1_5_non_skin_cancer_event_summary.csv")

library(dplyr)

event_tbl <- bind_rows(
  cvd_event,
  cancerall_event,
  skin_cancer_event,
  melanoma_event,
  other_skin_cancer_event,
  non_skin_cancer_event
)

library(dplyr)

n_cases_tbl <- event_tbl %>%
  dplyr::select(outcome, No_cases = N_event)

print(n_cases_tbl)

parse_hrci <- function(x){
  tibble(
    estimate = as.numeric(str_extract(x, "^[0-9.]+")),
    conf.low = as.numeric(str_match(x, "\\(([0-9.]+)-")[,2]),
    conf.high= as.numeric(str_match(x, "-([0-9.]+)\\)")[,2])
  )
}

hr_long <- final_tbl %>%
  pivot_longer(
    cols = c(crude_HR_medium, crude_HR_high, adj_HR_medium, adj_HR_high),
    names_to = c("model","metric","group"),
    names_sep = "_",
    values_to = "HR_CI"
  ) %>%
  dplyr::select(outcome, model, group, HR_CI)

p_long <- final_tbl %>%
  pivot_longer(
    cols = c(crude_p_medium, crude_p_high, adj_p_medium, adj_p_high),
    names_to = c("model","metric","group"),
    names_sep = "_",
    values_to = "p.value"
  ) %>%
  dplyr::select(outcome, model, group, p.value)

plot_df <- hr_long %>%
  left_join(p_long, by = c("outcome","model","group")) %>%
  mutate(tmp = map(HR_CI, parse_hrci)) %>%
  unnest(tmp) %>%
  mutate(
    model = recode(model, crude = "Crude", adj = "Adjusted"),
    model = factor(model, levels = c("Crude","Adjusted")),
    BEEM_group = ifelse(group == "medium", "Medium vs Low", "High vs Low"),
    BEEM_group = factor(BEEM_group, levels = c("Medium vs Low","High vs Low"))
  )

make_bmj_forest <- function(df, n_cases_tbl, file_png, file_pdf=NULL,
                            xlim.low=0.6, xlim.up=1.2){
  outcome_order <- c(
    "non_melanoma_skin_cancer_incidence",  # ✅ 改这个
    "melanoma_incidence",
    "non_skin_cancer_incidence",
    "cancer_incidence",
    "cvd_incidence"
  )
  
  
  if(!"No_cases" %in% names(n_cases_tbl)){
    possible_names <- c("n_cases", "cases", "N_cases", "No.of.cases")
    found <- possible_names[possible_names %in% names(n_cases_tbl)]
    if(length(found) == 1){
      n_cases_tbl <- n_cases_tbl %>%
        dplyr::rename(No_cases = !!found)
    } else {
      stop("❌ n_cases_tbl 没有找到 No_cases 列，请检查输入文件列名。")
    }
  }
  
  df2 <- df %>%
    dplyr::filter(outcome %in% outcome_order, model=="Adjusted") %>%
    dplyr::left_join(n_cases_tbl, by="outcome") %>%
    dplyr::mutate(
      estimate  = as.numeric(estimate),
      conf.low  = as.numeric(conf.low),
      conf.high = as.numeric(conf.high),
      outcome_label = factor(outcome,
                             levels = outcome_order,
                             labels = unname(outcome_labels_en[outcome_order])),
      HR_text = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
      term = ifelse(BEEM_group=="Medium vs Low", "Medium vs Low","High vs Low")
    )
  
  p_init <- ggplot(df2, aes(y=outcome_label, x=estimate, color=term)) +
    geom_point(position=position_dodge(width=0.7))
  pb <- ggplot_build(p_init)
  df2$y_dodged <- pb$data[[1]]$y
  
  df_outcome <- df2 %>%
    dplyr::group_by(outcome_label) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(show_label = as.character(outcome_label))
  
  df_cases <- df_outcome %>%
    dplyr::select(y_dodged, No_cases)
  
  p_outcome <- ggplot(df_outcome, aes(y=y_dodged, x=1, label=show_label)) +
    geom_text(hjust=0.5, size=3.5) +
    scale_y_continuous(limits=c(0, max(df2$y_dodged)+1), expand=expansion()) +
    theme_void()
  
  p_cases <- ggplot(df_cases, aes(y=y_dodged, x=1, label=No_cases)) +
    geom_text(size=3.5) +
    scale_y_continuous(limits=c(0, max(df2$y_dodged)+1), expand=expansion()) +
    theme_void()
  
  p_forest <- ggplot(df2, aes(x=estimate, y=y_dodged, color=term)) +
    geom_vline(xintercept=1, linetype=2) +
    geom_segment(aes(x=pmax(conf.low, xlim.low),
                     xend=pmin(conf.high, xlim.up),
                     y=y_dodged, yend=y_dodged)) +
    geom_point(size=2.6) +
    geom_segment(data=dplyr::filter(df2, conf.low < xlim.low),
                 aes(x=xlim.low*1.01, xend=xlim.low,
                     y=y_dodged, yend=y_dodged, color=term),
                 arrow=arrow(angle=45, length=unit(6,"pt"))) +
    geom_segment(data=dplyr::filter(df2, conf.high > xlim.up),
                 aes(x=xlim.up*0.99, xend=xlim.up,
                     y=y_dodged, yend=y_dodged, color=term),
                 arrow=arrow(angle=45, length=unit(6,"pt"))) +
    scale_x_log10(limits=c(xlim.low, xlim.up)) +   # ⬅️ 严格截断
    scale_y_continuous(limits=c(0, max(df2$y_dodged)+1), expand=expansion()) +
    scale_color_manual(values=c("Medium vs Low"="#1f78b4","High vs Low"="#e31a1c")) +
    labs(x="Hazard Ratio", y=NULL, color=NULL) +
    theme_minimal() +
    theme(legend.position="bottom",
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  p_hr <- ggplot(df2, aes(y=y_dodged, x=1, label=HR_text)) +
    geom_text(hjust=0.5, size=3.5) +   # ⬅️ 向中间偏移
    scale_y_continuous(limits=c(0, max(df2$y_dodged)+1), expand=expansion()) +
    theme_void()
  
  g_all <- plot_grid(
    p_outcome, p_cases, p_forest, p_hr,
    ncol=4, rel_widths=c(0.25,0.15,0.35,0.25), align="h"
  )
  
  header <- cowplot::plot_grid(
    ggplot() + annotate("text", x=1, y=1, label="Hospitalisations", size=4.2) + theme_void(),
    ggplot() + annotate("text", x=1, y=1, label="No. of cases", size=4.2) + theme_void(),
    ggplot() + annotate("text", x=1, y=1, label="Forest plot", size=4.2) + theme_void(),
    ggplot() + annotate("text", x=1, y=1, label="HR (95% CI)", size=4.2) + theme_void(),
    ncol=4, rel_widths=c(0.25,0.15,0.35,0.25)
  )
  
  g_final <- plot_grid(header, g_all, ncol=1, rel_heights=c(0.05,1))
  
  print(g_final)
  ggsave(file_png, g_final, width=12, height=8, dpi=300, bg="white")
  if(!is.null(file_pdf)) ggsave(file_pdf, g_final, width=12, height=8, bg="white")
}

n_cases_tbl <- n_cases_tbl
n_cases_tbl <- n_cases_tbl %>%
  mutate(outcome = recode(outcome,
                          cvd_incidence_narrow = "cvd_incidence",
                          all_cancer_incidence = "cancer_incidence"))

make_bmj_forest(
  plot_df,
  n_cases_tbl,
  file_png = "../result/del_skincancer_incidence_forest_bmj_style_split.png",
  file_pdf = "../result/del_skincancer_incidence_forest_bmj_style_split.pdf",
  xlim.low = 0.9,
  xlim.up  = 1.2
)

make_bmj_forest <- function(df, n_cases_tbl, file_png, file_pdf = NULL,
                            xlim.low = 0.6, xlim.up = 1.2) {
  outcome_order <- c(
    "non_melanoma_skin_cancer_incidence",  # ✅ incidence 图的顺序
    "melanoma_incidence",
    "non_skin_cancer_incidence",
    "cancer_incidence",
    "cvd_incidence"
  )
  
  if (!"No_cases" %in% names(n_cases_tbl)) {
    possible_names <- c("n_cases", "cases", "N_cases", "No.of.cases")
    found <- possible_names[possible_names %in% names(n_cases_tbl)]
    if (length(found) == 1) {
      n_cases_tbl <- n_cases_tbl %>%
        dplyr::rename(No_cases = !!found)
    } else {
      stop("❌ n_cases_tbl 没有找到 No_cases 列，请检查输入文件列名。")
    }
  }
  
  df2 <- df %>%
    dplyr::filter(outcome %in% outcome_order, model == "Adjusted") %>%
    dplyr::left_join(n_cases_tbl, by = "outcome") %>%
    dplyr::mutate(
      estimate      = as.numeric(estimate),
      conf.low      = as.numeric(conf.low),
      conf.high     = as.numeric(conf.high),
      outcome_label = factor(
        outcome,
        levels = outcome_order,
        labels = unname(outcome_labels_en[outcome_order])
      ),
      HR_text = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
      term    = ifelse(BEEM_group == "Medium vs Low",
                       "Medium vs Low", "High vs Low")
    )
  
  p_init <- ggplot(df2, aes(y = outcome_label, x = estimate, group = term)) +
    geom_point(position = position_dodge(width = 0.7))
  pb <- ggplot_build(p_init)
  df2$y_dodged <- pb$data[[1]]$y
  
  df_outcome <- df2 %>%
    dplyr::group_by(outcome_label) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(show_label = as.character(outcome_label))
  
  df_cases <- df_outcome %>%
    dplyr::select(y_dodged, No_cases)
  
  p_outcome <- ggplot(df_outcome, aes(y = y_dodged, x = 1, label = show_label)) +
    geom_text(hjust = 0.5, size = 4) +
    scale_y_continuous(limits = c(0, max(df2$y_dodged) + 1),
                       expand = expansion()) +
    theme_void()
  
  p_cases <- ggplot(df_cases, aes(y = y_dodged, x = 1, label = No_cases)) +
    geom_text(size = 4) +
    scale_y_continuous(limits = c(0, max(df2$y_dodged) + 1),
                       expand = expansion()) +
    theme_void()
  
  p_forest <- ggplot(df2, aes(x = estimate, y = y_dodged, shape = term)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_segment(
      aes(x = pmax(conf.low, xlim.low),
          xend = pmin(conf.high, xlim.up),
          y = y_dodged, yend = y_dodged),
      colour = "grey40"
    ) +
    geom_point(size = 2.6, colour = "grey10") +
    geom_segment(
      data = dplyr::filter(df2, conf.low < xlim.low),
      aes(x = xlim.low * 1.01, xend = xlim.low,
          y = y_dodged, yend = y_dodged),
      colour = "grey40",
      arrow  = arrow(angle = 45, length = unit(6, "pt"))
    ) +
    geom_segment(
      data = dplyr::filter(df2, conf.high > xlim.up),
      aes(x = xlim.up * 0.99, xend = xlim.up,
          y = y_dodged, yend = y_dodged),
      colour = "grey40",
      arrow  = arrow(angle = 45, length = unit(6, "pt"))
    ) +
    scale_x_log10(limits = c(xlim.low, xlim.up)) +
    scale_y_continuous(limits = c(0, max(df2$y_dodged) + 1),
                       expand = expansion()) +
    scale_shape_manual(values = c("Medium vs Low" = 16,   # 圆点
                                  "High vs Low"   = 15)) + # 方块
    labs(x = "Hazard Ratio", y = NULL, shape = NULL) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank()
    )
  
  p_hr <- ggplot(df2, aes(y = y_dodged, x = 1, label = HR_text)) +
    geom_text(hjust = 0.5, size = 4) +
    scale_y_continuous(limits = c(0, max(df2$y_dodged) + 1),
                       expand = expansion()) +
    theme_void()
  
  g_all <- plot_grid(
    p_outcome, p_cases, p_forest, p_hr,
    ncol = 4, rel_widths = c(0.25, 0.15, 0.35, 0.25), align = "h"
  )
  
  header <- cowplot::plot_grid(
    ggplot() + annotate("text", x = 1, y = 1, label = "Incidence", size = 4.2) + theme_void(),
    ggplot() + annotate("text", x = 1, y = 1, label = "No. of cases",      size = 4.2) + theme_void(),
    ggplot() + annotate("text", x = 1, y = 1, label = "Forest plot",       size = 4.2) + theme_void(),
    ggplot() + annotate("text", x = 1, y = 1, label = "HR (95% CI)",       size = 4.2) + theme_void(),
    ncol = 4, rel_widths = c(0.25, 0.15, 0.35, 0.25)
  )
  
  g_final <- plot_grid(header, g_all, ncol = 1, rel_heights = c(0.05, 1))
  
  print(g_final)
  ggsave(file_png, g_final, width = 12, height = 8, dpi = 300, bg = "white")
  if (!is.null(file_pdf)) {
    ggsave(file_pdf, g_final, width = 12, height = 8, bg = "white")
  }
}

n_cases_tbl <- n_cases_tbl %>%
  dplyr::mutate(
    outcome = dplyr::recode(
      outcome,
      cvd_incidence_narrow = "cvd_incidence",
      all_cancer_incidence = "cancer_incidence"
    )
  )

make_bmj_forest(
  plot_df,
  n_cases_tbl,
  file_png = "../result/grey_del_skincancer_incidence_forest_bmj_style_split.png",
  file_pdf = "../result/grey_del_skincancer_incidence_forest_bmj_style_split.pdf",
  xlim.low = 0.9,
  xlim.up  = 1.2
)
