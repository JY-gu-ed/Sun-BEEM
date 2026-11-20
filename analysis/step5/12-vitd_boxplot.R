if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(gtsummary)
library(dplyr)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step5/"))

getwd()

beem <- read.csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "washed_data", "4_2_1_beem_vitd_sunburn_social.csv"))

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

library(ggplot2)
library(ggsignif)
library(dplyr)

if(exists("beem")) {
  beem$BEEM <- factor(beem$BEEM, levels = c("low", "medium", "high"))
  beem_clean <- beem %>% filter(!is.na(vitd) & !is.na(BEEM))
} else {
  set.seed(123)
  beem_clean <- data.frame(
    BEEM = factor(rep(c("low", "medium", "high"), each=100), levels = c("low", "medium", "high")),
    vitd = c(rnorm(100, 45, 15), rnorm(100, 55, 15), rnorm(100, 60, 15))
  )
}

main_color <- "#C77C2E"  # 边框深色
fill_color <- "#E0A663"  # 填充浅色

p <- ggplot(beem_clean, aes(x = BEEM, y = vitd)) +
  
  geom_boxplot(
    color = main_color,          # 边框颜色
    fill = fill_color,           # 填充颜色
    alpha = 0.8,                 # 透明度
    width = 0.5,                 # 箱体宽度 (正常宽度，不再是0.15)
    size = 0.7,                  # 线条粗细
    outlier.colour = main_color, # 异常点颜色 (恢复显示异常点)
    outlier.fill = main_color,
    outlier.shape = 21,          # 实心空心圆
    outlier.size = 2             # 异常点大小
  ) +
  
  theme_classic() + 
  
  labs(
    title = "Serum 25(OH)D levels by Sun-BEEM category",
    y = "Serum 25-hydroxyvitamin D (nmol/L)",
    x = NULL
  ) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(color = "black", size = 11),
    axis.title.y = element_text(size = 12),
    axis.line = element_line(color = "black")
  ) +
  
  coord_cartesian(ylim = c(0, 90))

print(p)

ggsave("SunBEEM_Boxplot_Final.pdf", width = 8, height = 5)

library(ggplot2)
library(dplyr)

if(exists("beem")) {
  beem$BEEM <- factor(beem$BEEM, levels = c("low", "medium", "high"))
  beem_clean <- beem %>% filter(!is.na(vitd) & !is.na(BEEM))
}

main_color <- "#C77C2E"
fill_color <- "#E0A663"

p <- ggplot(beem_clean, aes(x = BEEM, y = vitd)) +
  
  geom_boxplot(
    color = main_color,
    fill = fill_color,
    alpha = 0.8,
    width = 0.5,
    size = 0.7,
    outlier.colour = main_color,
    outlier.fill = main_color,
    outlier.shape = 21,
    outlier.size = 2
  ) +
  
  theme_classic() + 
  
  labs(
    y = "Serum 25(OH)D (nmol/L)",
    x = NULL  # 不需要"BEEM"这个标题，因为刻度文字已经很清楚了
  ) +
  
  scale_x_discrete(
    labels = c("Low UV", "Medium UV", "High UV") 
  ) +
  
  theme(
    axis.text.x = element_text(color = "black", size = 12, face = "bold", margin = margin(t = 5)),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title.y = element_text(size = 12),
    axis.line = element_line(color = "black"),
    
    plot.margin = margin(10, 10, 10, 10)
  ) +
  
  coord_cartesian(ylim = c(0, 150))

print(p)

library(ggplot2)
library(dplyr)

if(!exists("beem")) {
  stop("请先加载 'beem' 数据框")
}

beem_clean <- beem %>%
  filter(!is.na(vitd) & !is.na(BEEM)) %>%
  mutate(BEEM = factor(BEEM, levels = c("low", "medium", "high")))

my_palette <- c("#8FC1B5", # 清新蓝绿色 (Teal)
                "#7BA0C0", # 柔和蓝色 (Muted Blue)
                "#E4A87C") # 暖橙/赤陶色 (Terracotta)

p <- ggplot(beem_clean, aes(x = BEEM, y = vitd, fill = BEEM)) +
  
  geom_violin(
    trim = FALSE,       # FALSE 表示保留数据的长尾，展示完整分布
    alpha = 0.8,        # 设置透明度，让颜色看起来更柔和
    color = NA,         # 去掉小提琴图的外部描边，看起来更干净
    scale = "width"     # 让各组最大宽度一致，便于视觉比较
  ) +
  
  geom_boxplot(
    width = 0.15,       # 宽度要非常窄，嵌在中间
    fill = "white",     # 内部填充纯白色，突出显示
    color = "#333333",  # 边框使用深灰色而非纯黑，更显质感
    size = 0.6,         # 线条粗细
    outlier.shape = 21, # 异常点形状：空心圆
    outlier.fill = "#333333", # 异常点填充深灰色
    outlier.size = 1.8,  # 异常点大小
    alpha = 0.9
  ) +
  
  scale_fill_manual(values = my_palette) +
  
  theme_classic() + # 使用经典主题，无网格背景
  
  labs(
    y = "Serum 25(OH)D (nmol/L)", # 规范的 Y 轴单位
    x = NULL # 移除 X 轴标题
  ) +
  
  scale_x_discrete(
    labels = c("Low UV", "Medium UV", "High UV") # 将小写标签替换为更正式的写法
  ) +
  
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(size = 12, color = "black", face = "bold", margin = margin(t = 10)),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 10)),
    axis.line = element_line(color = "black", size = 0.8), # 加粗坐标轴线
    
    plot.margin = margin(20, 20, 20, 20)
  ) +
  
  
  NULL # 占位符，方便注释上一行

print(p)
