if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step8/"))
getwd()

library(tidyverse)
library(ggplot2)
library(readr)

top15_cvd <- read.csv("results/vis/cvd/Top15_protective_mediators_cancer_mortality.csv")
View(top15_cvd)

top15_cancer <- read.csv("results/vis/cancer/Top15_protective_mediators_cancer_mortality.csv")
View(top15_cancer)
bm_cvd    <- top15_cvd$biomarker
bm_cancer <- top15_cancer$biomarker

overlap <- intersect(bm_cvd, bm_cancer)
cat("✅ 交集 biomarkers:\n", paste(overlap, collapse = ", "), "\n\n")

only_cvd <- setdiff(bm_cvd, bm_cancer)
only_cancer <- setdiff(bm_cancer, bm_cvd)

cat("✅ 仅 CVD Top15 有:\n", paste(only_cvd, collapse = ", "), "\n\n")
cat("✅ 仅 Cancer Top15 有:\n", paste(only_cancer, collapse = ", "), "\n")

cat("🔢 交集数量:", length(overlap), "\n")
cat("🔢 CVD 独有数量:", length(only_cvd), "\n")
cat("🔢 Cancer 独有数量:", length(only_cancer), "\n")
