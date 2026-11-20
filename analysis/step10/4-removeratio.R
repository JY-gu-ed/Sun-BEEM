if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())
library(tidyverse)   # 你也可以只加载 dplyr/readr/stringr
library(readr)
library(stringr)
library(dplyr)

f_groups <- file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step10/Nightingale_groups_no_to.csv")
f_lm     <- file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step10/results/2_nmr_UV_lm_trend_results_FDR.csv")

out_csv  <- file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step10/results/2_nmr_UV_lm_trend_results_FDR_no_to_annotated.csv")
out_removed_csv <- file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step10/results/2_nmr_UV_lm_removed_not_in_no_to.csv")

groups_no_to <- readr::read_csv(f_groups, show_col_types = FALSE) %>%
  dplyr::mutate(field_id = as.character(field_id)) %>%
  dplyr::distinct(field_id, .keep_all = TRUE)

lmres <- readr::read_csv(f_lm, show_col_types = FALSE) %>%
  dplyr::mutate(biomarker = as.character(biomarker))

lmres2 <- lmres %>%
  dplyr::mutate(
    field_id = stringr::str_match(biomarker, "^p(\\d+)_i\\d+$")[, 2],
    field_id = as.character(field_id)
  )

bad <- lmres2 %>%
  dplyr::filter(is.na(field_id)) %>%
  dplyr::distinct(biomarker)

if (nrow(bad) > 0) {
  warning("以下 biomarker 无法解析 field_id（不是 p####_i# 格式），将被剔除：\n",
          paste(bad$biomarker, collapse = ", "))
}

keep_ids <- groups_no_to %>%
  dplyr::pull(field_id)

lm_keep <- lmres2 %>%
  dplyr::filter(!is.na(field_id)) %>%
  dplyr::filter(field_id %in% keep_ids)

lm_removed <- lmres2 %>%
  dplyr::filter(!is.na(field_id)) %>%
  dplyr::filter(!(field_id %in% keep_ids))

lm_keep_anno <- lm_keep %>%
  dplyr::left_join(
    groups_no_to %>% dplyr::select(field_id, title, group, subgroup),
    by = "field_id"
  ) %>%
  dplyr::relocate(field_id, title, group, subgroup, .after = biomarker)

if (!dir.exists(dirname(out_csv))) dir.create(dirname(out_csv), recursive = TRUE)

readr::write_csv(lm_keep_anno, out_csv)
readr::write_csv(lm_removed, out_removed_csv)

message("Total in original lm: ", nrow(lmres))
message("Kept (no-to):         ", nrow(lm_keep_anno))
message("Removed (not no-to):  ", nrow(lm_removed))
message("✅ Saved annotated:    ", out_csv)
message("✅ Saved removed list: ", out_removed_csv)
