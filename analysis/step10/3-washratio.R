if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())
library(tidyverse)

in_txt  <- file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step10/codes_Nightingale_biomarker_groups.txt")
out_rds <- file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step10/biomarkers_fieldid_no_to.rds")
out_no_to_csv <- file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step10/Nightingale_groups_no_to.csv")
out_with_to_csv <- file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step10/Nightingale_groups_with_to.csv")

dict <- readr::read_tsv(in_txt, show_col_types = FALSE)

names(dict) <- tolower(names(dict))

need_cols <- c("field_id", "title")
if (!all(need_cols %in% names(dict))) {
  stop("找不到必需列 field_id/title。当前列名是：\n",
       paste(names(dict), collapse = ", "))
}

dict2 <- dict %>%
  mutate(
    title_lower = tolower(title),
    has_to = str_detect(title_lower, "\\bto\\b")
  )

dict_no_to   <- dict2 %>% filter(!has_to) %>% dplyr::select(-title_lower, -has_to)
dict_with_to <- dict2 %>% filter(has_to)  %>% dplyr::select(-title_lower, -has_to)

message("Total rows: ", nrow(dict2))
message("With 'to':  ", nrow(dict_with_to))
message("No 'to':    ", nrow(dict_no_to))

fieldid_no_to <- dict_no_to %>%
  transmute(field_id = as.character(field_id)) %>%
  distinct() %>%
  pull(field_id)

saveRDS(fieldid_no_to, out_rds)

readr::write_csv(dict_no_to, out_no_to_csv)
readr::write_csv(dict_with_to, out_with_to_csv)

message("Saved RDS: ", out_rds)
message("Saved no-to CSV: ", out_no_to_csv)
message("Saved with-to CSV: ", out_with_to_csv)
