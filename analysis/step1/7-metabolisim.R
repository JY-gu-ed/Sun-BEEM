if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(tidyverse)
library(ggplot2)

setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step1/"))

getwd()

library(tidyverse)
nmr <- read.csv(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "raw_data_3_8", "NMR_field_columns.csv"))

nmr <- nmr %>%
  select(-ends_with("i1")) 

num_cols <- setdiff(names(nmr)[sapply(nmr, is.numeric)], "eid")

nmr[num_cols] <- lapply(nmr[num_cols], function(x) as.vector(scale(log1p(x))))

readr::write_csv(nmr, "../washed_data/1_7_metabolomics_nmr.csv")
