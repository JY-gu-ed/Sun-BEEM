if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the project directory.")
setwd(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step10"))
dir.create(file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "step10/supplement"),
           recursive = TRUE,
           showWarnings = FALSE)

library(dplyr)
n_boot <- as.numeric(Sys.getenv("SUN_BEEM_N_BOOT", unset = "1000"))
seed <- as.numeric(Sys.getenv("SUN_BEEM_SEED", unset = "123"))
if (!is.finite(n_boot) || n_boot < 10 || n_boot != floor(n_boot) || n_boot > .Machine$integer.max ||
    !is.finite(seed) || seed < 0 || seed != floor(seed) || seed > .Machine$integer.max) stop("Invalid bootstrap settings.")
input_dir <- file.path("results/combine", paste0("bootstrap_", as.integer(n_boot), "_seed_", as.integer(seed)))
output_dir <- file.path(input_dir, "supplement")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cvd <- read.csv(
  file.path(input_dir, "cvd_mortality_mt_uv_mediation_significant_hits_annotated__keepTitlesIn_Nightingale_no_to__KEPT.csv")
)

head(cvd)
dim(cvd)
cancer <- read.csv(
  file.path(input_dir, "cancer_mt_uv_mediation_significant_hits_annotated__keepTitlesIn_Nightingale_no_to__KEPT.csv")
)

head(cancer)
dim(cancer)
cvd_protective_axb <- cvd %>%
  filter(
    !is.na(mediation_effect),
    !is.na(mediation_ci_low),
    !is.na(mediation_ci_high),
    mediation_effect < 0,
    mediation_ci_high < 0
  )

message("Metabolomic protective mediators (CVD): ",
        nrow(cvd_protective_axb))


cancer_protective_axb <- cancer %>%
  filter(
    !is.na(mediation_effect),
    !is.na(mediation_ci_low),
    !is.na(mediation_ci_high),
    mediation_effect < 0,
    mediation_ci_high < 0
  )

message("Metabolomic protective mediators (Cancer): ",
        nrow(cancer_protective_axb))

write.csv(
  cvd_protective_axb,
  file.path(output_dir, "Supplementary_Data_S3_Metabolomic_Protective_Mediators_CVD.csv"),
  row.names = FALSE
)

write.csv(
  cancer_protective_axb,
  file.path(output_dir, "Supplementary_Data_S4_Metabolomic_Protective_Mediators_Cancer.csv"),
  row.names = FALSE
)
