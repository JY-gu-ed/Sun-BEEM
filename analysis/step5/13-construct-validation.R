if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT")) || !nzchar(Sys.getenv("SUN_BEEM_CODE_ROOT"))) {
  stop("Set the project/code directories or use tools/run_analysis.R.")
}
suppressPackageStartupMessages({library(data.table); library(broom)})
source(file.path(Sys.getenv("SUN_BEEM_CODE_ROOT"), "analysis/shared/construct_validation.R"))
root <- Sys.getenv("SUN_BEEM_PROJECT_ROOT")
d <- as.data.frame(fread(file.path(root, "washed_data/4_2_2_ana_for_mortality.csv")))
metadata_path <- Sys.getenv("SUN_BEEM_BASELINE_METADATA")
candidates <- if (nzchar(metadata_path)) metadata_path else file.path(root, "raw_data_3_8", c("request_an_columns.csv", "socia_field_columns.csv"))
metadata <- NULL
for (path in candidates) {
  if (file.exists(path) && all(c("eid", "p54_i0") %in% names(fread(path, nrows = 0)))) {
    metadata <- fread(path, select = c("eid", "p54_i0"))
    break
  }
}
if (is.null(metadata)) stop("Provide baseline assessment centre field p54_i0 via SUN_BEEM_BASELINE_METADATA.")
stopifnot(!anyDuplicated(metadata$eid), !anyDuplicated(d$eid))
d$assessment_centre <- metadata$p54_i0[match(d$eid, metadata$eid)]
d$season <- baseline_season(d$date_assess)
out <- file.path(root, "result/construct_validation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
results <- list()
for (response in c("vitd", "sunburn_child")) {
  fit <- fit_construct_comparison(d, response)
  row <- as.data.table(broom::tidy(fit, conf.int = TRUE))
  row <- row[grepl("^BEEM", term)]
  row[, `:=`(response = response, N = nobs(fit),
             adjustment = "age_group + sex + BMI_category + season + assessment_centre")]
  results[[response]] <- row
}
fwrite(rbindlist(results), file.path(out, "category_adjusted_comparisons.csv"))
components <- c("bi_outdoor", "bi_solarium_use", "bi_uv_protection", "bi_swr2008")
contrasts <- list()
for (response in c("vitd", "sunburn_child")) for (component in components) {
  z <- d[complete.cases(d[, c(response, component)]), c(response, component), drop = FALSE]
  fit <- lm(reformulate(component, response = response), data = z)
  row <- as.data.table(broom::tidy(fit, conf.int = TRUE))[term == component]
  row[, `:=`(response = response, component = component, N = nobs(fit))]
  contrasts[[paste(response, component)]] <- row
}
fwrite(rbindlist(contrasts), file.path(out, "component_comparisons.csv"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
