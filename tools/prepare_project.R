project_root <- Sys.getenv("SUN_BEEM_PROJECT_ROOT")
if (!nzchar(project_root)) stop("Set SUN_BEEM_PROJECT_ROOT to the project directory.")
project_root <- normalizePath(project_root, mustWork = TRUE)
folders <- c("raw_data_3_8", "washed_data", "result", "results",
             paste0("step", c(1:8, 10)),
             "step8/results", "step8/results/combine", "step10/results",
             "elife_revision/results", "JID_R1_20261104/results")
for (folder in folders) {
  dir.create(file.path(project_root, folder), recursive = TRUE, showWarnings = FALSE)
}
message("Project directories prepared.")
