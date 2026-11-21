args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript tools/run_analysis.R analysis/path/to/script.R")
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
code_root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)))
Sys.setenv(SUN_BEEM_CODE_ROOT = code_root)
relative <- gsub("\\\\", "/", args[[1]])
if (!startsWith(relative, "analysis/") || any(strsplit(relative, "/", fixed = TRUE)[[1]] == "..")) {
  stop("Choose an analysis script within this repository.")
}
script <- normalizePath(file.path(code_root, relative), mustWork = TRUE)
if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT.")
setwd(normalizePath(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), mustWork = TRUE))
source(script, local = new.env(parent = globalenv()))
