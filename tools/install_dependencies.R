# Run from the Sun-BEEM code directory in the intended R environment.
packages <- readLines("requirements-R.txt", warn = FALSE)
packages <- trimws(packages[nzchar(trimws(packages))])
bioc <- c("clusterProfiler", "enrichplot", "org.Hs.eg.db")
available <- rownames(installed.packages())
cran_missing <- setdiff(packages, c(available, bioc, "grid"))
if (length(cran_missing)) install.packages(cran_missing, repos = "https://cloud.r-project.org")
bioc_missing <- setdiff(intersect(packages, bioc), rownames(installed.packages()))
if (length(bioc_missing)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)
}
