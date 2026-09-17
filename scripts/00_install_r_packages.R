options(repos = c(CRAN = "https://cloud.r-project.org"))
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
options(download.file.method = "curl", download.file.extra = "--ssl-no-revoke -L")

cran_packages <- c(
  "Seurat",
  "SeuratObject",
  "data.table",
  "dplyr",
  "ggplot2",
  "patchwork",
  "remotes",
  "SoupX",
  "harmony"
)

bioc_packages <- c(
  "SingleCellExperiment",
  "scDblFinder",
  "DESeq2",
  "edgeR",
  "limma",
  "infercnv",
  "GEOquery",
  "org.Hs.eg.db",
  "celldex",
  "SingleR"
)

installed <- rownames(installed.packages())
missing_cran <- setdiff(cran_packages, installed)
if (length(missing_cran) > 0) {
  install.packages(missing_cran, dependencies = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

installed <- rownames(installed.packages())
missing_bioc <- setdiff(bioc_packages, installed)
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message("Package installation step finished.")
