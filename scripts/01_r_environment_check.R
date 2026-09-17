dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

required_cran <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "data.table",
  "dplyr",
  "ggplot2",
  "patchwork",
  "remotes",
  "SoupX",
  "harmony"
)

required_bioc <- c(
  "SingleCellExperiment",
  "scDblFinder",
  "DESeq2",
  "edgeR",
  "limma",
  "infercnv",
  "GEOquery",
  "celldex",
  "SingleR"
)

all_required <- c(required_cran, required_bioc)
installed <- rownames(installed.packages())
status <- data.frame(
  package = all_required,
  source = c(rep("CRAN/GitHub", length(required_cran)), rep("Bioconductor", length(required_bioc))),
  installed = all_required %in% installed,
  stringsAsFactors = FALSE
)

print(status, row.names = FALSE)

missing <- status$package[!status$installed]
if (length(missing) > 0) {
  message("\nMissing packages:")
  message(paste(missing, collapse = ", "))
  message("\nInstall only after confirming network access and disk/time budget.")
} else {
  message("\nAll required packages are installed.")
}

sessionInfo()
