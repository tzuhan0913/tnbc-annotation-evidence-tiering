# =============================================================================
# Tier 3 re-run: pathfindR active-subnetwork enrichment on the harmonised labels.
#
# Uses exactly the parameters of the original analysis (Functions.R / run_pathfindR_dtls):
#   pin_name_path = "STRING", min_gset_size = 15, max_gset_size = 1000,
#   p_val_threshold = 0.05, enrichment_threshold = 0.05, iterations = 50,
#   convert2alias = FALSE, databases GO-CC / GO-BP / GO-MF / KEGG.
#
# What changes is the INPUT gene list:
#   "patient"  (default, recommended) - the Tier 2 patient-level genes
#                (edgeR QL paired model, FDR < 0.05 and |logFC| > 0.5; ~830 genes).
#                This keeps Tier 3 from inheriting a cell-level gene list at all.
#   "cell"     - the Tier 1 cell-level genes under the harmonised labels, capped at
#                TOP_N_CELL by adjusted P value so the run stays tractable.
#   "legacy"   - the original 869-gene list, for a like-for-like comparison.
#
# Requirements: pathfindR (and a working Java, which pathfindR's active-subnetwork
# search needs), dplyr. Expect roughly 10-40 min per database.
#
#   setwd("C:/Users/tzuha/Documents/投稿文章(JS)")
#   source("scripts/run_tier3_pathfindr.R")
# =============================================================================
INPUT_SETS  <- c("patient", "cell")     # add "legacy" to also rerun the original list
TOP_N_CELL  <- 1000                     # cap for the cell-level list
DBS         <- c("KEGG", "GO-BP", "GO-CC", "GO-MF")
ITER        <- 50
N_PROC      <- 4                        # raising this too far makes pathfindR flaky
# pathfindR's active-subnetwork search shells out to Java with the output path on the command
# line.  The project folder is "投稿文章(JS)" - non-ASCII plus parentheses - which is exactly the
# kind of path that breaks such calls on Windows (the original analysis used an absolute
# D:/Downloads path for the same reason).  So pathfindR runs in a plain-ASCII working folder
# and only the result tables are written back into the project.
WORK_DIR <- file.path(if (nzchar(Sys.getenv("LOCALAPPDATA"))) Sys.getenv("LOCALAPPDATA") else tempdir(),
                      "pathfindr_run")

suppressPackageStartupMessages({ library(data.table); library(dplyr) })
find_root <- function(marker) {
  ## marker must be a folder the scripts never create themselves: an earlier version keyed on
  ## "analysis/...", which the buggy run had already created under scripts/, so it kept matching there.
  p <- getwd()
  for (i in 0:4) {
    if (dir.exists(file.path(p, marker))) return(normalizePath(p, winslash = "/", mustWork = FALSE))
    p <- file.path(p, "..")
  }
  stop("Cannot find '", marker, "' from ", getwd(),
       '\n  setwd("C:/Users/tzuha/Documents/投稿文章(JS)") and source the script again.')
}
root <- find_root(file.path("data", "GSE161529"))
cat("project root:", root, "\n")
indir  <- file.path(root, "analysis", "epithelial_only_pseudobulk_R")
outdir <- file.path(root, "analysis", "tier3_pathfindr_reanalysis")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(WORK_DIR, showWarnings = FALSE, recursive = TRUE)
cat("pathfindR working folder:", WORK_DIR, "\n")
if (grepl("[^ -~]", WORK_DIR) || grepl("[()]", WORK_DIR))
  warning("WORK_DIR still contains non-ASCII characters or parentheses; ",
          "pathfindR's Java step may fail. Set WORK_DIR to something like C:/pathfindr_run.")

if (!requireNamespace("pathfindR", quietly = TRUE))
  stop('pathfindR is not installed.  install.packages("pathfindR")  (it also needs Java).')
library(pathfindR)
cat("pathfindR version:", as.character(packageVersion("pathfindR")), "\n")
## The active-subnetwork search is a bundled Java archive.  Antivirus heuristics (Avast's
## IDP.HELU.JarMal* in particular) quarantine unsigned JARs, which makes runs fail - or worse,
## silently return fewer subnetworks.  Fail loudly here rather than hours later.
jdir <- system.file("java", package = "pathfindR")
jars <- if (nzchar(jdir)) list.files(jdir, pattern = "\\.jar$") else character(0)
cat("pathfindR Java payload:", if (length(jars)) paste(jars, collapse = ", ") else "MISSING", "\n")
if (!length(jars))
  stop("pathfindR's ActiveSubnetworkSearch.jar is missing from\n  ", jdir,
       "\nIt has most likely been quarantined by antivirus software. Restore it, add the\n",
       "pathfindR library folder and ", WORK_DIR, " to the antivirus exclusion list,\n",
       'then reinstall with install.packages("pathfindR", type = "binary") and run again.')
if (system("java -version") != 0) stop("java is not callable from R; check JAVA_HOME and PATH.")

# ---- build the input gene tables -------------------------------------------
make_input <- function(which_set) {
  if (which_set == "patient") {
    f <- file.path(indir, "tier2a_paired_edgeR_QL_min10.csv")
    if (!file.exists(f)) stop("missing ", f)
    d <- fread(f)
    d <- d[FDR < 0.05 & abs(logFC) > 0.5]
    data.frame(Gene.symbol = d$gene, logFC = d$logFC, p_val = d$FDR)
  } else if (which_set == "cell") {
    f <- file.path(indir, "tier1_cell_level_findmarkers.csv")
    if (!file.exists(f)) stop("missing ", f)
    d <- fread(f)
    d <- d[p_val_adj < 0.05 & abs(avg_log2FC) > 0.5][order(p_val)]
    if (nrow(d) > TOP_N_CELL) d <- d[seq_len(TOP_N_CELL)]
    # pathfindR rejects p = 0, which Seurat returns after Bonferroni underflow
    p <- pmax(d$p_val_adj, .Machine$double.xmin)
    data.frame(Gene.symbol = d$gene, logFC = d$avg_log2FC, p_val = p)
  } else {
    f <- file.path(root, "inputs", "liao_master_zip_20260707",
                   "05-家嫺論文補分析", "Data", "data2.csv")
    if (!file.exists(f)) stop("missing the original 869-gene list at ", f)
    d <- fread(f)
    setnames(d, trimws(names(d)))
    data.frame(Gene.symbol = trimws(d$Gene), logFC = d$avg_log2FC,
               p_val = pmax(d$p_val_adj, .Machine$double.xmin))
  }
}

all_terms <- list()
failed <- character(0)
for (set in INPUT_SETS) {
  inp <- make_input(set)
  inp <- inp[!is.na(inp$Gene.symbol) & inp$Gene.symbol != "" &
             is.finite(inp$logFC) & is.finite(inp$p_val) & inp$p_val > 0 & inp$p_val <= 1, ]
  inp <- inp[!duplicated(inp$Gene.symbol), ]
  cat("\n=====", set, "input:", nrow(inp), "genes ",
      sprintf("(%d up / %d down)", sum(inp$logFC > 0), sum(inp$logFC < 0)), "=====\n")
  fwrite(inp, file.path(outdir, sprintf("pathfindr_input_%s.csv", set)))

  for (db in DBS) {
    od <- file.path(WORK_DIR, sprintf("run_%s_%s", set, gsub("[^A-Za-z0-9]", "", db)))
    if (dir.exists(od)) unlink(od, recursive = TRUE)
    cat("--- ", set, "/", db, " ---\n")
    run_one <- function(np, dir_) pathfindR::run_pathfindR(
      input = inp, gene_sets = db, pin_name_path = "STRING",
      p_val_threshold = 0.05, enrichment_threshold = 0.05,
      min_gset_size = 15, max_gset_size = 1000, convert2alias = FALSE,
      plot_enrichment_chart = FALSE, iterations = ITER,
      output_dir = dir_, n_processes = np)
    res <- try(run_one(N_PROC, od), silent = FALSE)
    if (inherits(res, "try-error")) {
      cat("   parallel run failed; retrying", set, "/", db, "single-threaded\n")
      od <- paste0(od, "_serial")
      if (dir.exists(od)) unlink(od, recursive = TRUE)
      res <- try(run_one(1, od), silent = FALSE)
    }
    if (inherits(res, "try-error") || is.null(res) || nrow(res) == 0) {
      cat("   *** NO RESULT for", set, db, "- this database is missing from the summary ***\n")
      failed <<- c(failed, paste(set, db)); next
    }
    res$database <- db; res$input_set <- set
    all_terms[[paste(set, db)]] <- res
    fwrite(res, file.path(outdir, sprintf("enriched_%s_%s.csv", set, gsub("[^A-Za-z0-9]", "", db))))
    cat("   ", nrow(res), "enriched terms (raw pathfindR output in", od, "); top 5:\n")
    print(head(res[, c("Term_Description", "Fold_Enrichment", "lowest_p")], 5))
  }
}

if (length(failed))
  warning("these input-set / database combinations produced no result: ",
          paste(failed, collapse = "; "),
          " - do not treat the summary below as complete until they are re-run.")
comb <- rbindlist(all_terms, fill = TRUE)
fwrite(comb, file.path(outdir, "pathfindr_all_enriched_terms_reanalysis.csv"))

# ---- did the claims of the original analysis survive? -----------------------
claims <- data.table(
  claim = c("Cell cycle", "Oocyte meiosis", "p53 signaling pathway", "ECM-receptor interaction",
            "Focal adhesion", "chromosome segregation", "mitotic cell cycle", "mitotic cytokinesis",
            "G2/M transition", "kinetochore", "spindle", "midbody", "centrosome",
            "microtubule binding", "protein serine/threonine kinase activity", "protein kinase binding"),
  database = c(rep("KEGG", 5), rep("GO-BP", 4), rep("GO-CC", 4), rep("GO-MF", 3)))
chk <- rbindlist(lapply(INPUT_SETS, function(set) {
  cc <- copy(claims)[, input_set := set]
  cc[, found := mapply(function(p, db) {
        s <- comb[input_set == set & database == db]
        if (!nrow(s)) return(NA)
        any(grepl(p, s$Term_Description, ignore.case = TRUE))
      }, claim, database)]
  cc
}))
fwrite(chk, file.path(outdir, "original_pathway_claims_recheck.csv"))
cat("\n===== do the original Figure 4 pathway claims still hold? =====\n"); print(chk)
cat("\nwritten to", outdir, "\n")
sessionInfo()
