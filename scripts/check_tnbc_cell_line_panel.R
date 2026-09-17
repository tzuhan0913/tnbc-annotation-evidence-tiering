# =============================================================================
# Curate the 24-line TNBC panel that Tier 5 rests on.
#
# The panel in analysis/external_drug_validation/inputs/tnbc_cell_line_panel_template.csv
# is marked "manual_common_TNBC_panel_needs_confirmation" on every row, so it has never
# been verified.  Two independent checks are run here:
#
#   A. IDENTITY - Cellosaurus REST API: accession, species, disease term, and whether the
#      line is flagged as problematic (misidentified / contaminated).  Needs internet.
#   B. RECEPTOR STATUS - ER / PR / HER2 from DepMap expression, which is what actually
#      decides whether a line is triple-negative.  Needs a small DepMap download; see below.
#
# DepMap download (2 minutes, tiny file):
#   1. go to https://depmap.org/portal/data_page/?tab=customDownloads
#   2. Genes: ESR1, PGR, ERBB2   |   Datasets: Expression (log2 TPM+1)
#   3. Cell lines: all   ->  Download
#   4. save the CSV into  data_external/depmap/   (any filename)
# If that folder is empty the script still runs part A and tells you what is missing.
#
#   setwd("C:/Users/tzuha/Documents/投稿文章(JS)")
#   source("scripts/check_tnbc_cell_line_panel.R")
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
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
root   <- find_root(file.path("data", "GSE161529"))
cat("project root:", root, "\n")
panelf <- file.path(root, "analysis", "external_drug_validation", "inputs",
                    "tnbc_cell_line_panel_template.csv")
outdir <- file.path(root, "analysis", "tier5_gdsc_reanalysis")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
panel  <- fread(panelf)
setnames(panel, 1, "cell_line_name")
cat("panel lines:", nrow(panel), "\n")
norm <- function(x) gsub("[^A-Z0-9]+", "", toupper(as.character(x)))

# ---- A. Cellosaurus identity check -----------------------------------------
cello <- data.table(cell_line_name = panel$cell_line_name, accession = NA_character_,
                    cello_name = NA_character_, species = NA_character_,
                    disease = NA_character_, problematic = NA_character_,
                    cellosaurus_url = NA_character_)
have_json <- requireNamespace("jsonlite", quietly = TRUE)
if (!have_json) cat('\n(install.packages("jsonlite") to enable the Cellosaurus check)\n')
for (i in seq_len(nrow(panel))) {
  nm <- panel$cell_line_name[i]
  if (!have_json) next
  q <- utils::URLencode(sprintf('id:"%s"', nm), reserved = TRUE)
  url <- sprintf("https://api.cellosaurus.org/search/cell-line?q=%s&format=json&rows=5", q)
  j <- tryCatch(jsonlite::fromJSON(url, simplifyVector = FALSE), error = function(e) NULL)
  cl <- tryCatch(j[["Cellosaurus"]][["cell-line-list"]], error = function(e) NULL)
  if (is.null(cl) || !length(cl)) { cat(sprintf("  %-12s  not found via the API\n", nm)); next }
  x <- cl[[1]]
  getf <- function(node, key) {
    v <- tryCatch(node[[key]], error = function(e) NULL)
    if (is.null(v)) return(NA_character_)
    paste(unlist(lapply(v, function(e) if (is.list(e)) unlist(e)[1] else e)), collapse = "; ")
  }
  ac <- tryCatch(unlist(x[["accession-list"]])[["value"]], error = function(e) NA_character_)
  if (is.na(ac)) ac <- tryCatch(unlist(x[["accession-list"]])[1], error = function(e) NA_character_)
  cello[i, `:=`(accession = ac,
                cello_name = getf(x, "name-list"),
                species    = getf(x, "species-list"),
                disease    = getf(x, "disease-list"),
                problematic = as.character(!is.null(x[["problematic"]]) && isTRUE(x[["problematic"]])),
                cellosaurus_url = if (!is.na(ac)) paste0("https://www.cellosaurus.org/", ac) else NA_character_)]
  cat(sprintf("  %-12s  %-12s  %s\n", nm, ac, substr(cello$disease[i], 1, 60)))
}

# ---- B. receptor status from DepMap expression ------------------------------
## Search the whole of data_external, not just data_external/depmap: the custom download lands
## wherever the browser was pointed, and its file name varies with the release.
ddir <- file.path(root, "data_external")
dir.create(file.path(ddir, "depmap"), showWarnings = FALSE, recursive = TRUE)
files <- list.files(ddir, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE,
                    recursive = TRUE)
expr <- NULL
for (f in files) {
  d <- tryCatch(suppressWarnings(fread(f, nrows = 5)), error = function(e) NULL)
  if (is.null(d) || !length(grep("ESR1", names(d), ignore.case = TRUE))) next
  expr <- fread(f); cat("\nusing DepMap file:", sub(root, "", f, fixed = TRUE), "\n"); break
}
res <- merge(panel, cello, by = "cell_line_name", all.x = TRUE, sort = FALSE)
if (is.null(expr)) {
  cat("\nNo DepMap expression file found in", ddir, "\n",
      "-> download ESR1 / PGR / ERBB2 expression as described in the header and re-run.\n")
  res[, `:=`(ESR1 = NA_real_, PGR = NA_real_, ERBB2 = NA_real_, receptor_call = NA_character_)]
} else {
  gcol <- function(g) { h <- grep(paste0("^", g, "\\b|^", g, " "), names(expr), ignore.case = TRUE)
                        if (length(h)) h[1] else NA_integer_ }
  ## DepMap custom downloads identify rows by ModelID (ACH-...) and, when the metadata box is
  ## ticked, also by one or more name columns.  Rather than guess which column holds the cell
  ## line name, score every character column by how many panel lines it matches and take the best.
  panel_key <- unique(norm(panel$cell_line_name))
  cand <- names(expr)[sapply(expr, function(c) is.character(c) || is.factor(c))]
  if (!length(cand)) cand <- names(expr)[1]
  hits <- sapply(cand, function(cc) sum(unique(norm(expr[[cc]])) %in% panel_key))
  key  <- cand[which.max(hits)]
  cat("cell-line key column:", key, "- matches", max(hits), "of", length(panel_key), "panel lines\n")
  if (max(hits) == 0)
    warning("No column in the DepMap file matches the panel cell-line names. Re-download with\n",
            '  "Add cell line metadata to download" ticked, so the file carries CellLineName.')
  expr[, .k := norm(get(key))]
  idx <- sapply(c("ESR1", "PGR", "ERBB2"), gcol)
  cat("matched expression columns:", paste(names(expr)[idx[!is.na(idx)]], collapse = ", "), "\n")
  ## A missing receptor is not a passing receptor.  The custom-download form silently drops a gene
  ## if it was typed on the same line as another, so say so instead of calling lines TN by default.
  if (any(is.na(idx)))
    warning("the DepMap file has no column for: ", paste(names(idx)[is.na(idx)], collapse = ", "),
            ".\n  Lines cannot be excluded on that receptor - re-download with one gene per line ",
            "and re-run before treating the curated panel as final.")
  e <- expr[, c(".k", names(expr)[idx[!is.na(idx)]]), with = FALSE]
  setnames(e, c(".k", names(idx)[!is.na(idx)]))
  res[, .k := norm(cell_line_name)]
  res <- merge(res, e, by = ".k", all.x = TRUE, sort = FALSE)
  # log2(TPM+1) thresholds, calibrated on this release rather than taken from habit.
  # ESR1: the ER+ reference lines (MCF7 5.62, T47D 5.64, BT474 4.43) sit well above every line
  #       in the panel (max 2.67), so 3 separates them cleanly.
  # ERBB2: the 14 HER2-amplified breast lines in DepMap all sit at >= 9.65 (HCC1954 11.94,
  #       SKBR3 10.80, BT474 10.47) while HER2-negative MCF7 is 5.40 and the all-line 95th
  #       percentile is 6.5.  A threshold of 7, which is often quoted, falls inside the
  #       non-amplified range and misclassifies ordinary TNBC lines; 9 sits in the empty gap.
  #       Lines between 7 and 9 (MDA-MB-453 at 8.06) are flagged as borderline, not excluded.
  THR <- c(ESR1 = 3, PGR = 2, ERBB2 = 9)
  ERBB2_BORDERLINE <- 7
  if (!"ERBB2" %in% names(res)) res[, ERBB2 := NA_real_]
  if (!"PGR"   %in% names(res)) res[, PGR   := NA_real_]
  if (!"ESR1"  %in% names(res)) res[, ESR1  := NA_real_]
  res[, receptors_assessed := paste(names(idx)[!is.na(idx)], collapse = "/")]
  res[, receptor_call := fifelse(is.na(ESR1) & is.na(PGR) & is.na(ERBB2), "not in DepMap",
       fifelse((!is.na(ESR1) & ESR1 >= THR["ESR1"]) | (!is.na(PGR) & PGR >= THR["PGR"]) |
               (!is.na(ERBB2) & ERBB2 >= THR["ERBB2"]), "NOT triple-negative", "triple-negative"))]
  res[, .k := NULL]
}

if ("ERBB2" %in% names(res))
  res[, erbb2_note := fifelse(!is.na(ERBB2) & ERBB2 >= ERBB2_BORDERLINE & ERBB2 < THR["ERBB2"],
                              "borderline HER2 (7-9); kept, but test Tier 5 without it", "")]
res[, flag := fifelse(!is.na(problematic) & problematic == "TRUE", "Cellosaurus: problematic",
      fifelse(!is.na(receptor_call) & receptor_call == "NOT triple-negative", "receptor-positive",
      fifelse(!is.na(receptor_call) & receptor_call == "not in DepMap", "no DepMap data", "")))]
res[, include := flag == ""]
fwrite(res, file.path(outdir, "tnbc_panel_curation.csv"))
cat("\n===== panel curation =====\n")
print(res[, .(cell_line_name, accession,
              ESR1 = round(ESR1, 2), PGR = round(PGR, 2), ERBB2 = round(ERBB2, 2),
              receptor_call, flag, include,
              note = if ("erbb2_note" %in% names(res)) erbb2_note else "")])
if ("erbb2_note" %in% names(res) && any(nzchar(res$erbb2_note)))
  cat("\nBorderline HER2 (kept in the panel, but worth a Tier 5 sensitivity run without them): ",
      paste(res$cell_line_name[nzchar(res$erbb2_note)], collapse = ", "), "\n", sep = "")
cat(sprintf("\n%d of %d lines pass both checks.\n", sum(res$include, na.rm = TRUE), nrow(res)))
## The curated panel is written EVERY time, even when nothing is excluded.  An earlier version
## wrote it only on failure, so a stale file from a previous run with a different threshold stayed
## on disk and was silently picked up by Tier 5 - exactly the kind of invisible inherited input
## this project exists to complain about.
pin <- file.path(root, "analysis", "external_drug_validation", "inputs")
keep <- function(f, d, label) {
  fwrite(d[, .(cell_line_name, subtype = "TNBC", curation_status = label)], file.path(pin, f))
  cat(sprintf("  %-38s %2d lines\n", f, nrow(d)))
}
cat("panel files written to", pin, ":\n")
keep("tnbc_cell_line_panel_curated.csv", res[include == TRUE], "verified_cellosaurus_depmap")
## The strict panel additionally drops lines whose ERBB2 sits in the 7-9 borderline band, where the
## HER2 call is genuinely contested.  Tier 5 is run on both; if the two disagree, that disagreement
## is a result and belongs in the paper, not something to resolve by picking the nicer panel.
strict <- res[include == TRUE &
              (!("erbb2_note" %in% names(res)) | !nzchar(erbb2_note))]
keep("tnbc_cell_line_panel_strict.csv", strict, "verified_cellosaurus_depmap_strict_HER2")
if (any(!res$include, na.rm = TRUE))
  cat("Excluded outright: ", paste(res$cell_line_name[!res$include], collapse = ", "), "\n", sep = "")
cat('\nNext: run scripts/run_tier5_gdsc.R twice -\n',
    '  PANEL <- "full"    (', sum(res$include, na.rm = TRUE), ' lines) - the primary analysis\n',
    '  PANEL <- "strict"  (', nrow(strict), ' lines) - the sensitivity analysis\n', sep = "")
sessionInfo()
