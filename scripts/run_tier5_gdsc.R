# =============================================================================
# Tier 5 re-run: GDSC pharmacogenomic cross-reference for the harmonised mitotic axis.
#
# Two things change relative to the original target-axis analysis:
#   1. the axis is the full set of twelve mitotic regulators supported at the patient
#      level, not only AURKA / AURKB / PLK1;
#   2. the TNBC-favourable direction is calibrated against a background of all other
#      GDSC compounds, and tested, instead of being reported as a bare count.
#      Without that background, "7 of 12 compounds look TNBC-favourable" is
#      uninterpretable: roughly half of any compound set will, by construction.
#
# Inputs already present in data_external/gdsc:
#   sanger-dose-response.csv, Cell_Lines_Details.xlsx, screened_compounds_rel_8.5.csv
# This script does NOT depend on Tier 4.  The target-axis analysis - which is the part that
# changes in R12 - uses only GDSC's own TARGET annotation, so it can be run before or without
# run_tier4_dgidb.R.  If analysis/tier4_dgidb_reanalysis/candidate_drugs_for_external_validation.csv
# happens to exist, an extra exact-drug section is appended; if not, that section is skipped
# and everything else still runs.
#
#   setwd("C:/Users/tzuha/Documents/投稿文章(JS)")
#   source("scripts/run_tier5_gdsc.R")
# =============================================================================
MITOTIC <- c("PLK1", "AURKA", "AURKB", "TTK", "NEK2", "MELK", "CDK1", "CCNB1",
             "CCNB2", "GTSE1", "TOP2A", "MKI67")
# GDSC writes targets with its own vocabulary; these alternates map onto the same axis
AXIS_PATTERNS <- c("AURKA", "AURKB", "AURORA", "PLK1", "PLK ", "^PLK", "TTK", "MPS1",
                   "NEK2", "MELK", "\\bCDK1\\b", "TOP2A", "TOP2")
MIN_PROFILES <- 3     # a compound needs at least this many TNBC-panel profiles to count

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
root <- find_root(file.path("data", "GSE161529"))
cat("project root:", root, "\n")
ddir   <- file.path(root, "data_external", "gdsc")
outdir <- file.path(root, "analysis", "tier5_gdsc_reanalysis")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
for (f in c("sanger-dose-response.csv", "Cell_Lines_Details.xlsx", "screened_compounds_rel_8.5.csv"))
  if (!file.exists(file.path(ddir, f))) stop("missing ", file.path(ddir, f))
if (!requireNamespace("readxl", quietly = TRUE)) stop('install.packages("readxl")')

norm <- function(x) gsub("[^A-Z0-9]+", "", toupper(as.character(x)))

resp  <- fread(file.path(ddir, "sanger-dose-response.csv"))
comp  <- fread(file.path(ddir, "screened_compounds_rel_8.5.csv"))
cells <- as.data.table(readxl::read_excel(file.path(ddir, "Cell_Lines_Details.xlsx")))
# Which TNBC cell-line panel to use.  Run the script once with each value: the axis result is
# sensitive to two contested lines, so both numbers belong in the paper.
#   PANEL <- "full"    the lines that pass the Cellosaurus and ER/PR/HER2 checks (24)
#   PANEL <- "strict"  the same, minus the lines whose ERBB2 sits in the 7-9 borderline band
#                      (HCC38 and MDA-MB-453) - the sensitivity analysis
# Both files are written by scripts/check_tnbc_cell_line_panel.R.  Outputs from the strict run are
# suffixed "_strictpanel" so neither run overwrites the other.
# Set PANEL in the console before sourcing to override the default, e.g.
#   PANEL <- "strict"; source("scripts/run_tier5_gdsc.R")
if (!exists("PANEL")) PANEL <- "full"
stopifnot(PANEL %in% c("full", "strict"))
pdir  <- file.path(root, "analysis", "external_drug_validation", "inputs")
pfile <- if (PANEL == "strict") {
  file.path(pdir, "tnbc_cell_line_panel_strict.csv")
} else if (file.exists(file.path(pdir, "tnbc_cell_line_panel_curated.csv"))) {
  file.path(pdir, "tnbc_cell_line_panel_curated.csv")
} else {
  file.path(pdir, "tnbc_cell_line_panel_template.csv")
}
if (!file.exists(pfile))
  stop("panel file not found: ", pfile, "\n  run scripts/check_tnbc_cell_line_panel.R first.")
SFX <- if (PANEL == "strict") "_strictpanel" else ""
out <- function(f) file.path(outdir, sub("\\.csv$", paste0(SFX, ".csv"), f))
cat("TNBC panel:", PANEL, "->", basename(pfile), "\n")
panel <- fread(pfile)
setnames(panel, 1, "cell_line_name")
cat("dose-response rows:", nrow(resp), "| compounds:", nrow(comp), "| cell lines:", nrow(cells), "\n")

cat("dose-response columns:", paste(names(resp), collapse = ", "), "\n")
cn <- names(cells)
setnames(cells, cn[grep("^COSMIC", cn)][1], "COSMIC_ID")
setnames(cells, cn[grep("^Sample Name", cn)][1], "CELL_LINE_ANNOT")
setnames(cells, cn[grep("TCGA", cn)][1], "TCGA_LABEL")
setnames(cells, cn[grep("Tissue", cn)][1], "GDSC_TISSUE_1")
cells[, COSMIC_ID := suppressWarnings(as.numeric(COSMIC_ID))]
resp[, COSMIC_ID := suppressWarnings(as.numeric(COSMIC_ID))]
for (v in intersect(c("TCGA_LABEL", "GDSC_TISSUE_1", "TARGET", "TARGET_PATHWAY", "is_axis"), names(resp)))
  resp[, (v) := NULL]                       # avoid .x/.y collisions on the merges below
resp <- merge(resp, cells[, .(COSMIC_ID, CELL_LINE_ANNOT, TCGA_LABEL, GDSC_TISSUE_1)],
              by = "COSMIC_ID", all.x = TRUE)
resp[, cell_name := if ("CELL_LINE_NAME" %in% names(resp)) as.character(CELL_LINE_NAME) else NA_character_]
resp[is.na(cell_name) | cell_name == "", cell_name := as.character(CELL_LINE_ANNOT)]

panel_norm <- unique(norm(panel$cell_line_name))
resp[, is_tnbc_panel := norm(cell_name) %in% panel_norm]
resp[, is_breast := grepl("BRCA", as.character(TCGA_LABEL), ignore.case = TRUE) |
                    grepl("breast", as.character(GDSC_TISSUE_1), ignore.case = TRUE)]
cat("TNBC-panel cell lines matched in GDSC:",
    uniqueN(resp[is_tnbc_panel == TRUE]$cell_name), "of", length(panel_norm), "\n")

comp[, TARGET := as.character(TARGET)]
comp[, target_up := toupper(ifelse(is.na(TARGET), "", TARGET))]
comp[, is_axis := Reduce(`|`, lapply(AXIS_PATTERNS, function(p) grepl(p, target_up)))]
comp[, axis_subclass := fifelse(!is_axis, NA_character_,
        fifelse(grepl("TOP2", target_up) & !grepl("AURK|PLK|TTK|MPS1|NEK2|MELK|CDK1", target_up),
                "topoisomerase II (cytotoxic)", "targeted mitotic / cell-cycle"))]
axis_ids <- unique(comp[is_axis == TRUE]$DRUG_ID)
cat("GDSC compounds whose TARGET mentions the mitotic axis:", length(axis_ids), "\n")
resp <- merge(resp, comp[, .(DRUG_ID, TARGET, TARGET_PATHWAY, is_axis, axis_subclass)], by = "DRUG_ID", all.x = TRUE)
resp[is.na(is_axis), is_axis := FALSE]

# ---- per-compound TNBC vs the rest -----------------------------------------
med <- function(x) { v <- suppressWarnings(as.numeric(x)); v <- v[is.finite(v)]
                     if (length(v)) median(v) else NA_real_ }
wilx <- function(a, b) {
  a <- suppressWarnings(as.numeric(a)); b <- suppressWarnings(as.numeric(b))
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < MIN_PROFILES || length(b) < MIN_PROFILES) return(NA_real_)
  suppressWarnings(wilcox.test(a, b, alternative = "less")$p.value)
}
per_drug <- resp[, {
  t <- .SD[is_tnbc_panel == TRUE]; o <- .SD[is_tnbc_panel == FALSE]
  .(targets = paste(sort(unique(na.omit(TARGET))), collapse = "; "),
    pathway = paste(sort(unique(na.omit(TARGET_PATHWAY))), collapse = "; "),
    is_axis = any(is_axis), subclass = paste(sort(unique(na.omit(axis_subclass))), collapse = "; "), n_all = .N, n_tnbc = nrow(t), n_other = nrow(o),
    med_z_tnbc = med(t$Z_SCORE_PUBLISHED), med_z_other = med(o$Z_SCORE_PUBLISHED),
    med_auc_tnbc = med(t$AUC_PUBLISHED), med_auc_other = med(o$AUC_PUBLISHED),
    p_z_tnbc_lower = wilx(t$Z_SCORE_PUBLISHED, o$Z_SCORE_PUBLISHED))
}, by = .(DRUG_ID, DRUG_NAME)]
per_drug[, z_favours_tnbc := med_z_tnbc < med_z_other]
per_drug[, auc_favours_tnbc := med_auc_tnbc < med_auc_other]
ev <- per_drug[n_tnbc >= MIN_PROFILES & !is.na(z_favours_tnbc)]
fwrite(per_drug[order(-is_axis, p_z_tnbc_lower)], out("gdsc_per_compound_summary.csv"))
fwrite(resp[is_axis == TRUE], out("gdsc_axis_cell_line_profiles.csv"))

# ---- collapse repeated screens to the molecule ------------------------------
# GDSC screens some compounds more than once under different DRUG_IDs (tozasertib, BI-2536 and
# alisertib are all screened twice in the axis set).  Counting those screens separately treats a
# technical replicate as an independent compound - the same pseudoreplication this paper is about -
# and it inflates both the candidate set and the background.  The unit of analysis here is therefore
# the MOLECULE.  Three collapse rules are reported, because there is no single obviously correct one:
#   any    - favourable if any screen of that molecule favoured the TNBC panel (most generous)
#   all    - favourable only if every screen did (most conservative)
#   meanz  - favourable if the mean of (median TNBC z - median other z) across screens is negative
# The primary analysis is "meanz"; the other two are the sensitivity analysis.
ev[, molkey := norm(DRUG_NAME)]
ev[, dz := med_z_tnbc - med_z_other]
molt <- ev[, .(n_screens = .N, is_axis = any(is_axis),
              fav_any  = any(z_favours_tnbc), fav_all = all(z_favours_tnbc),
              fav_meanz = mean(dz, na.rm = TRUE) < 0,
              mean_dz = mean(dz, na.rm = TRUE),
              discordant = uniqueN(z_favours_tnbc) > 1), by = molkey]
cat(sprintf("\nevaluable screens: %d -> molecules: %d (%d screened more than once, %d discordant)\n",
            nrow(ev), nrow(molt), sum(molt$n_screens > 1), sum(molt$discordant)))
fwrite(molt[order(-is_axis, mean_dz)], out("gdsc_molecule_level_calls.csv"))

calib <- function(sel, label, rule) {
  f <- molt[[paste0("fav_", rule)]]
  k_s <- sum(f[sel]); n_s <- sum(sel)
  k_b <- sum(f[!sel]); n_b <- sum(!sel)
  if (!n_s || !n_b) return(NULL)
  bg_rate <- k_b / n_b
  ftg <- fisher.test(matrix(c(k_s, n_s - k_s, k_b, n_b - k_b), nrow = 2), alternative = "greater")
  ft2 <- fisher.test(matrix(c(k_s, n_s - k_s, k_b, n_b - k_b), nrow = 2))
  mw  <- suppressWarnings(wilcox.test(molt$mean_dz[sel], molt$mean_dz[!sel], alternative = "less"))
  data.table(set = label, collapse_rule = rule, n_molecules = n_s, favourable = k_s,
             fraction = round(k_s / n_s, 3), background_n = n_b, background_favourable = k_b,
             background_fraction = round(bg_rate, 3),
             binom_vs_half = signif(binom.test(k_s, n_s, 0.5, alternative = "greater")$p.value, 3),
             binom_vs_background = signif(binom.test(k_s, n_s, bg_rate, alternative = "greater")$p.value, 3),
             fisher_one_sided = signif(ftg$p.value, 3), fisher_two_sided = signif(ft2$p.value, 3),
             wilcoxon_on_z = signif(mw$p.value, 3))
}
sets <- list("mitotic target axis" = molt$is_axis)
res <- rbindlist(lapply(names(sets), function(s)
        rbindlist(lapply(c("meanz", "any", "all"), function(r) calib(sets[[s]], s, r)))), fill = TRUE)
fwrite(res, out("gdsc_axis_vs_background_calibration.csv"))
cat("\n===== mitotic axis vs the rest of GDSC (molecule level) =====\n"); print(res)
cat("\nNOTE: the primary row is collapse_rule = 'meanz'. If the three rules disagree about\n",
    "significance, report the range and describe the result as borderline, not as supported.\n")
ax <- ev[is_axis == TRUE]; bg <- ev[is_axis == FALSE]
cat("\n--- evaluable axis compounds ---\n")
print(ax[order(p_z_tnbc_lower), .(DRUG_NAME, targets, subclass, n_tnbc,
       med_z_tnbc = round(med_z_tnbc, 3), med_z_other = round(med_z_other, 3),
       p = signif(p_z_tnbc_lower, 3))])

# ---- optional: exact candidate drugs from the Tier 4 re-run -----------------
cf <- file.path(root, "analysis", "tier4_dgidb_reanalysis", "candidate_drugs_for_external_validation.csv")
if (file.exists(cf)) {
  cand <- fread(cf)
  cand[, key := norm(drug_name)]
  per_drug[, key := norm(DRUG_NAME)]
  hit <- merge(cand, per_drug, by = "key", all.x = TRUE)
  hit[, in_gdsc := !is.na(DRUG_ID)]
  fwrite(hit, out("gdsc_exact_candidate_drug_match.csv"))
  cat(sprintf("\nexact-drug coverage: %d of %d Tier 4 candidates found in GDSC; %d evaluable in the TNBC panel\n",
              sum(hit$in_gdsc), nrow(hit), sum(hit$in_gdsc & hit$n_tnbc >= MIN_PROFILES, na.rm = TRUE)))
  print(hit[in_gdsc == TRUE, .(drug_name, n_tnbc, med_z_tnbc = round(med_z_tnbc, 3),
                               med_z_other = round(med_z_other, 3), z_favours_tnbc)])

  # --- the same calibration for the Tier 4 drug list, at the molecule level ---
  # This is the comparison that decides whether the drug list, as opposed to the target axis,
  # survives external data.  It must use the same unit (molecule) and the same background as above.
  cand_mol <- unique(cand$key)
  sel <- molt$molkey %in% cand_mol
  appr <- unique(cand[approved == TRUE & anti_neoplastic == TRUE]$key)
  sel2 <- molt$molkey %in% appr
  res4 <- rbindlist(c(
    lapply(c("meanz", "any", "all"), function(r) calib(sel,  "Tier 4 drug list", r)),
    lapply(c("meanz", "any", "all"), function(r) calib(sel2, "Tier 4, approved & anti-neoplastic", r))
  ), fill = TRUE)
  fwrite(rbindlist(list(res, res4), fill = TRUE),
         out("tier5_target_axis_vs_drug_list.csv"))
  cat("\n===== target axis vs the Tier 4 drug list, molecule level =====\n")
  print(rbindlist(list(res, res4), fill = TRUE)[, .(set, collapse_rule, n_molecules, favourable,
                                                    fraction, background_fraction,
                                                    fisher_one_sided, fisher_two_sided)])

  # per-compound multiplicity correction, one test per molecule
  pc <- ev[molkey %in% cand_mol, .(p = min(p_z_tnbc_lower, na.rm = TRUE)), by = .(molkey, DRUG_NAME)]
  pc <- pc[is.finite(p)][order(p)]
  pc[, BH := p.adjust(p, "BH")]
  fwrite(pc, out("tier5_candidate_per_compound_BH.csv"))
  cat(sprintf("\nper-molecule BH across %d candidate molecules; %d at BH < 0.05\n",
              nrow(pc), sum(pc$BH < 0.05)))
  print(head(pc, 5))
} else {
  cat("\n(no Tier 4 candidate list yet - run scripts/run_tier4_dgidb.R first for the exact-drug part)\n")
}
cat("\nwritten to", outdir, "\n")
sessionInfo()
