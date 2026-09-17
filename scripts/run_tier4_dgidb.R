# =============================================================================
# Tier 4 re-run: DGIdb drug-gene matching on the harmonised labels.
#
# The original Table 1 came from a DGIdb v4.4 query on the 869-gene cell-level list.
# This script repeats the matching on the Tier 2 patient-level gene list (default), so
# no drug hypothesis inherits a cell-level gene list, and reports which of the original
# 27 compounds survive.
#
# Data: DGIdb publishes a bulk interactions table.  The script looks for a local copy at
#   data_external/dgidb/interactions.tsv
# and, if it is absent, tries to download one.  If the download fails, download it by hand
# from https://www.dgidb.org/downloads (the "interactions" TSV) and put it at that path,
# then run this script again.  Record the version you used - it is not the 2021 v4.4 release.
#
#   setwd("C:/Users/tzuha/Documents/投稿文章(JS)")
#   source("scripts/run_tier4_dgidb.R")
# =============================================================================
GENE_SET  <- "patient"     # "patient" (Tier 2 genes) | "cell" (Tier 1 genes) | "legacy"
MIN_GENES <- 2             # a compound is kept when it interacts with >= this many genes
TOP_N_CELL <- 1000

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
indir  <- file.path(root, "analysis", "epithelial_only_pseudobulk_R")
outdir <- file.path(root, "analysis", "tier4_dgidb_reanalysis")
ddir   <- file.path(root, "data_external", "dgidb")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(ddir,   showWarnings = FALSE, recursive = TRUE)

# ---- 1. the DGIdb interactions table ---------------------------------------
ipath <- file.path(ddir, "interactions.tsv")
if (!file.exists(ipath)) {
  urls <- c("https://dgidb.org/data/latest/interactions.tsv",
            "https://www.dgidb.org/data/latest/interactions.tsv",
            "https://dgidb.org/data/monthly_tsvs/2022-Feb/interactions.tsv")
  for (u in urls) {
    cat("trying", u, "\n")
    ok <- tryCatch({ download.file(u, ipath, mode = "wb", quiet = FALSE); file.exists(ipath) &&
                     file.size(ipath) > 100000 }, error = function(e) FALSE)
    if (isTRUE(ok)) { cat("downloaded from", u, "\n"); break }
    if (file.exists(ipath)) unlink(ipath)
  }
}
if (!file.exists(ipath))
  stop("Could not obtain the DGIdb interactions table.\n",
       "Download it by hand from https://www.dgidb.org/downloads and save it as:\n  ", ipath)
dg <- fread(ipath, sep = "\t", quote = "")
cat("DGIdb rows:", nrow(dg), "| columns:", paste(names(dg), collapse = ", "), "\n")
writeLines(c(paste("file:", ipath),
             paste("downloaded/used:", format(Sys.time())),
             paste("rows:", nrow(dg)),
             paste("columns:", paste(names(dg), collapse = ", "))),
           file.path(outdir, "dgidb_source_version.txt"))

pick <- function(cands) { h <- intersect(cands, names(dg)); if (length(h)) h[1] else NA_character_ }
c_gene  <- pick(c("gene_name", "gene_claim_name", "gene_concept_id"))
c_drug  <- pick(c("drug_name", "drug_claim_name", "drug_claim_primary_name"))
c_score <- pick(c("interaction_score", "interaction_group_score"))
c_appr  <- pick(c("approved"))
c_anti  <- pick(c("anti_neoplastic"))
c_type  <- pick(c("interaction_type", "interaction_types"))
if (is.na(c_gene) || is.na(c_drug)) stop("Unrecognised DGIdb layout: ", paste(names(dg), collapse = ", "))
cat("using columns -> gene:", c_gene, "| drug:", c_drug, "| score:", c_score, "\n")

# ---- 2. the query gene list ------------------------------------------------
if (GENE_SET == "patient") {
  d <- fread(file.path(indir, "tier2a_paired_edgeR_QL_min10.csv"))[FDR < 0.05 & abs(logFC) > 0.5]
  genes <- data.table(gene = d$gene, effect = d$logFC, padj = d$FDR)
} else if (GENE_SET == "cell") {
  d <- fread(file.path(indir, "tier1_cell_level_findmarkers.csv"))[p_val_adj < 0.05 & abs(avg_log2FC) > 0.5][order(p_val)]
  if (nrow(d) > TOP_N_CELL) d <- d[seq_len(TOP_N_CELL)]
  genes <- data.table(gene = d$gene, effect = d$avg_log2FC, padj = d$p_val_adj)
} else {
  d <- fread(file.path(root, "inputs", "liao_master_zip_20260707",
                       "05-家嫺論文補分析", "Data", "data2.csv"))
  setnames(d, trimws(names(d)))
  genes <- data.table(gene = trimws(d$Gene), effect = d$avg_log2FC, padj = d$p_val_adj)
}
genes <- genes[!is.na(gene) & gene != ""][!duplicated(gene)]
cat("query genes (", GENE_SET, "):", nrow(genes), "\n")
fwrite(genes, file.path(outdir, sprintf("query_genes_%s.csv", GENE_SET)))

# ---- 3. match --------------------------------------------------------------
sel <- c(c_gene, c_drug, c_score, c_appr, c_anti, c_type)
sel <- sel[!is.na(sel)]
ed  <- unique(dg[get(c_gene) %in% genes$gene, ..sel])
setnames(ed, c_gene, "gene"); setnames(ed, c_drug, "drug")
## DGIdb ships rows whose drug_name is the literal string "NULL".  In this query they amount to 91
## gene interactions, which is more than any real compound, so an unfiltered multi-target ranking
## puts a non-existent drug at the top of the candidate list.  Drop them and say how many.
n_null <- sum(toupper(trimws(ed$drug)) %in% c("NULL", "NA", "NONE"))
if (n_null) cat("dropping", n_null, "edges whose drug name is NULL/NA - a DGIdb artefact\n")
ed <- ed[!is.na(drug) & trimws(drug) != "" &
         !(toupper(trimws(drug)) %in% c("NULL", "NA", "NONE"))]
if (!is.na(c_score)) setnames(ed, c_score, "interaction_score") else ed[, interaction_score := NA_real_]
ed[, interaction_score := suppressWarnings(as.numeric(interaction_score))]
if (!is.na(c_appr)) setnames(ed, c_appr, "approved") else ed[, approved := NA]
if (!is.na(c_anti)) setnames(ed, c_anti, "anti_neoplastic") else ed[, anti_neoplastic := NA]
ed <- merge(ed, genes, by = "gene", all.x = TRUE)
cat("drug-gene edges:", nrow(ed), "| distinct drugs:", uniqueN(ed$drug), "\n")

drug <- ed[, .(n_genes = uniqueN(gene),
               genes = paste(sort(unique(gene)), collapse = "; "),
               max_score = suppressWarnings(max(interaction_score, na.rm = TRUE)),
               sum_score = sum(interaction_score, na.rm = TRUE),
               approved = any(tolower(as.character(approved)) %in% c("true", "yes", "1")),
               anti_neoplastic = any(tolower(as.character(anti_neoplastic)) %in% c("true", "yes", "1")),
               mean_effect = mean(effect, na.rm = TRUE)), by = drug]
drug[!is.finite(max_score), max_score := NA_real_]
multi <- drug[n_genes >= MIN_GENES][order(-approved, -anti_neoplastic, -n_genes, -sum_score)]
fwrite(ed[order(drug, gene)], file.path(outdir, sprintf("dgidb_edges_%s.csv", GENE_SET)))
fwrite(drug[order(-n_genes, -sum_score)], file.path(outdir, sprintf("dgidb_drug_summary_%s.csv", GENE_SET)))
fwrite(multi, file.path(outdir, sprintf("dgidb_multitarget_%s.csv", GENE_SET)))
cat("multi-target compounds (>=", MIN_GENES, "genes):", nrow(multi), "\n")
print(head(multi[, .(drug, n_genes, approved, anti_neoplastic, genes)], 25))

# ---- 4. what happened to the original 27 compounds? ------------------------
oldf <- file.path(root, "analysis", "external_drug_validation", "inputs",
                  "dgidb_candidate_drugs_for_external_validation.csv")
if (file.exists(oldf)) {
  old <- fread(oldf)
  nm  <- intersect(c("drug_name", "drug", "Drug", "Drug name"), names(old))[1]
  oldd <- toupper(trimws(old[[nm]]))
  newd <- toupper(trimws(multi$drug))
  cmp <- data.table(original_drug = sort(unique(oldd)))
  cmp[, retained_in_reanalysis := original_drug %in% newd]
  cmp[, still_in_dgidb_at_all := original_drug %in% toupper(trimws(drug$drug))]
  fwrite(cmp, file.path(outdir, "original_27_compounds_recheck.csv"))
  cat("\n===== original candidate compounds =====\n")
  cat(sprintf("%d of %d retained as multi-target hits under the %s gene set\n",
              sum(cmp$retained_in_reanalysis), nrow(cmp), GENE_SET))
  print(cmp)
} else {
  cat("\n(original candidate-drug list not found at", oldf, "- skipping the comparison)\n")
}

# ---- 5. hand-off for Tier 5 ------------------------------------------------
fwrite(multi[, .(drug_name = drug, n_genes, genes, approved, anti_neoplastic, max_score, sum_score)],
       file.path(outdir, "candidate_drugs_for_external_validation.csv"))
cat("\nwritten to", outdir, "\n")
cat("next: send analysis/tier4_dgidb_reanalysis/candidate_drugs_for_external_validation.csv for the GDSC step.\n")
sessionInfo()
