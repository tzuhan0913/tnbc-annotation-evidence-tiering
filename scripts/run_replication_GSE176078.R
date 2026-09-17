# =============================================================================
# Compact replication in GSE176078 (Wu et al. 2021)
#
# Purpose:
#   Replicate the methodological risk tested in the BioData Mining manuscript:
#   malignant-cell annotation is a zeroth evidence tier in single-cell
#   drug-target studies.  This is a compact sanity check, not a new biological
#   study.
#
# Analyses performed, in order:
#   A. Annotation agreement between authors' epithelial labels and an inferCNV
#      label set derived here.
#   B. Unit of replication under both label sets: cell-level FindMarkers and
#      patient-level edgeR quasi-likelihood pseudo-bulk.
#   C. G2/M program versus single-gene stability under both label sets.
#   D. Detection-rate diagnostic for the top 200 cell-level genes.
#
# Explicitly not performed:
#   No pathway enrichment, no DGIdb, no GDSC, no perturbation signatures.
#
# Expected input:
#   Download and extract the GEO processed scRNA-seq supplement:
#     https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
#
#   The extracted directory should contain:
#     count_matrix_sparse.mtx
#     count_matrix_barcodes.tsv
#     count_matrix_genes.tsv
#     metadata.csv
#
#   Place the extracted files under one of:
#     data/GSE176078/Wu_etal_2021_BRCA_scRNASeq/
#     data/GSE176078/
#     data_external/GSE176078/Wu_etal_2021_BRCA_scRNASeq/
#     data_external/GSE176078/
#
# Run from the project root:
#   setwd("C:/Users/tzuha/Documents/投稿文章(JS)")
#   source("scripts/run_replication_GSE176078.R")
# =============================================================================

MINC <- 10
REF_CAP_PER_PATIENT_CELLTYPE <- 400
INFERCNV_CUTOFF <- 0.1
INFERCNV_ANALYSIS_MODE <- "samples"
INFERCNV_CLUSTER_BY_GROUPS <- TRUE
INFERCNV_DENOISE <- TRUE
INFERCNV_HMM <- FALSE
SKIP_SUBCLUSTER <- TRUE
CNV_CORR_MALIGNANT_THRESHOLD <- 0.5
DETECTION_RATE_LOW_THRESHOLD <- 0.01
TOP_N_DETECTION <- 200

KEY_MITOTIC <- c("PLK1", "AURKA", "AURKB", "TTK", "NEK2", "MELK", "CDK1",
                 "CCNB1", "CCNB2", "GTSE1", "TOP2A", "MKI67")
KEY_IG <- c("IGKC", "IGLC2", "IGHG1", "IGHA1")
KEY <- c(KEY_MITOTIC, KEY_IG)

G2M <- c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80",
         "CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","PIMREG",
         "FAM64A","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11",
         "ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","JPT1","HN1",
         "CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5",
         "CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN",
         "LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5",
         "CENPA")

# Fallback only if the deposited metadata does not include a clinical-subtype
# column.  These are the ten TNBC samples listed for GSE176078/Wu 2021 in GEO
# and widely reused downstream; the script records when this fallback is used.
FALLBACK_TNBC_SAMPLES <- c("CID3946", "CID3963", "CID44041", "CID4465",
                           "CID4495", "CID44971", "CID44991", "CID4513",
                           "CID4515", "CID4523")

suppressPackageStartupMessages({
  library(Matrix)
})

find_root <- function(marker) {
  p <- getwd()
  for (i in 0:5) {
    if (dir.exists(file.path(p, marker))) {
      return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
    p <- file.path(p, "..")
  }
  stop("Cannot find '", marker, "' from ", getwd(),
       '\n  setwd("C:/Users/tzuha/Documents/投稿文章(JS)") and source again.')
}

root <- find_root(file.path("data", "GSE161529"))
outdir <- file.path(root, "analysis", "replication_GSE176078")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

fail_step <- function(step, message) {
  fn <- file.path(outdir, paste0("FAILED_", gsub("[^A-Za-z0-9]+", "_", step), ".txt"))
  writeLines(c(paste0("Step: ", step), paste0("Message: ", message)), fn)
  stop(step, ": ", message, call. = FALSE)
}

need_pkg <- function(pkg, step) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    fail_step(step, paste0("Required R package is not installed or not visible in .libPaths(): ", pkg))
  }
}

for (pkg in c("data.table", "Seurat", "edgeR", "limma", "infercnv")) {
  need_pkg(pkg, "package check")
}

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(edgeR)
  library(limma)
  library(infercnv)
})

set.seed(2026)
cat("project root:", root, "\n")
cat("output dir:", outdir, "\n")

write_audit_rules <- function() {
  rules <- data.table(
    item = c(
      "analysis_order",
      "label_set_1",
      "label_set_2",
      "infercnv_cutoff",
      "infercnv_analysis_mode",
      "infercnv_cluster_by_groups",
      "infercnv_denoise",
      "infercnv_HMM",
      "infercnv_skip_subcluster",
      "infercnv_reference",
      "reference_cap",
      "cnv_score",
      "malignant_threshold",
      "non_malignant_threshold",
      "tier2_counting_unit",
      "tier2_min_cells_per_patient_condition",
      "repeated_measurements",
      "excluded_tiers"
    ),
    value = c(
      "A annotation agreement -> B cell-level/Tier2 -> C G2M program -> D detection diagnostics",
      "Authors' published epithelial labels: Cancer Epithelial=malignant, Normal Epithelial=non_malignant",
      "inferCNV-derived labels using a pre-specified rule before any differential expression",
      as.character(INFERCNV_CUTOFF),
      INFERCNV_ANALYSIS_MODE,
      as.character(INFERCNV_CLUSTER_BY_GROUPS),
      as.character(INFERCNV_DENOISE),
      as.character(INFERCNV_HMM),
      paste0(as.character(SKIP_SUBCLUSTER),
             " (uses inferCNV resume checkpoints to skip STEP 15 sample dendrogram subclustering; ",
             "this is only for tumor_subclusters/heatmap ordering when HMM=FALSE and no_plot=TRUE)"),
      "immune/stromal cells from the same TNBC tumours",
      paste0(REF_CAP_PER_PATIENT_CELLTYPE, " cells per patient x cell type"),
      "mean squared denoised CNV deviation per cell",
      paste0("CNV correlation >= ", CNV_CORR_MALIGNANT_THRESHOLD,
             " and CNV score >= reference median"),
      "CNV correlation < reference 99th percentile and CNV score < reference 99th percentile; otherwise unresolved",
      "patient",
      as.character(MINC),
      "collapse to the measured entity before testing",
      "no pathway enrichment, no DGIdb, no GDSC, no perturbation-signature analysis"
    )
  )
  fwrite(rules, file.path(outdir, "analysis_order_and_thresholds.csv"))
}
write_audit_rules()

candidate_data_dirs <- function(root) {
  unique(c(
    file.path(root, "data", "GSE176078", "Wu_etal_2021_BRCA_scRNASeq"),
    file.path(root, "data", "GSE176078"),
    file.path(root, "data_external", "GSE176078", "Wu_etal_2021_BRCA_scRNASeq"),
    file.path(root, "data_external", "GSE176078")
  ))
}

locate_wu_dir <- function(root) {
  needed <- c("count_matrix_sparse.mtx", "count_matrix_barcodes.tsv",
              "count_matrix_genes.tsv", "metadata.csv")
  dirs <- candidate_data_dirs(root)
  hit <- dirs[vapply(dirs, function(d) all(file.exists(file.path(d, needed))), logical(1))]
  if (length(hit)) return(hit[[1]])
  archive_hits <- unlist(lapply(c(file.path(root, "data", "GSE176078"),
                                  file.path(root, "data_external", "GSE176078"),
                                  file.path(root, "data"),
                                  file.path(root, "data_external")),
                                function(d) if (dir.exists(d)) {
                                  list.files(d, pattern = "GSE176078.*BRCA_scRNASeq.*tar.gz$",
                                             recursive = TRUE, full.names = TRUE)
                                } else character()))
  msg <- paste0(
    "Cannot find extracted GSE176078 scRNA-seq files. Expected files: ",
    paste(needed, collapse = ", "),
    ". Looked in: ", paste(dirs, collapse = " | "), "."
  )
  if (length(archive_hits)) {
    msg <- paste0(msg, " Found archive(s) not extracted: ",
                  paste(archive_hits, collapse = " | "),
                  ". Extract one archive so the four files are directly available.")
  } else {
    msg <- paste0(
      msg,
      " Download from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
    )
  }
  fail_step("data availability verification", msg)
}

data_dir <- locate_wu_dir(root)
cat("GSE176078 data dir:", data_dir, "\n")

read_wu_counts_metadata <- function(data_dir) {
  mtx_file <- file.path(data_dir, "count_matrix_sparse.mtx")
  barcode_file <- file.path(data_dir, "count_matrix_barcodes.tsv")
  gene_file <- file.path(data_dir, "count_matrix_genes.tsv")
  meta_file <- file.path(data_dir, "metadata.csv")
  for (f in c(mtx_file, barcode_file, gene_file, meta_file)) {
    if (!file.exists(f)) fail_step("data availability verification", paste("Missing", f))
  }
  counts <- as(readMM(mtx_file), "CsparseMatrix")
  barcodes <- fread(barcode_file, header = FALSE)$V1
  genes <- fread(gene_file, header = FALSE)
  gene_symbols <- genes[[1]]
  rownames(counts) <- make.unique(as.character(gene_symbols))
  colnames(counts) <- as.character(barcodes)
  meta <- fread(meta_file)
  if (!("cell" %in% names(meta))) {
    first <- names(meta)[1]
    if (first %in% c("V1", "X", "", "...1") || any(as.character(meta[[first]]) %in% colnames(counts))) {
      setnames(meta, first, "cell")
    }
  }
  if (!("cell" %in% names(meta))) {
    rn <- rownames(as.data.frame(meta))
    if (length(rn) == nrow(meta) && any(rn %in% colnames(counts))) {
      meta[, cell := rn]
    } else {
      fail_step("metadata verification", "Cannot identify a cell-barcode column in metadata.csv")
    }
  }
  meta[, cell := as.character(cell)]
  common <- intersect(colnames(counts), meta$cell)
  if (!length(common)) {
    fail_step("metadata verification", "No overlap between count matrix barcodes and metadata cells")
  }
  counts <- counts[, common, drop = FALSE]
  meta <- meta[match(common, cell)]
  if (!identical(meta$cell, colnames(counts))) {
    fail_step("metadata verification", "Metadata rows could not be aligned to count columns")
  }
  list(counts = counts, meta = meta)
}

dat <- read_wu_counts_metadata(data_dir)
counts <- dat$counts
meta <- dat$meta
rm(dat)
cat("loaded counts:", nrow(counts), "genes x", ncol(counts), "cells\n")

column_examples <- function(x, n = 6) {
  v <- unique(as.character(x))
  v <- v[!is.na(v) & v != ""]
  paste(head(v, n), collapse = "; ")
}

metadata_columns <- data.table(
  column = names(meta),
  n_unique = vapply(meta, function(x) length(unique(as.character(x))), integer(1)),
  examples = vapply(meta, column_examples, character(1))
)
metadata_columns[, candidate_role := fifelse(grepl("celltype|cell.type|major|minor|annotation", column, ignore.case = TRUE),
                                             "cell_type_or_annotation",
                                      fifelse(grepl("subtype|state|pam50|clinical|diagnosis|treatment", column, ignore.case = TRUE),
                                             "subtype_or_state",
                                      fifelse(grepl("orig|sample|patient|donor|cid", column, ignore.case = TRUE),
                                             "sample_or_patient", "")))]
fwrite(metadata_columns, file.path(outdir, "metadata_columns_available.csv"))

choose_col <- function(meta, exact = character(), regex = character(), min_unique = 2) {
  for (x in exact) if (x %in% names(meta)) return(x)
  if (length(regex)) {
    hit <- names(meta)[Reduce(`|`, lapply(regex, function(r) grepl(r, names(meta), ignore.case = TRUE)))]
    hit <- hit[vapply(hit, function(h) length(unique(as.character(meta[[h]]))) >= min_unique, logical(1))]
    if (length(hit)) return(hit[[1]])
  }
  NA_character_
}

patient_col <- choose_col(meta, exact = c("orig.ident", "sample", "patient", "patient_id", "Sample", "Patient"),
                          regex = c("orig", "sample", "patient", "donor", "cid"))
if (is.na(patient_col)) {
  cid_cols <- names(meta)[vapply(meta, function(x) any(grepl("^CID[0-9A-Za-z-]+$", as.character(x))), logical(1))]
  if (length(cid_cols)) patient_col <- cid_cols[[1]]
}
if (is.na(patient_col)) fail_step("metadata verification", "Cannot identify patient/sample column")

celltype_col <- choose_col(meta, exact = c("celltype_major", "cell_type_major", "celltype", "cell_type",
                                           "major_celltype", "major_cell_type"),
                           regex = c("celltype", "cell.type", "major"))
if (is.na(celltype_col)) {
  fail_step("metadata verification",
            "Cannot identify authors' epithelial/malignant cell-type labels. Required design depends on these labels.")
}

state_cols <- names(meta)[grepl("subtype|state|pam50|scsubtype|cancer|annotation|celltype|major|minor",
                                names(meta), ignore.case = TRUE)]
fwrite(data.table(role = c("patient_col", "celltype_col", rep("possible_state_or_annotation_col", length(state_cols))),
                  column = c(patient_col, celltype_col, state_cols)),
       file.path(outdir, "metadata_columns_used_or_reported.csv"))

meta[, patient := as.character(get(patient_col))]
meta[, author_celltype := as.character(get(celltype_col))]

clinical_col <- NA_character_
clinical_candidates <- names(meta)[grepl("subtype|clinical|diagnosis|tnbc|er|her2|receptor", names(meta), ignore.case = TRUE)]
for (cc in clinical_candidates) {
  vals <- toupper(as.character(meta[[cc]]))
  if (any(grepl("TNBC|TRIPLE", vals))) {
    clinical_col <- cc
    break
  }
}

if (!is.na(clinical_col)) {
  meta[, clinical_subtype_for_filter := as.character(get(clinical_col))]
  tnbc_cells <- meta$cell[grepl("TNBC|TRIPLE", toupper(meta$clinical_subtype_for_filter))]
  tnbc_filter_source <- paste0("metadata column: ", clinical_col)
} else {
  tnbc_cells <- meta$cell[meta$patient %in% FALLBACK_TNBC_SAMPLES]
  tnbc_filter_source <- paste0("fallback sample list: ", paste(FALLBACK_TNBC_SAMPLES, collapse = ", "))
}
if (!length(tnbc_cells)) {
  fail_step("TNBC sample verification", "No TNBC cells identified from metadata or fallback sample list")
}

meta[, is_tnbc := cell %in% tnbc_cells]
fwrite(data.table(tnbc_filter_source = tnbc_filter_source,
                  n_tnbc_cells = sum(meta$is_tnbc),
                  n_tnbc_patients = uniqueN(meta$patient[meta$is_tnbc])),
       file.path(outdir, "tnbc_filter_source.csv"))

tnbc_patients <- sort(unique(meta$patient[meta$is_tnbc]))
cat("TNBC patients:", paste(tnbc_patients, collapse = ", "), "\n")

is_author_malignant <- grepl("^Cancer Epithelial$|Cancer.*Epithelial|Malignant|Tumou?r",
                             meta$author_celltype, ignore.case = TRUE)
is_author_nonmalig <- grepl("^Normal Epithelial$|Normal.*Epithelial|Non.?malignant.*Epithelial|Epithelial.*Normal",
                            meta$author_celltype, ignore.case = TRUE)
author_epi <- meta$is_tnbc & (is_author_malignant | is_author_nonmalig)
if (!any(author_epi & is_author_malignant)) {
  fail_step("authors malignant label verification",
            "Authors' malignant/cancer epithelial labels are not recoverable from metadata")
}
if (!any(author_epi & is_author_nonmalig)) {
  fail_step("authors non-malignant epithelial label verification",
            "Authors' non-malignant/normal epithelial labels are not recoverable from metadata")
}

meta[, author_label := NA_character_]
meta[author_epi & is_author_malignant, author_label := "malignant"]
meta[author_epi & is_author_nonmalig, author_label := "non_malignant"]

sample_list <- meta[is_tnbc == TRUE, .(
  n_cells_total = .N,
  n_author_epithelial = sum(!is.na(author_label)),
  n_author_malignant = sum(author_label == "malignant", na.rm = TRUE),
  n_author_non_malignant = sum(author_label == "non_malignant", na.rm = TRUE),
  author_celltype_examples = paste(head(unique(author_celltype), 10), collapse = "; ")
), by = patient][order(patient)]
sample_list[, author_malignant_fraction := n_author_malignant / (n_author_malignant + n_author_non_malignant)]
fwrite(sample_list, file.path(outdir, "sample_list_and_author_epithelial_counts.csv"))

# ---- A. Annotation agreement ------------------------------------------------
cat("\nA. annotation agreement\n")
gene_pos_file <- file.path(root, "analysis", "epithelial_only_pseudobulk_20260912", "gencode_v19_gene_pos.txt")
if (!file.exists(gene_pos_file)) {
  fail_step("inferCNV setup", paste("Missing gene order file:", gene_pos_file))
}

epi_cells <- meta$cell[author_epi]
reference_meta <- meta[is_tnbc == TRUE & is.na(author_label)]
reference_meta <- reference_meta[!grepl("Epithelial|Cancer|Tumou?r|Malignant", author_celltype, ignore.case = TRUE)]
if (!nrow(reference_meta)) {
  fail_step("inferCNV reference selection",
            "No immune/stromal reference cells found among TNBC tumours after excluding epithelial cells")
}
reference_meta[, ref_group := paste(patient, author_celltype, sep = "__")]
set.seed(2026)
ref_cells <- reference_meta[, {
  if (.N > REF_CAP_PER_PATIENT_CELLTYPE) .SD[sample(.N, REF_CAP_PER_PATIENT_CELLTYPE)] else .SD
}, by = ref_group]$cell

filter_counts <- data.table(
  filter = c("all_loaded_cells", "tnbc_cells", "authors_epithelial_cells",
             "infercnv_reference_cells_before_cap", "infercnv_reference_cells_after_cap"),
  n_cells = c(ncol(counts), sum(meta$is_tnbc), length(epi_cells), nrow(reference_meta), length(ref_cells))
)
fwrite(filter_counts, file.path(outdir, "filter_cell_counts.csv"))

run_infercnv_labels <- function(counts, meta, epi_cells, ref_cells, outdir, gene_pos_file) {
  use_cells <- c(epi_cells, ref_cells)
  sub_counts <- counts[, use_cells, drop = FALSE]
  ann <- data.frame(
    group = c(paste0("epi_", meta$patient[match(epi_cells, meta$cell)]),
              paste0("ref_", make.names(meta$author_celltype[match(ref_cells, meta$cell)]))),
    row.names = use_cells
  )
  ann_file <- file.path(outdir, "infercnv_annotations_GSE176078.txt")
  write.table(ann, ann_file, sep = "\t", quote = FALSE, col.names = FALSE)
  ref_groups <- unique(ann$group[startsWith(ann$group, "ref_")])
  if (!length(ref_groups)) fail_step("inferCNV setup", "No reference groups were written")
  cat("infercnv version:", as.character(packageVersion("infercnv")), "\n")
  ico <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = sub_counts,
    annotations_file = ann_file,
    delim = "\t",
    gene_order_file = gene_pos_file,
    ref_group_names = ref_groups
  )
  icdir <- file.path(outdir, "infercnv")
  dir.create(icdir, showWarnings = FALSE, recursive = TRUE)
  if (SKIP_SUBCLUSTER && INFERCNV_ANALYSIS_MODE == "samples" && !INFERCNV_HMM) {
    f14 <- list.files(icdir, pattern = "^14_invert_log_transform.*\\.infercnv_obj$",
                      full.names = TRUE)
    if (length(f14) > 1) {
      fail_step("inferCNV skip subcluster",
                paste("Multiple STEP 14 checkpoints found:", paste(basename(f14), collapse = ", ")))
    }
    if (length(f14) == 1) {
      f15 <- file.path(icdir, sub("^14_invert_log_transform", "15_no_subclustering",
                                  basename(f14)))
      if (!file.exists(f15)) {
        cat("SKIP_SUBCLUSTER: writing STEP 15 checkpoint from", basename(f14), "\n")
        ok <- file.copy(f14, f15, overwrite = FALSE)
        if (!ok || !file.exists(f15)) {
          fail_step("inferCNV skip subcluster",
                    paste("Could not create", f15, "from", f14))
        }
      } else {
        cat("SKIP_SUBCLUSTER: STEP 15 checkpoint already exists:", basename(f15), "\n")
      }
    }
  }
  ico <- infercnv::run(
    ico,
    cutoff = INFERCNV_CUTOFF,
    out_dir = icdir,
    cluster_by_groups = INFERCNV_CLUSTER_BY_GROUPS,
    denoise = INFERCNV_DENOISE,
    HMM = INFERCNV_HMM,
    analysis_mode = INFERCNV_ANALYSIS_MODE,
    num_threads = 4,
    no_plot = TRUE,
    output_format = NA
  )
  cnvm <- ico@expr.data - 1
  is_ref <- colnames(cnvm) %in% ref_cells
  cnv_score <- colMeans(cnvm^2)
  ref_score <- cnv_score[is_ref]
  refprof <- cnvm[, is_ref, drop = FALSE]
  epi_prof <- cnvm[, !is_ref, drop = FALSE]
  epi_names <- colnames(epi_prof)
  lab <- data.table(
    cell = epi_names,
    patient = meta$patient[match(epi_names, meta$cell)],
    cnv_score = as.numeric(cnv_score[epi_names]),
    cnv_corr = NA_real_,
    corr_thr99 = NA_real_,
    score_ref_median = as.numeric(stats::median(ref_score, na.rm = TRUE)),
    score_ref_p99 = as.numeric(stats::quantile(ref_score, 0.99, na.rm = TRUE))
  )
  safe_cor <- function(mat, vec) {
    suppressWarnings(as.numeric(cor(mat, vec, use = "pairwise.complete.obs")))
  }
  for (p in unique(lab$patient)) {
    idx <- which(lab$patient == p)
    scores <- lab$cnv_score[idx]
    hi <- idx[scores > lab$score_ref_p99[idx][1]]
    if (length(hi) < 20) {
      hi <- idx[order(-scores)][seq_len(min(length(idx), max(20, floor(length(idx) / 10))))]
    }
    cons <- rowMeans(epi_prof[, hi, drop = FALSE])
    lab$cnv_corr[idx] <- safe_cor(epi_prof[, idx, drop = FALSE], cons)
    ref_corr <- safe_cor(refprof, cons)
    lab$corr_thr99[idx] <- if (any(is.finite(ref_corr))) {
      as.numeric(stats::quantile(ref_corr, 0.99, na.rm = TRUE))
    } else {
      Inf
    }
  }
  lab[, infercnv_label := "unresolved"]
  lab[cnv_corr >= CNV_CORR_MALIGNANT_THRESHOLD & cnv_score >= score_ref_median,
      infercnv_label := "malignant"]
  lab[cnv_corr < corr_thr99 & cnv_score < score_ref_p99,
      infercnv_label := "non_malignant"]
  lab[]
}

infer_labs <- run_infercnv_labels(counts, meta, epi_cells, ref_cells, outdir, gene_pos_file)
fwrite(infer_labs, file.path(outdir, "annotation_infercnv_cell_labels.csv"))

labels <- meta[match(infer_labs$cell, cell), .(cell, patient, author_label)]
labels <- cbind(labels, infer_labs[, .(cnv_score, cnv_corr, corr_thr99, infercnv_label)])
fwrite(labels, file.path(outdir, "annotation_cell_level_label_sets.csv"))

conf_tab <- table(author = labels$author_label, infercnv = labels$infercnv_label)
conf <- as.data.table(as.matrix(conf_tab), keep.rownames = "author")
setcolorder(conf, "author")
fwrite(conf, file.path(outdir, "annotation_confusion_matrix.csv"))

resolved <- labels[author_label %in% c("malignant", "non_malignant") &
                     infercnv_label %in% c("malignant", "non_malignant")]
kappa2 <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  if (n == 0) return(NA_real_)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / (n^2)
  if (isTRUE(all.equal(pe, 1))) return(NA_real_)
  (po - pe) / (1 - pe)
}
kappa <- kappa2(resolved$author_label, resolved$infercnv_label)

ann_by_patient <- labels[, .(
  author_malignant = sum(author_label == "malignant", na.rm = TRUE),
  author_non_malignant = sum(author_label == "non_malignant", na.rm = TRUE),
  author_unresolved = sum(is.na(author_label)),
  infercnv_malignant = sum(infercnv_label == "malignant", na.rm = TRUE),
  infercnv_non_malignant = sum(infercnv_label == "non_malignant", na.rm = TRUE),
  infercnv_unresolved = sum(infercnv_label == "unresolved", na.rm = TRUE)
), by = patient][order(patient)]
ann_by_patient[, author_malignant_fraction := author_malignant / (author_malignant + author_non_malignant)]
ann_by_patient[, infercnv_malignant_fraction :=
                 infercnv_malignant / (infercnv_malignant + infercnv_non_malignant)]
ann_by_patient[, disagreement_pp := 100 * abs(author_malignant_fraction - infercnv_malignant_fraction)]
fwrite(ann_by_patient, file.path(outdir, "annotation_agreement_per_patient.csv"))

ann_summary <- data.table(
  metric = c("n_author_epithelial_cells", "n_resolved_by_both", "cohen_kappa_resolved",
             "max_patient_disagreement_pp", "n_author_unresolved_in_epithelial_set",
             "n_infercnv_unresolved_in_epithelial_set"),
  value = c(nrow(labels), nrow(resolved), kappa,
            max(ann_by_patient$disagreement_pp, na.rm = TRUE),
            sum(is.na(labels$author_label)), sum(labels$infercnv_label == "unresolved"))
)
fwrite(ann_summary, file.path(outdir, "annotation_agreement_summary.csv"))

# ---- B. Unit of replication under both label sets ---------------------------
cat("\nB. unit of replication under both label sets\n")

make_label_table <- function(source) {
  if (source == "authors") {
    z <- meta[author_epi, .(cell, patient, label = author_label)]
  } else if (source == "infercnv") {
    z <- labels[infercnv_label %in% c("malignant", "non_malignant"),
                .(cell, patient, label = infercnv_label)]
  } else {
    stop("unknown source")
  }
  z <- z[!is.na(label)]
  z[, label := factor(label, levels = c("non_malignant", "malignant"))]
  z
}

fit_label_set <- function(label_source, counts, outdir) {
  lab <- make_label_table(label_source)
  if (!nrow(lab)) fail_step(paste("label set", label_source), "No labelled epithelial cells")
  if (!all(c("malignant", "non_malignant") %in% as.character(lab$label))) {
    fail_step(paste("label set", label_source),
              "Both malignant and non_malignant epithelial cells are required")
  }
  lab <- lab[match(intersect(colnames(counts), lab$cell), cell)]
  sub_counts <- counts[, lab$cell, drop = FALSE]
  if (!identical(colnames(sub_counts), lab$cell)) {
    fail_step(paste("label set", label_source), "Counts and labels are not aligned")
  }
  fwrite(lab, file.path(outdir, paste0("labels_", label_source, ".csv")))

  so <- CreateSeuratObject(sub_counts)
  so <- NormalizeData(so, verbose = FALSE)
  Idents(so) <- lab$label
  t1 <- FindMarkers(so, ident.1 = "malignant", ident.2 = "non_malignant",
                    test.use = "wilcox", logfc.threshold = 0, min.pct = 0.01)
  t1$gene <- rownames(t1)
  t1$FDR_BH <- p.adjust(t1$p_val, method = "BH")
  t1 <- as.data.table(t1)
  ord <- order(t1$p_val_adj, -abs(t1$avg_log2FC), na.last = TRUE)
  t1[, tier1_rank := NA_integer_]
  t1$tier1_rank[ord] <- seq_len(nrow(t1))
  setcolorder(t1, c("gene", "tier1_rank"))
  fwrite(t1, file.path(outdir, paste0("B_tier1_FindMarkers_", label_source, ".csv")))

  grp <- paste(lab$patient, lab$label, sep = "__")
  ncell <- table(grp)
  okg <- names(ncell)[ncell >= MINC]
  removed_groups <- data.table(
    patient_condition = names(ncell),
    n_cells = as.integer(ncell),
    kept_for_pseudobulk = names(ncell) %in% okg
  )
  fwrite(removed_groups, file.path(outdir, paste0("B_pseudobulk_group_filter_", label_source, ".csv")))
  if (!length(okg)) fail_step(paste("pseudo-bulk", label_source), "No patient x condition groups pass MINC")
  pb <- sapply(okg, function(g) Matrix::rowSums(sub_counts[, grp == g, drop = FALSE]))
  if (is.null(dim(pb))) pb <- matrix(pb, ncol = 1, dimnames = list(rownames(sub_counts), okg))
  info <- data.table(sample = okg,
                     patient = sub("__.*", "", okg),
                     condition = sub(".*__", "", okg),
                     n_cells = as.integer(ncell[okg]))
  fwrite(info, file.path(outdir, paste0("B_pseudobulk_samples_", label_source, ".csv")))
  paired <- intersect(info$patient[info$condition == "malignant"],
                      info$patient[info$condition == "non_malignant"])
  if (length(paired) < 3) {
    fail_step(paste("pseudo-bulk", label_source),
              paste0("Fewer than 3 paired patients pass MINC=", MINC,
                     ". Paired patients: ", paste(paired, collapse = ", ")))
  }
  sel <- which(info$patient %in% paired)
  inf <- as.data.frame(info[sel])
  inf$condition <- factor(inf$condition, levels = c("non_malignant", "malignant"))
  inf$patient <- factor(inf$patient)
  design <- model.matrix(~ patient + condition, inf)
  y <- DGEList(pb[, sel, drop = FALSE], group = inf$condition)
  keep <- filterByExpr(y, design) | rownames(y) %in% intersect(KEY, rownames(y))
  filter_genes <- data.table(
    n_genes_total = nrow(y),
    n_genes_filterByExpr_or_key = sum(keep),
    n_genes_removed = nrow(y) - sum(keep)
  )
  fwrite(filter_genes, file.path(outdir, paste0("B_filterByExpr_counts_", label_source, ".csv")))
  y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  res <- glmQLFTest(fit, coef = "conditionmalignant")
  t2 <- as.data.table(topTags(res, n = Inf)$table, keep.rownames = "gene")
  fwrite(t2, file.path(outdir, paste0("B_tier2_edgeR_QL_", label_source, ".csv")))

  key_dt <- data.table(gene = KEY)
  key_dt <- merge(key_dt, t1[, .(gene, tier1_rank, avg_log2FC, pct.1, pct.2,
                                p_val, p_val_adj, FDR_BH)],
                  by = "gene", all.x = TRUE)
  key_dt <- merge(key_dt, t2[, .(gene, tier2_logFC = logFC, tier2_PValue = PValue,
                                tier2_FDR = FDR)],
                  by = "gene", all.x = TRUE)
  key_dt[, gene_class := fifelse(gene %in% KEY_MITOTIC, "mitotic", "immunoglobulin")]
  key_dt[, label_source := label_source]
  key_dt[, status := fifelse(is.na(tier1_rank) & !(gene %in% rownames(counts)), "missing_gene",
                      fifelse(is.na(tier2_FDR), "not_tested_or_filtered", "tested"))]
  fwrite(key_dt, file.path(outdir, paste0("B_key_gene_ranks_and_patient_FDR_", label_source, ".csv")))

  list(label_source = label_source, labels = lab, t1 = t1, t2 = t2, key = key_dt,
       y = y, design = design, info = inf, paired = paired)
}

fit_authors <- fit_label_set("authors", counts, outdir)
fit_infercnv <- fit_label_set("infercnv", counts, outdir)

question_b <- rbindlist(list(fit_authors$key, fit_infercnv$key), fill = TRUE)
fwrite(question_b, file.path(outdir, "B_key_gene_comparison_both_label_sets.csv"))

# ---- C. Program versus single gene -----------------------------------------
cat("\nC. program versus single gene\n")

program_test <- function(fit_obj) {
  label_source <- fit_obj$label_source
  v <- voom(fit_obj$y, fit_obj$design)
  g2 <- intersect(G2M, rownames(v))
  if (length(g2) < 10) {
    fail_step(paste("G2M program", label_source),
              paste("Too few G2/M genes in the fitted expression matrix:", length(g2)))
  }
  idx <- ids2indices(list(G2M = G2M), rownames(v))
  cam <- camera(v, idx, fit_obj$design, contrast = ncol(fit_obj$design))
  fr <- fry(v, idx, fit_obj$design, contrast = ncol(fit_obj$design))
  prog <- data.table(
    label_source = label_source,
    test = c("camera", "fry"),
    n_genes = c(cam$NGenes[1], fr$NGenes[1]),
    direction = c(cam$Direction[1], fr$Direction[1]),
    p_value = c(cam$PValue[1], fr$PValue[1])
  )
  fwrite(prog, file.path(outdir, paste0("C_g2m_camera_fry_", label_source, ".csv")))
  Z <- t(scale(t(v$E[g2, , drop = FALSE])))
  sc <- colMeans(Z, na.rm = TRUE)
  inf <- fit_obj$info
  paired <- fit_obj$paired
  dif <- sapply(paired, function(p) {
    m <- which(inf$patient == p & inf$condition == "malignant")
    n <- which(inf$patient == p & inf$condition == "non_malignant")
    if (length(m) != 1 || length(n) != 1) return(NA_real_)
    sc[m] - sc[n]
  })
  tt <- if (sum(is.finite(dif)) > 2) t.test(dif[is.finite(dif)])$p.value else NA_real_
  per_patient <- data.table(
    label_source = label_source,
    patient = names(dif),
    g2m_malignant_minus_non_malignant = as.numeric(dif),
    n_g2m_genes = length(g2),
    paired_t_p = tt,
    positive = as.numeric(dif) > 0
  )
  fwrite(per_patient, file.path(outdir, paste0("C_g2m_per_patient_", label_source, ".csv")))
  data.table(
    label_source = label_source,
    g2m_genes_tested = length(g2),
    g2m_positive_patients = sum(per_patient$positive, na.rm = TRUE),
    g2m_paired_patients = nrow(per_patient),
    g2m_paired_t_p = tt,
    camera_p = cam$PValue[1],
    fry_p = fr$PValue[1],
    mitotic_genes_patient_FDR_lt_0_05 = sum(fit_obj$key$gene_class == "mitotic" &
                                             fit_obj$key$tier2_FDR < 0.05, na.rm = TRUE),
    immunoglobulin_genes_patient_FDR_lt_0_05 = sum(fit_obj$key$gene_class == "immunoglobulin" &
                                                    fit_obj$key$tier2_FDR < 0.05, na.rm = TRUE)
  )
}

program_summary <- rbindlist(list(program_test(fit_authors), program_test(fit_infercnv)))
fwrite(program_summary, file.path(outdir, "C_program_vs_single_gene_summary.csv"))

# ---- D. Detection-rate diagnostic ------------------------------------------
cat("\nD. detection-rate diagnostic\n")

is_uncharacterised <- function(gene) {
  grepl("^AC[0-9]", gene) | grepl("^AL[0-9]", gene) | grepl("^LINC", gene) |
    grepl("-AS[0-9]$", gene)
}

detection_diagnostic <- function(fit_obj) {
  label_source <- fit_obj$label_source
  t1 <- copy(fit_obj$t1)
  t1 <- t1[order(tier1_rank)]
  top <- head(t1$gene, TOP_N_DETECTION)
  all_tnbc_cells <- meta$cell[meta$is_tnbc]
  all_tnbc_cells <- intersect(all_tnbc_cells, colnames(counts))
  dr <- Matrix::rowMeans(counts[top, all_tnbc_cells, drop = FALSE] > 0)
  det <- data.table(
    label_source = label_source,
    gene = top,
    tier1_rank = match(top, t1$gene),
    detection_rate_all_tnbc_cells = as.numeric(dr[top]),
    uncharacterised_locus = is_uncharacterised(top)
  )
  det <- merge(det, t1[, .(gene, avg_log2FC, pct.1, pct.2, p_val_adj, FDR_BH)],
               by = "gene", all.x = TRUE)
  setorder(det, tier1_rank)
  fwrite(det, file.path(outdir, paste0("D_detection_rate_top200_", label_source, ".csv")))
  data.table(
    label_source = label_source,
    n_top_genes = nrow(det),
    median_detection_rate = median(det$detection_rate_all_tnbc_cells, na.rm = TRUE),
    fraction_under_1pct = mean(det$detection_rate_all_tnbc_cells < DETECTION_RATE_LOW_THRESHOLD, na.rm = TRUE),
    fraction_uncharacterised_loci = mean(det$uncharacterised_locus, na.rm = TRUE)
  )
}

detection_summary <- rbindlist(list(detection_diagnostic(fit_authors),
                                    detection_diagnostic(fit_infercnv)))
fwrite(detection_summary, file.path(outdir, "D_detection_rate_summary.csv"))

# ---- Headline summary and report -------------------------------------------
cat("\nWriting headline summary and report\n")

max_disagree <- max(ann_by_patient$disagreement_pp, na.rm = TRUE)

headline_annotation <- ann_by_patient[, .(
  row_type = "annotation_patient",
  patient,
  gene = NA_character_,
  gene_class = NA_character_,
  author_malignant_fraction,
  infercnv_malignant_fraction,
  disagreement_pp,
  cohen_kappa_resolved = kappa,
  authors_patient_FDR = NA_real_,
  infercnv_patient_FDR = NA_real_
)]

key_wide <- merge(
  fit_authors$key[, .(gene, gene_class, authors_patient_FDR = tier2_FDR)],
  fit_infercnv$key[, .(gene, infercnv_patient_FDR = tier2_FDR)],
  by = "gene",
  all = TRUE
)
headline_genes <- key_wide[, .(
  row_type = "key_gene",
  patient = NA_character_,
  gene,
  gene_class,
  author_malignant_fraction = NA_real_,
  infercnv_malignant_fraction = NA_real_,
  disagreement_pp = NA_real_,
  cohen_kappa_resolved = kappa,
  authors_patient_FDR,
  infercnv_patient_FDR
)]

headline <- rbindlist(list(headline_annotation, headline_genes), fill = TRUE)
fwrite(headline, file.path(outdir, "headline_summary.csv"))

leading_gene_sentence <- function() {
  auth_sig <- fit_authors$key[gene_class == "mitotic" & tier2_FDR < 0.05, gene]
  inf_sig <- fit_infercnv$key[gene_class == "mitotic" & tier2_FDR < 0.05, gene]
  ig_auth <- fit_authors$key[gene_class == "immunoglobulin" & tier2_FDR < 0.05, gene]
  ig_inf <- fit_infercnv$key[gene_class == "immunoglobulin" & tier2_FDR < 0.05, gene]
  changed <- !(setequal(auth_sig, inf_sig) && setequal(ig_auth, ig_inf))
  paste0(
    "The leading-gene conclusion ",
    if (changed) "differs" else "does not materially differ",
    " between label sets: authors label set has ",
    length(auth_sig), "/12 mitotic genes and ", length(ig_auth),
    "/4 immunoglobulin genes with patient-level FDR < 0.05; inferCNV label set has ",
    length(inf_sig), "/12 mitotic genes and ", length(ig_inf),
    "/4 immunoglobulin genes with patient-level FDR < 0.05."
  )
}

program_sentence <- function() {
  s <- program_summary
  paste0(
    "The G2/M program result is summarized at the patient level under both label sets: ",
    paste(paste0(s$label_source, " ", s$g2m_positive_patients, "/", s$g2m_paired_patients,
                 " patients positive; camera P=", signif(s$camera_p, 3),
                 "; fry P=", signif(s$fry_p, 3)), collapse = "; "),
    ". Compare this with the single-gene counts in the key-gene tables rather than treating cells as replicates."
  )
}

detection_sentence <- function() {
  paste0(
    "Among the top ", TOP_N_DETECTION, " cell-level genes, the low-detection pattern is: ",
    paste(paste0(detection_summary$label_source, " ",
                 round(100 * detection_summary$fraction_under_1pct, 1),
                 "% under 1% detection and ",
                 round(100 * detection_summary$fraction_uncharacterised_loci, 1),
                 "% uncharacterised loci"), collapse = "; "),
    "."
  )
}

report <- c(
  "# GSE176078 compact replication report",
  "",
  "This analysis is a compact replication of a methodological risk, not a new biological study. It does not run pathway enrichment, drug-gene matching, pharmacogenomic validation, or perturbation-signature analysis.",
  "",
  "## Data and labels",
  paste0("Input directory: `", data_dir, "`."),
  paste0("TNBC filter source: ", tnbc_filter_source, "."),
  paste0("Authors' epithelial label column: `", celltype_col, "`. Patient/sample column: `", patient_col, "`."),
  paste0("Reported annotation/state columns: ", paste(state_cols, collapse = ", "), "."),
  "",
  "## A. Annotation agreement",
  paste0("The two label sets disagree by up to ", round(max_disagree, 1),
         " percentage points across patients."),
  paste0("Cohen's kappa over cells resolved by both label sets is ", signif(kappa, 3), "."),
  paste0("inferCNV leaves ", sum(labels$infercnv_label == "unresolved"),
         " authors-epithelial cells unresolved; authors' malignant/non-malignant epithelial labels leave ",
         sum(is.na(labels$author_label)), " unresolved within the epithelial comparison set."),
  "",
  "## B. Unit of replication",
  leading_gene_sentence(),
  "",
  "## C. Program versus single gene",
  program_sentence(),
  "",
  "## D. Detection-rate diagnostic",
  detection_sentence(),
  "",
  "## Reproducibility",
  "See `analysis_order_and_thresholds.csv` for the pre-specified rules and `sessionInfo.txt` for package versions."
)
writeLines(report, file.path(outdir, "REPORT.md"))

pkg_versions <- data.table(
  package = c("R", "Matrix", "data.table", "Seurat", "edgeR", "limma", "infercnv"),
  version = c(R.version.string,
              as.character(packageVersion("Matrix")),
              as.character(packageVersion("data.table")),
              as.character(packageVersion("Seurat")),
              as.character(packageVersion("edgeR")),
              as.character(packageVersion("limma")),
              as.character(packageVersion("infercnv")))
)
fwrite(pkg_versions, file.path(outdir, "package_versions.csv"))
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nDone. Results written to:", outdir, "\n")
