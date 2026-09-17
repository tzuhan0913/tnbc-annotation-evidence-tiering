# =============================================================================
# Epithelial-only malignancy calling + harmonised Tier 1 / Tier 2 analysis
# (R / Bioconductor version for the final manuscript numbers)
#
# LABEL_SOURCE controls where the malignant / non-malignant labels come from:
#   "python"   - reuse the labels produced by the Python preview pipeline
#                (analysis/epithelial_only_pseudobulk_20260912/epithelial_labels.csv).
#                Fastest, and the DE numbers below are then produced by edgeR in R.
#   "builtin"  - re-derive the labels here with a self-contained implementation of the
#                inferCNV core steps (no infercnv / JAGS needed).  Use this to confirm
#                the labels independently.
#   "infercnv" - use the infercnv package (needs JAGS installed first:
#                https://sourceforge.net/projects/mcmc-jags/  then
#                install.packages("rjags"); BiocManager::install("infercnv")).
#
# Pipeline (identical in all three cases):
#   1. GEO GSE161529 TNBC raw matrices, restricted to the official post-QC cells
#   2. CNV inference on tumour epithelial cells (official clusters 0/2) against pooled
#      non-epithelial reference cells (clusters 1 T/NK, 3 myeloid, 5 fibroblast,
#      8 endothelial; max 400 cells per patient x cell type):
#      median-library normalisation, gene filter (mean >= 0.1), log2(x+1),
#      subtraction of the mean of the reference cell-type means, clipping at +/-3,
#      101-gene pyramidal smoothing within chromosome, per-cell median centring,
#      subtraction of the smoothed reference mean, and a 1.5-SD reference noise filter
#   3. Per cell: CNV score (mean squared denoised CNV) and CNV correlation with the
#      patient's consensus tumour profile.  Pre-specified rule:
#        malignant     : CNV correlation >= 0.5 and CNV score >= reference median
#        non-malignant : CNV correlation < reference 99th pct and CNV score < ref 99th pct
#        unresolved    : otherwise (excluded from DE)
#   4. Tier 1: Seurat FindMarkers (Wilcoxon), malignant vs non-malignant epithelial cells
#   5. Tier 2: edgeR QL paired pseudo-bulk (~ patient + condition), unpaired sensitivity,
#      leave-one-patient-out, and G2/M camera / fry gene-set tests
#
# Run from the project root ("投稿文章(JS)") in RStudio:
#   setwd("C:/Users/tzuha/Documents/投稿文章(JS)")
#   source("scripts/run_epithelial_only_pseudobulk.R")
# Needs ~6-8 GB free RAM (54,820-cell sparse matrix); Tier 1 FindMarkers takes ~10-30 min.
# =============================================================================
LABEL_SOURCE   <- "infercnv"  # "infercnv" (manuscript numbers) | "builtin" | "python"
INFERCNV_PLOT  <- FALSE       # TRUE also renders the inferCNV heatmap (slow on 34k cells; useful for Figure 2)
MINC           <- 10          # minimum cells per patient x condition pseudo-bulk (10 and 20 reported)
MAX_EPI        <- Inf         # cap epithelial cells per patient BEFORE CNV inference (e.g. 3000) if
                              # infercnv runs out of memory; cells not scored are dropped from Tier 1/2
SKIP_SUBCLUSTER <- TRUE       # inferCNV STEP 15 runs hclust(parallelDist(...)) on every sample group.
                              # With analysis_mode="samples" its only product is the heatmap dendrogram
                              # (@tumor_subclusters), which nothing downstream of denoising uses when
                              # HMM=FALSE and no_plot=TRUE, but on a >10,000-cell group it can take many
                              # hours.  TRUE writes the step-15 checkpoint from the step-14 object so
                              # inferCNV's own resume_mode skips it; the per-cell CNV matrix is identical.
                              # Set FALSE (and INFERCNV_PLOT <- TRUE) if you want the clustered heatmap.

suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(Seurat); library(edgeR); library(limma)
})
set.seed(2026)
root    <- getwd()
## tolerate being sourced from the project root or from scripts/ : find the folder holding data/GSE161529
if (!dir.exists(file.path(root, "data", "GSE161529"))) {
  up <- normalizePath(file.path(root, ".."), winslash = "/", mustWork = FALSE)
  if (dir.exists(file.path(up, "data", "GSE161529"))) root <- up
}
if (!dir.exists(file.path(root, "data", "GSE161529")))
  stop("Cannot find data/GSE161529 - set the working directory to the project root:\n",
       '  setwd("C:/Users/tzuha/Documents/投稿文章(JS)")')
cat("project root:", root, "\n")
rawdir  <- file.path(root, "data", "GSE161529", "raw")
pydir   <- file.path(root, "analysis", "epithelial_only_pseudobulk_20260912")
outdir  <- file.path(root, "analysis", "epithelial_only_pseudobulk_R")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

samples <- c("GSM4909281_TN-MH0126" = "TN_0126", "GSM4909282_TN-MH0135" = "TN_0135",
             "GSM4909283_TN-SH0106" = "TN_0106", "GSM4909284_TN-MH0114-T2" = "TN_0114_T2",
             "GSM4909285_TN-B1-MH4031" = "TN_B1_4031", "GSM4909286_TN-B1-MH0131" = "TN_B1_0131",
             "GSM4909287_TN-B1-Tum0554" = "TN_B1_0554", "GSM4909288_TN-B1-MH0177" = "TN_B1_0177")
KEY <- c("PLK1","AURKA","AURKB","TTK","NEK2","MELK","CDK1","CCNB1","CCNB2","GTSE1","TOP2A","MKI67",
         "IGKC","IGLC2","IGHG1","IGHA1","SPP1")
G2M <- c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2","CKS1B","MKI67",
         "TMPO","CENPF","TACC3","PIMREG","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11",
         "ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","JPT1","HN1","CDC20","TTK","CDC25C","KIF2C",
         "RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR",
         "CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA")

# ---- 1. load matrices (official post-QC cells only) -------------------------
feat <- fread(file.path(root, "data", "GSE161529", "GSE161529_features.tsv.gz"), header = FALSE)
meta <- fread(file.path(root, "analysis", "reanalysis", "tnbc_object_inspection", "TNBC_metadata.csv"))
setnames(meta, 1, "cell")
gsym <- make.unique(feat$V2)

mats <- list()
for (pre in names(samples)) {
  nm <- samples[[pre]]
  m  <- as(readMM(file.path(rawdir, paste0(pre, "-matrix.mtx.gz"))), "CsparseMatrix")
  bc <- fread(file.path(rawdir, paste0(pre, "-barcodes.tsv.gz")), header = FALSE)$V1
  dimnames(m) <- list(gsym, paste0(nm, "_", bc))
  mats[[nm]] <- m[, colnames(m) %in% meta$cell, drop = FALSE]
  cat(nm, ncol(mats[[nm]]), "cells\n")
}
counts <- do.call(cbind, mats); rm(mats); gc()
meta   <- meta[match(colnames(counts), meta$cell)]
ctype  <- c(`0` = "epithelial", `2` = "epithelial", `1` = "T_NK", `3` = "myeloid", `5` = "fibroblast",
            `8` = "endothelial", `4` = "plasma", `6` = "B", `7` = "B_myeloid_mixed")[as.character(meta$seurat_clusters)]
patient <- meta$group
print(table(patient, ctype))

# ---- 2. labels --------------------------------------------------------------
epi     <- which(ctype == "epithelial")
if (is.finite(MAX_EPI)) {
  epi <- unlist(lapply(split(epi, patient[epi]),
                       function(ix) if (length(ix) > MAX_EPI) sample(ix, MAX_EPI) else ix), use.names = FALSE)
  cat("epithelial cells after capping at", MAX_EPI, "per patient:", length(epi), "\n")
}
refpool <- which(ctype %in% c("T_NK", "myeloid", "fibroblast", "endothelial"))
ref     <- unlist(lapply(split(refpool, paste(patient[refpool], ctype[refpool])),
                         function(ix) if (length(ix) > 400) sample(ix, 400) else ix), use.names = FALSE)

## helper: inferCNV core steps, implemented directly (no infercnv / JAGS dependency)
cnv_profiles <- function(counts, epi, ref, gene_pos_file, win = 101) {
  pos <- fread(gene_pos_file, header = FALSE, col.names = c("sym", "chr", "start", "end"))
  pos <- pos[chr %in% paste0("chr", 1:22)][!duplicated(sym)][sym %in% rownames(counts)]
  pos[, chrn := as.integer(sub("chr", "", chr))]
  setorder(pos, chrn, start)
  use <- c(epi, ref)
  lib <- Matrix::colSums(counts[, use, drop = FALSE]); med <- median(lib)
  Nm  <- counts[pos$sym, use, drop = FALSE] %*% Diagonal(x = med / lib)
  keep <- Matrix::rowMeans(Nm) >= 0.1
  pos <- pos[keep]; Nm <- Nm[keep, , drop = FALSE]
  cat("genes used for CNV:", nrow(Nm), "\n")

  ## sparse smoothing operator: 101-gene pyramidal kernel, within chromosome
  half <- win %/% 2
  kern <- c(1:(half + 1), half:1)
  ii <- jj <- xx <- vector("list", 22)
  off <- 0
  for (cc in 1:22) {
    idx <- which(pos$chrn == cc); n <- length(idx); if (n == 0) next
    k <- if (n >= win) kern else { h <- (n - 1) %/% 2; c(seq_len(h + 1), rev(seq_len(h))) }
    kh <- length(k) %/% 2
    I <- J <- X <- numeric(0)
    for (d in seq_along(k)) {           # offset d - kh - 1
      o <- d - kh - 1
      a <- pmax(1, 1 - o):pmin(n, n - o)
      I <- c(I, a); J <- c(J, a + o); X <- c(X, rep(k[d], length(a)))
    }
    S <- sparseMatrix(i = I, j = J, x = X, dims = c(n, n))
    S <- Diagonal(x = 1 / Matrix::rowSums(S)) %*% S   # edge renormalisation (rows sum to 1)
    S <- as(S, "TsparseMatrix")
    ii[[cc]] <- S@i + 1 + off; jj[[cc]] <- S@j + 1 + off; xx[[cc]] <- S@x
    off <- off + n
  }
  SM <- sparseMatrix(i = unlist(ii), j = unlist(jj), x = unlist(xx), dims = c(nrow(Nm), nrow(Nm)))

  isref   <- seq_along(use) > length(epi)
  ref_ct  <- ctype[use[isref]]
  Xr      <- log2(as.matrix(Nm[, isref, drop = FALSE]) + 1)
  ref_mean <- rowMeans(sapply(unique(ref_ct), function(c0) rowMeans(Xr[, ref_ct == c0, drop = FALSE])))
  col_medians <- function(M) apply(M, 2, median)
  proc <- function(X) {                    # X: genes x cells (dense)
    X <- X - ref_mean
    X[X >  3] <-  3; X[X < -3] <- -3
    X <- as.matrix(SM %*% X)               # 101-gene pyramidal smoothing along genes
    sweep(X, 2, col_medians(X), "-")
  }
  Sr <- proc(Xr); ref_s <- rowMeans(Sr); Sr <- Sr - ref_s
  lo <- mean(Sr) - 1.5 * sd(Sr); hi <- mean(Sr) + 1.5 * sd(Sr)
  dn <- function(S) { S[S > lo & S < hi] <- 0; S }
  score <- setNames(numeric(length(use)), colnames(Nm))
  score[isref] <- colMeans(dn(Sr)^2)
  prof <- matrix(0, nrow(Nm), sum(!isref))
  ec <- which(!isref)
  for (s in seq(1, length(ec), by = 2000)) {
    sel <- ec[s:min(s + 1999, length(ec))]
    S <- proc(log2(as.matrix(Nm[, sel, drop = FALSE]) + 1)) - ref_s
    score[sel] <- colMeans(dn(S)^2); prof[, match(sel, ec)] <- S
    cat("  CNV chunk", s, "\n")
  }
  list(score = score, prof = prof, refprof = Sr, cells = colnames(Nm), isref = isref)
}

label_cells <- function(score, prof, refprof, cells, isref, pat) {
  s50 <- quantile(score[isref], .50); s99 <- quantile(score[isref], .99)
  out <- data.frame(cell = cells[!isref], patient = pat[!isref], cnv_score = score[!isref],
                    cnv_corr = NA_real_, corr_thr99 = NA_real_)
  for (p in unique(out$patient)) {
    e  <- which(out$patient == p)
    hi <- if (sum(score[!isref][e] > s99) >= 20) e[score[!isref][e] > s99] else e[order(-score[!isref][e])][seq_len(max(20, length(e) %/% 10))]
    cons <- rowMeans(prof[, hi, drop = FALSE])
    out$cnv_corr[e]   <- as.numeric(cor(prof[, e, drop = FALSE], cons))
    out$corr_thr99[e] <- quantile(as.numeric(cor(refprof, cons)), .99)
  }
  out$label <- "unresolved"
  out$label[out$cnv_corr >= 0.5 & out$cnv_score >= s50] <- "malignant"
  out$label[out$cnv_corr < out$corr_thr99 & out$cnv_score < s99] <- "non_malignant"
  out
}

if (LABEL_SOURCE == "python") {
  lab <- fread(file.path(pydir, "epithelial_labels.csv"))
  setnames(lab, 1, "cell")
  lab <- as.data.frame(lab[, .(cell, patient, cnv_score, cnv_corr, label)])
} else if (LABEL_SOURCE == "builtin") {
  cn  <- cnv_profiles(counts, epi, ref, file.path(pydir, "gencode_v19_gene_pos.txt"))
  lab <- label_cells(cn$score, cn$prof, cn$refprof, cn$cells, cn$isref, patient[c(epi, ref)])
} else {
  library(infercnv)
  ann <- data.frame(group = c(paste0("epi_", patient[epi]), paste0("ref_", ctype[ref])),
                    row.names = colnames(counts)[c(epi, ref)])
  af <- file.path(outdir, "infercnv_annotations.txt")
  write.table(ann, af, sep = "\t", quote = FALSE, col.names = FALSE)
  ico <- CreateInfercnvObject(counts[, c(epi, ref)], annotations_file = af, delim = "\t",
                              gene_order_file = file.path(pydir, "gencode_v19_gene_pos.txt"),
                              ref_group_names = unique(ann$group[startsWith(ann$group, "ref_")]))
  cat("infercnv version:", as.character(packageVersion("infercnv")), "\n")
  icdir <- file.path(outdir, "infercnv")
  if (SKIP_SUBCLUSTER && !INFERCNV_PLOT && dir.exists(icdir)) {
    f14 <- list.files(icdir, pattern = "^14_invert_log_transform.*\\.infercnv_obj$", full.names = TRUE)
    if (length(f14) == 1) {
      f15 <- file.path(icdir, sub("^14_invert_log_transform", "15_no_subclustering", basename(f14)))
      if (!file.exists(f15)) {
        cat("SKIP_SUBCLUSTER: writing the STEP 15 checkpoint from", basename(f14), "\n")
        file.copy(f14, f15)          # partition_method='none' only fills @tumor_subclusters
      }
    }
  }
  ico <- infercnv::run(ico, cutoff = 0.1, out_dir = icdir, cluster_by_groups = TRUE,
                       denoise = TRUE, HMM = FALSE, analysis_mode = "samples", num_threads = 4,
                       no_plot = !INFERCNV_PLOT, output_format = if (INFERCNV_PLOT) "png" else NA)
  cnvm  <- ico@expr.data - 1
  isref <- colnames(cnvm) %in% colnames(counts)[ref]
  lab <- label_cells(colMeans(cnvm^2), cnvm[, !isref, drop = FALSE], cnvm[, isref, drop = FALSE],
                     colnames(cnvm), isref, patient[match(colnames(cnvm), colnames(counts))])
}
fwrite(lab, file.path(outdir, sprintf("epithelial_labels_%s.csv", LABEL_SOURCE)))
print(table(lab$patient, lab$label))

# ---- 3. Tier 1: cell-level FindMarkers on the harmonised labels -------------
keepc <- lab$cell[lab$label != "unresolved"]
so <- CreateSeuratObject(counts[, keepc, drop = FALSE])
so <- NormalizeData(so, verbose = FALSE)
Idents(so) <- factor(lab$label[match(keepc, lab$cell)])
t1 <- FindMarkers(so, ident.1 = "malignant", ident.2 = "non_malignant", test.use = "wilcox",
                  logfc.threshold = 0, min.pct = 0.01)
t1$gene <- rownames(t1); t1$FDR_BH <- p.adjust(t1$p_val, "BH")
fwrite(t1, file.path(outdir, "tier1_cell_level_findmarkers.csv"))
cat("Tier 1 |log2FC|>0.5 & Bonferroni<0.05:", sum(t1$p_val_adj < .05 & abs(t1$avg_log2FC) > .5), "\n")

# ---- 4. Tier 2: edgeR QL pseudo-bulk ---------------------------------------
grp   <- paste(lab$patient, lab$label, sep = "__")[match(keepc, lab$cell)]
ncell <- table(grp); okg <- names(ncell)[ncell >= MINC]
pb    <- sapply(okg, function(g) Matrix::rowSums(counts[, keepc[grp == g], drop = FALSE]))
info  <- data.frame(sample = okg, patient = sub("__.*", "", okg), condition = sub(".*__", "", okg),
                    n_cells = as.integer(ncell[okg]))
fwrite(info, file.path(outdir, sprintf("pseudobulk_samples_min%d.csv", MINC)))
paired <- intersect(info$patient[info$condition == "malignant"], info$patient[info$condition == "non_malignant"])
cat("paired patients:", paste(paired, collapse = ", "), "\n")

fit_qlf <- function(sel, paired_design = TRUE) {
  inf <- info[sel, ]
  inf$condition <- factor(inf$condition, levels = c("non_malignant", "malignant"))
  inf$patient   <- factor(inf$patient)
  design <- if (paired_design) model.matrix(~ patient + condition, inf) else model.matrix(~ condition, inf)
  y    <- DGEList(pb[, sel], group = inf$condition)
  keep <- filterByExpr(y, design) | rownames(y) %in% KEY
  y    <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])    # calcNormFactors() in edgeR < 4.0
  y    <- estimateDisp(y, design)
  fit  <- glmQLFit(y, design, robust = TRUE)
  res  <- glmQLFTest(fit, coef = "conditionmalignant")
  list(tab = topTags(res, n = Inf)$table, y = y, design = design)
}
full <- fit_qlf(which(info$patient %in% paired), TRUE)
fwrite(data.table(gene = rownames(full$tab), full$tab), file.path(outdir, sprintf("tier2a_paired_edgeR_QL_min%d.csv", MINC)))
unp <- fit_qlf(seq_len(nrow(info)), FALSE)
fwrite(data.table(gene = rownames(unp$tab), unp$tab), file.path(outdir, sprintf("tier2b_unpaired_edgeR_QL_min%d.csv", MINC)))
loso <- rbindlist(lapply(paired, function(p) {
  r <- fit_qlf(which(info$patient %in% setdiff(paired, p)), TRUE)$tab
  data.table(dropped = p, gene = KEY, logFC = r[KEY, "logFC"], FDR = r[KEY, "FDR"])
}))
fwrite(loso, file.path(outdir, sprintf("tier2a_loso_key_genes_min%d.csv", MINC)))

# ---- 5. G2/M programme: camera / fry ---------------------------------------
v   <- voom(full$y, full$design)
idx <- ids2indices(list(G2M = G2M), rownames(v))
cam <- camera(v, idx, full$design, contrast = ncol(full$design))
fr  <- fry(v, idx, full$design, contrast = ncol(full$design))
fwrite(data.table(test = c("camera", "fry"), rbind(cam[, c("NGenes","Direction","PValue")],
                                                   fr[, c("NGenes","Direction","PValue")])),
       file.path(outdir, sprintf("g2m_camera_fry_min%d.csv", MINC)))

cat("\n--- key genes, patient level (edgeR QL) ---\n")
print(full$tab[intersect(KEY, rownames(full$tab)), c("logFC", "PValue", "FDR")])
print(cam); print(fr)
sessionInfo()
