# =============================================================================
# Tier 2 only: edgeR QL pseudo-bulk at several minimum-cell thresholds.
#
# Reads the malignancy labels that run_epithelial_only_pseudobulk.R already wrote
# (analysis/epithelial_only_pseudobulk_R/epithelial_labels_infercnv.csv), so it does NOT
# re-run inferCNV and does NOT re-run Seurat FindMarkers.  Takes 2-4 minutes.
#
# Produces, for each threshold in MINC_SET:
#   pseudobulk_samples_min<M>.csv          pseudo-bulk sample table
#   tier2a_paired_edgeR_QL_min<M>.csv      paired model, ~ patient + condition
#   tier2b_unpaired_edgeR_QL_min<M>.csv    unpaired sensitivity model
#   tier2a_loso_key_genes_min<M>.csv       leave-one-patient-out, key genes
#   g2m_camera_fry_min<M>.csv              camera / fry on the G2/M set
#   g2m_per_patient_min<M>.csv             per-patient G2/M score difference + paired t
# and one combined summary:
#   tier2_threshold_sensitivity.csv        key genes side by side across thresholds
#
# Run from the project root:
#   setwd("C:/Users/tzuha/Documents/投稿文章(JS)")
#   source("scripts/run_tier2_sensitivity.R")
# =============================================================================
MINC_SET <- c(10, 20)

suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(edgeR); library(limma)
})
set.seed(2026)
root <- getwd()
if (!dir.exists(file.path(root, "data", "GSE161529"))) {
  up <- normalizePath(file.path(root, ".."), winslash = "/", mustWork = FALSE)
  if (dir.exists(file.path(up, "data", "GSE161529"))) root <- up
}
if (!dir.exists(file.path(root, "data", "GSE161529")))
  stop('Cannot find data/GSE161529 - setwd("C:/Users/tzuha/Documents/投稿文章(JS)") first.')
cat("project root:", root, "\n")
rawdir <- file.path(root, "data", "GSE161529", "raw")
outdir <- file.path(root, "analysis", "epithelial_only_pseudobulk_R")
labf   <- file.path(outdir, "epithelial_labels_infercnv.csv")
if (!file.exists(labf)) stop("Missing ", labf, " - run scripts/run_epithelial_only_pseudobulk.R first.")

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

# ---- labels, then only the cells they name ---------------------------------
lab <- fread(labf)
lab <- lab[label %in% c("malignant", "non_malignant")]
cat("labelled epithelial cells:", nrow(lab), "\n"); print(table(lab$patient, lab$label))

feat <- fread(file.path(root, "data", "GSE161529", "GSE161529_features.tsv.gz"), header = FALSE)
gsym <- make.unique(feat$V2)
mats <- list()
for (pre in names(samples)) {
  nm <- samples[[pre]]
  m  <- as(readMM(file.path(rawdir, paste0(pre, "-matrix.mtx.gz"))), "CsparseMatrix")
  bc <- fread(file.path(rawdir, paste0(pre, "-barcodes.tsv.gz")), header = FALSE)$V1
  dimnames(m) <- list(gsym, paste0(nm, "_", bc))
  mats[[nm]] <- m[, colnames(m) %in% lab$cell, drop = FALSE]     # only the labelled cells
  cat(nm, ncol(mats[[nm]]), "cells\n")
}
counts <- do.call(cbind, mats); rm(mats); invisible(gc())
lab <- lab[match(colnames(counts), cell)]
stopifnot(identical(lab$cell, colnames(counts)))

fit_qlf <- function(pb, info, sel, paired_design) {
  inf <- info[sel, ]
  inf$condition <- factor(inf$condition, levels = c("non_malignant", "malignant"))
  inf$patient   <- factor(inf$patient)
  design <- if (paired_design) model.matrix(~ patient + condition, inf) else model.matrix(~ condition, inf)
  y    <- DGEList(pb[, sel, drop = FALSE], group = inf$condition)
  keep <- filterByExpr(y, design) | rownames(y) %in% KEY
  y    <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  y    <- estimateDisp(y, design)
  fit  <- glmQLFit(y, design, robust = TRUE)
  res  <- glmQLFTest(fit, coef = "conditionmalignant")
  list(tab = topTags(res, n = Inf)$table, y = y, design = design, info = inf)
}

allkey <- list()
for (MINC in MINC_SET) {
  cat("\n==================== MINC =", MINC, "====================\n")
  grp   <- paste(lab$patient, lab$label, sep = "__")
  ncell <- table(grp); okg <- names(ncell)[ncell >= MINC]
  pb    <- sapply(okg, function(g) Matrix::rowSums(counts[, grp == g, drop = FALSE]))
  info  <- data.frame(sample = okg, patient = sub("__.*", "", okg), condition = sub(".*__", "", okg),
                      n_cells = as.integer(ncell[okg]))
  fwrite(info, file.path(outdir, sprintf("pseudobulk_samples_min%d.csv", MINC)))
  print(info)
  paired <- intersect(info$patient[info$condition == "malignant"],
                      info$patient[info$condition == "non_malignant"])
  cat("paired patients (", length(paired), "):", paste(paired, collapse = ", "), "\n")
  if (length(paired) < 3) { cat("fewer than 3 paired patients - skipping this threshold\n"); next }

  full <- fit_qlf(pb, info, which(info$patient %in% paired), TRUE)
  fwrite(data.table(gene = rownames(full$tab), full$tab),
         file.path(outdir, sprintf("tier2a_paired_edgeR_QL_min%d.csv", MINC)))
  unp <- fit_qlf(pb, info, seq_len(nrow(info)), FALSE)
  fwrite(data.table(gene = rownames(unp$tab), unp$tab),
         file.path(outdir, sprintf("tier2b_unpaired_edgeR_QL_min%d.csv", MINC)))
  cat("paired: FDR<0.05 & |logFC|>0.5 :", sum(full$tab$FDR < .05 & abs(full$tab$logFC) > .5), "genes\n")

  loso <- rbindlist(lapply(paired, function(p) {
    r <- fit_qlf(pb, info, which(info$patient %in% setdiff(paired, p)), TRUE)$tab
    data.table(dropped = p, gene = KEY, logFC = r[KEY, "logFC"], FDR = r[KEY, "FDR"])
  }))
  fwrite(loso, file.path(outdir, sprintf("tier2a_loso_key_genes_min%d.csv", MINC)))
  losohit <- loso[gene %in% KEY, .(loso_sig = sum(FDR < 0.05, na.rm = TRUE), loso_n = .N), by = gene]

  # ---- G2/M: camera / fry, plus the per-patient paired score ----
  v   <- voom(full$y, full$design)
  idx <- ids2indices(list(G2M = G2M), rownames(v))
  cam <- camera(v, idx, full$design, contrast = ncol(full$design))
  fr  <- fry(v, idx, full$design, contrast = ncol(full$design))
  fwrite(data.table(test = c("camera", "fry"),
                    rbind(cam[, c("NGenes", "Direction", "PValue")], fr[, c("NGenes", "Direction", "PValue")])),
         file.path(outdir, sprintf("g2m_camera_fry_min%d.csv", MINC)))
  g2 <- intersect(G2M, rownames(v))
  Z  <- t(scale(t(v$E[g2, , drop = FALSE])))
  sc <- colMeans(Z, na.rm = TRUE)
  inf <- full$info
  dif <- sapply(paired, function(p)
    sc[which(inf$patient == p & inf$condition == "malignant")] -
    sc[which(inf$patient == p & inf$condition == "non_malignant")])
  tt <- if (length(dif) > 2) t.test(dif)$p.value else NA_real_
  fwrite(data.table(patient = paired, g2m_diff = as.numeric(dif),
                    n_g2m_genes = length(g2), paired_t_p = tt,
                    positive = sprintf("%d/%d", sum(dif > 0), length(dif))),
         file.path(outdir, sprintf("g2m_per_patient_min%d.csv", MINC)))
  cat(sprintf("G2/M: %d/%d patients positive, paired t p = %.4g | camera p = %.3g | fry p = %.3g\n",
              sum(dif > 0), length(dif), tt, cam$PValue[1], fr$PValue[1]))

  k <- full$tab[intersect(KEY, rownames(full$tab)), c("logFC", "PValue", "FDR")]
  cat("\n--- key genes (MINC =", MINC, ") ---\n"); print(k)
  allkey[[as.character(MINC)]] <- data.table(
    MINC = MINC, n_paired_patients = length(paired), gene = rownames(k),
    logFC = k$logFC, FDR = k$FDR)[losohit, on = "gene"]
}

comb <- rbindlist(allkey, fill = TRUE)
wide <- dcast(comb, gene ~ MINC, value.var = c("logFC", "FDR", "loso_sig", "loso_n"))
wide <- wide[order(match(gene, KEY))]
fwrite(wide, file.path(outdir, "tier2_threshold_sensitivity.csv"))
cat("\n===== threshold sensitivity (key genes) =====\n"); print(wide)
cat("\nwritten to", outdir, "\n")
sessionInfo()
