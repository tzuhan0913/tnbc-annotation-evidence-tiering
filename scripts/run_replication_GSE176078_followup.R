# =============================================================================
# GSE176078 replication, follow-up analysis.
#
# The first run (scripts/run_replication_GSE176078.R) found that the two label sets
# agree on 99.7% of cells yet give 675 vs 52 patient-level significant genes.  That
# comparison is confounded: the two label sets do NOT pair the same patients.
#   authors' labels pair  CID4465, CID44971, CID44991, CID4515
#   inferCNV labels pair  CID44971, CID44991, CID4513, CID4515
# So "labels" and "which patients are comparable at all" are mixed together.
#
# This script separates them.  For each label set it fits the same paired edgeR QL
# model on three patient sets:
#   all        every patient with >= MINC cells in both conditions (reproduces run 1)
#   intersect  only patients paired under BOTH label sets   <- the decisive comparison
#   intersect  again at MINC = 20                            <- group-size sensitivity
# If the difference survives "intersect", it is attributable to the labels.
# If it does not, the finding is that annotation decides which patients are
# comparable, which is a different and equally reportable claim.  Either way the
# script prints which of the two it is.
#
#   setwd("C:/Users/tzuha/Documents/投稿文章(JS)")
#   source("scripts/run_replication_GSE176078_followup.R")
#
# Needs: Matrix, data.table, edgeR, limma.  Expect 10-25 min, mostly matrix reading.
# =============================================================================
MINC        <- 10          # primary minimum cells per patient x condition
MINC_STRICT <- 20          # sensitivity
CHUNK       <- 2e7         # non-zero entries read per pass; lower it if memory is tight

suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(edgeR); library(limma)
})
set.seed(2026)

find_root <- function(marker) {
  p <- getwd()
  for (i in 0:4) {
    if (dir.exists(file.path(p, marker))) return(normalizePath(p, winslash = "/", mustWork = FALSE))
    p <- file.path(p, "..")
  }
  stop("Cannot find '", marker, "' from ", getwd(),
       '\n  setwd("C:/Users/tzuha/Documents/投稿文章(JS)") and source the script again.')
}
root   <- find_root(file.path("data", "GSE176078"))
cat("project root:", root, "\n")
ddir   <- file.path(root, "data", "GSE176078", "Wu_etal_2021_BRCA_scRNASeq")
outdir <- file.path(root, "analysis", "replication_GSE176078")
stopifnot(dir.exists(ddir), dir.exists(outdir))

KEY <- c("PLK1","AURKA","AURKB","TTK","NEK2","MELK","CDK1","CCNB1","CCNB2","GTSE1","TOP2A","MKI67",
         "IGKC","IGLC2","IGHG1","IGHA1")
MITOTIC <- KEY[1:12]; IG <- KEY[13:16]
G2M <- c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2","CKS1B","MKI67",
         "TMPO","CENPF","TACC3","PIMREG","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11",
         "ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","JPT1","HN1","CDC20","TTK","CDC25C","KIF2C",
         "RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR",
         "CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA")

# ---- 1. labels written by the first run -------------------------------------
lab_a <- fread(file.path(outdir, "labels_authors.csv"))
lab_i <- fread(file.path(outdir, "labels_infercnv.csv"))
for (d in list(lab_a, lab_i)) stopifnot(all(c("cell", "patient", "label") %in% names(d)))
lab_a <- lab_a[label %in% c("malignant", "non_malignant")]
lab_i <- lab_i[label %in% c("malignant", "non_malignant")]
cat("labels: authors", nrow(lab_a), "cells | inferCNV", nrow(lab_i), "cells\n")
want <- union(lab_a$cell, lab_i$cell)
cat("union of labelled cells:", length(want), "\n")

# ---- 2. per-patient agreement, including balanced accuracy -------------------
# kappa alone is misleading when one class dominates (93% of these cells are malignant),
# so sensitivity, specificity and balanced accuracy are reported alongside it.
mg <- merge(lab_a[, .(cell, patient, author = label)],
            lab_i[, .(cell, infercnv = label)], by = "cell")
agree_one <- function(d) {
  a <- sum(d$author == "malignant"     & d$infercnv == "malignant")
  b <- sum(d$author == "malignant"     & d$infercnv == "non_malignant")
  cc <- sum(d$author == "non_malignant" & d$infercnv == "malignant")
  dd <- sum(d$author == "non_malignant" & d$infercnv == "non_malignant")
  n <- a + b + cc + dd
  po <- (a + dd) / n
  pe <- ((a + b) * (a + cc) + (cc + dd) * (b + dd)) / n^2
  sens <- if ((a + b) > 0) a / (a + b) else NA_real_
  spec <- if ((cc + dd) > 0) dd / (cc + dd) else NA_real_
  data.table(n_resolved_both = n, prevalence_malignant = round((a + b) / n, 4),
             raw_agreement = round(po, 4),
             cohen_kappa = round(if (pe < 1) (po - pe) / (1 - pe) else NA_real_, 4),
             sensitivity_malignant = round(sens, 4), specificity_non_malignant = round(spec, 4),
             balanced_accuracy = round(mean(c(sens, spec), na.rm = TRUE), 4))
}
agr <- rbindlist(c(list(cbind(patient = "ALL", agree_one(mg))),
                   lapply(split(mg, mg$patient), function(d) cbind(patient = d$patient[1], agree_one(d)))))
fwrite(agr, file.path(outdir, "A2_agreement_balanced_accuracy.csv"))
cat("\n===== A2: agreement with balanced accuracy =====\n"); print(agr)

# ---- 3. counts for the labelled cells only ----------------------------------
bc <- fread(file.path(ddir, "count_matrix_barcodes.tsv"), header = FALSE)$V1
gn <- fread(file.path(ddir, "count_matrix_genes.tsv"), header = FALSE)$V1
cat("\nmatrix dims from tsv:", length(gn), "genes x", length(bc), "cells\n")
keep_col <- which(bc %in% want)
stopifnot(length(keep_col) > 0)
newcol <- integer(length(bc)); newcol[keep_col] <- seq_along(keep_col)
cat("columns kept:", length(keep_col), "of", length(bc), "\n")

mtx <- file.path(ddir, "count_matrix_sparse.mtx")
hdr <- readLines(mtx, n = 20)
dl  <- which(!startsWith(hdr, "%"))[1]                 # first non-comment line = dims
dims <- as.numeric(strsplit(trimws(hdr[dl]), "\\s+")[[1]])
cat("mtx header dims:", dims[1], "x", dims[2], "| nnz:", dims[3], "\n")
stopifnot(dims[1] == length(gn), dims[2] == length(bc))

## Read the triplets in chunks and keep only the wanted columns, so peak memory stays
## bounded instead of holding all ~10^8 non-zero entries at once.
ii <- jj <- integer(0); xx <- numeric(0)
skip0 <- dl; done <- 0
repeat {
  part <- fread(mtx, skip = skip0 + done, nrows = CHUNK, header = FALSE,
                colClasses = c("integer", "integer", "numeric"), showProgress = FALSE)
  if (!nrow(part)) break
  sel <- newcol[part$V2] > 0L
  if (any(sel)) {
    ii <- c(ii, part$V1[sel]); jj <- c(jj, newcol[part$V2[sel]]); xx <- c(xx, part$V3[sel])
  }
  done <- done + nrow(part)
  cat("  read", format(done, big.mark = ","), "of", format(dims[3], big.mark = ","),
      "entries; kept", format(length(xx), big.mark = ","), "\n")
  if (nrow(part) < CHUNK) break
}
stopifnot(done == dims[3])
counts <- sparseMatrix(i = ii, j = jj, x = xx,
                       dims = c(length(gn), length(keep_col)),
                       dimnames = list(make.unique(gn), bc[keep_col]))
rm(ii, jj, xx); invisible(gc())
cat("counts:", nrow(counts), "genes x", ncol(counts), "cells\n")

# ---- 4. paired edgeR QL, three patient sets per label set -------------------
pseudobulk <- function(lab, minc) {
  l <- lab[cell %in% colnames(counts)]
  grp <- paste(l$patient, l$label, sep = "__")
  idx <- match(l$cell, colnames(counts))
  ok  <- names(table(grp))[table(grp) >= minc]
  pb  <- sapply(ok, function(g) Matrix::rowSums(counts[, idx[grp == g], drop = FALSE]))
  info <- data.frame(sample = ok, patient = sub("__.*", "", ok), condition = sub(".*__", "", ok),
                     n_cells = as.integer(table(grp)[ok]))
  list(pb = pb, info = info)
}
fit_paired <- function(pb, info, patients) {
  sel <- which(info$patient %in% patients)
  inf <- info[sel, ]
  inf$condition <- factor(inf$condition, levels = c("non_malignant", "malignant"))
  inf$patient   <- factor(inf$patient)
  design <- model.matrix(~ patient + condition, inf)
  y    <- DGEList(pb[, sel, drop = FALSE], group = inf$condition)
  keep <- filterByExpr(y, design) | rownames(y) %in% c(KEY, G2M)
  y    <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  y    <- estimateDisp(y, design)
  res  <- glmQLFTest(glmQLFit(y, design, robust = TRUE), coef = "conditionmalignant")
  list(tab = topTags(res, n = Inf)$table, y = y, design = design, info = inf)
}
paired_of <- function(info, minc) {
  i <- info[info$n_cells >= minc, ]
  intersect(i$patient[i$condition == "malignant"], i$patient[i$condition == "non_malignant"])
}

SETS <- list(authors = lab_a, infercnv = lab_i)
PB   <- lapply(SETS, pseudobulk, minc = MINC)
PB20 <- lapply(SETS, pseudobulk, minc = MINC_STRICT)
pair10 <- lapply(PB,   function(z) paired_of(z$info, MINC))
pair20 <- lapply(PB20, function(z) paired_of(z$info, MINC_STRICT))
inter10 <- Reduce(intersect, pair10)
inter20 <- Reduce(intersect, pair20)
cat("\npaired patients at MINC =", MINC, ":\n")
for (s in names(pair10)) cat("  ", s, ":", paste(pair10[[s]], collapse = ", "), "\n")
cat("  intersection:", paste(inter10, collapse = ", "), "\n")
cat("paired patients at MINC =", MINC_STRICT, ": intersection =",
    if (length(inter20)) paste(inter20, collapse = ", ") else "(none)", "\n")

MODELS <- list()
for (s in names(SETS)) {
  MODELS[[paste(s, "all",       sep = "|")]] <- list(z = PB[[s]],   p = pair10[[s]], minc = MINC)
  if (length(inter10) >= 2)
    MODELS[[paste(s, "intersect", sep = "|")]] <- list(z = PB[[s]], p = inter10,     minc = MINC)
  if (length(inter20) >= 2)
    MODELS[[paste(s, "intersect_min20", sep = "|")]] <- list(z = PB20[[s]], p = inter20, minc = MINC_STRICT)
}

summ <- list(); keytabs <- list(); fits <- list()
for (nm in names(MODELS)) {
  m <- MODELS[[nm]]
  if (length(m$p) < 2) { cat("skipping", nm, "- fewer than 2 paired patients\n"); next }
  f <- fit_paired(m$z$pb, m$z$info, m$p)
  fits[[nm]] <- f
  tab <- f$tab
  sizes <- f$info[order(f$info$patient, f$info$condition), ]
  k <- tab[intersect(KEY, rownames(tab)), c("logFC", "PValue", "FDR"), drop = FALSE]
  keytabs[[nm]] <- data.table(model = nm, gene = rownames(k),
                              gene_class = ifelse(rownames(k) %in% MITOTIC, "mitotic", "immunoglobulin"),
                              logFC = k$logFC, PValue = k$PValue, FDR = k$FDR)
  summ[[nm]] <- data.table(
    model = nm, label_set = sub("\\|.*", "", nm), patient_set = sub(".*\\|", "", nm),
    minc = m$minc, n_paired_patients = length(m$p), patients = paste(sort(m$p), collapse = ";"),
    min_group_cells = min(sizes$n_cells), median_group_cells = median(sizes$n_cells),
    n_genes_tested = nrow(tab), n_FDR05 = sum(tab$FDR < 0.05),
    n_FDR05_lfc05 = sum(tab$FDR < 0.05 & abs(tab$logFC) > 0.5),
    mitotic_FDR05 = sum(k[MITOTIC[MITOTIC %in% rownames(k)], "FDR"] < 0.05, na.rm = TRUE),
    ig_FDR05 = sum(k[IG[IG %in% rownames(k)], "FDR"] < 0.05, na.rm = TRUE))
  fwrite(data.table(gene = rownames(tab), tab), file.path(outdir,
        sprintf("E_tier2_%s.csv", gsub("[|]", "_", nm))))
  fwrite(as.data.table(sizes), file.path(outdir,
        sprintf("E_groups_%s.csv", gsub("[|]", "_", nm))))
}
S <- rbindlist(summ); fwrite(S, file.path(outdir, "E_model_summary.csv"))
K <- rbindlist(keytabs); fwrite(K, file.path(outdir, "E_key_genes_by_model.csv"))
cat("\n===== E: model summary =====\n"); print(S)

# ---- 5. the decisive contrast: same patients, different labels --------------
verdict <- "inconclusive - the intersection model could not be fitted"
if (all(c("authors|intersect", "infercnv|intersect") %in% names(fits))) {
  ta <- fits[["authors|intersect"]]$tab; ti <- fits[["infercnv|intersect"]]$tab
  g  <- intersect(rownames(ta), rownames(ti))
  cmp <- data.table(gene = g, logFC_authors = ta[g, "logFC"], FDR_authors = ta[g, "FDR"],
                    logFC_infercnv = ti[g, "logFC"], FDR_infercnv = ti[g, "FDR"])
  cmp[, `:=`(sig_authors = FDR_authors < 0.05, sig_infercnv = FDR_infercnv < 0.05)]
  cmp[, discordant := sig_authors != sig_infercnv]
  fwrite(cmp[order(-discordant, FDR_authors)],
         file.path(outdir, "E_intersection_gene_comparison.csv"))
  na <- sum(cmp$sig_authors); ni <- sum(cmp$sig_infercnv); nd <- sum(cmp$discordant)
  r  <- suppressWarnings(cor(cmp$logFC_authors, cmp$logFC_infercnv, method = "pearson"))
  sd_ <- cmp[discordant == TRUE]
  cat(sprintf("\n===== same %d patients, two label sets =====\n", length(inter10)))
  cat(sprintf("genes FDR<0.05: authors %d | inferCNV %d | discordant %d\n", na, ni, nd))
  cat(sprintf("logFC correlation between label sets: r = %.3f\n", r))
  if (nrow(sd_))
    cat(sprintf("of the discordant genes, %.1f%% share the sign of logFC; median |dlogFC| = %.2f\n",
                100 * mean(sd_$logFC_authors * sd_$logFC_infercnv > 0),
                median(abs(sd_$logFC_authors - sd_$logFC_infercnv))))
  a_all <- S[model == "authors|all"]$n_FDR05; i_all <- S[model == "infercnv|all"]$n_FDR05
  gap_all <- abs(a_all - i_all); gap_int <- abs(na - ni)
  verdict <- if (gap_int >= 0.5 * gap_all)
    sprintf(paste("LABELS: the gap persists on identical patients (%d vs %d genes, versus %d vs %d",
                  "when each label set chooses its own patients)."), na, ni, a_all, i_all)
  else
    sprintf(paste("PATIENT SET: on identical patients the gap largely closes (%d vs %d genes, versus",
                  "%d vs %d when each label set chooses its own patients). The annotation acts mainly",
                  "by deciding which patients are comparable at all."), na, ni, a_all, i_all)
  cat("\nVERDICT: ", verdict, "\n", sep = "")
}

# ---- 6. G2/M program on the same models ------------------------------------
g2m <- list()
for (nm in names(fits)) {
  f <- fits[[nm]]
  v <- voom(f$y, f$design)
  idx <- ids2indices(list(G2M = G2M), rownames(v))
  cam <- camera(v, idx, f$design, contrast = ncol(f$design))
  fr  <- fry(v, idx, f$design, contrast = ncol(f$design))
  gg <- intersect(G2M, rownames(v))
  Z  <- t(scale(t(v$E[gg, , drop = FALSE]))); sc <- colMeans(Z, na.rm = TRUE)
  pts <- unique(f$info$patient)
  dif <- sapply(pts, function(p)
    sc[which(f$info$patient == p & f$info$condition == "malignant")] -
    sc[which(f$info$patient == p & f$info$condition == "non_malignant")])
  g2m[[nm]] <- data.table(model = nm, n_g2m_genes = length(gg), n_patients = length(pts),
                          positive = sum(dif > 0), camera_p = cam$PValue[1], fry_p = fr$PValue[1],
                          paired_t_p = if (length(dif) > 2) t.test(dif)$p.value else NA_real_)
}
G <- rbindlist(g2m); fwrite(G, file.path(outdir, "E_g2m_by_model.csv"))
cat("\n===== E: G2/M program by model =====\n"); print(G)

writeLines(c("GSE176078 follow-up", paste("run:", format(Sys.time())),
             paste("MINC:", MINC, "| strict:", MINC_STRICT),
             paste("intersection patients:", paste(inter10, collapse = ", ")),
             paste("VERDICT:", verdict)),
           file.path(outdir, "E_VERDICT.txt"))
cat("\nwritten to", outdir, "\n")
sessionInfo()
