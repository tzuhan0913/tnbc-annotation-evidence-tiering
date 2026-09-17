# GSE176078 compact replication report

This analysis is a compact replication of a methodological risk, not a new biological study. It does not run pathway enrichment, drug-gene matching, pharmacogenomic validation, or perturbation-signature analysis.

## Data and labels
Input directory: `C:/Users/tzuha/Documents/投稿文章(JS)/data/GSE176078/Wu_etal_2021_BRCA_scRNASeq`.
TNBC filter source: metadata column: subtype.
Authors' epithelial label column: `celltype_major`. Patient/sample column: `orig.ident`.
Reported annotation/state columns: subtype, celltype_subset, celltype_minor, celltype_major.

## A. Annotation agreement
The two label sets disagree by up to 20 percentage points across patients.
Cohen's kappa over cells resolved by both label sets is 0.974.
inferCNV leaves 971 authors-epithelial cells unresolved; authors' malignant/non-malignant epithelial labels leave 0 unresolved within the epithelial comparison set.

## B. Unit of replication
The leading-gene conclusion differs between label sets: authors label set has 12/12 mitotic genes and 0/4 immunoglobulin genes with patient-level FDR < 0.05; inferCNV label set has 3/12 mitotic genes and 0/4 immunoglobulin genes with patient-level FDR < 0.05.

## C. Program versus single gene
The G2/M program result is summarized at the patient level under both label sets: authors 4/4 patients positive; camera P=1.41e-28; fry P=0.00213; infercnv 4/4 patients positive; camera P=7.51e-30; fry P=0.00124. Compare this with the single-gene counts in the key-gene tables rather than treating cells as replicates.

## D. Detection-rate diagnostic
Among the top 200 cell-level genes, the low-detection pattern is: authors 9.5% under 1% detection and 2% uncharacterised loci; infercnv 8% under 1% detection and 1.5% uncharacterised loci.

## Reproducibility
See `analysis_order_and_thresholds.csv` for the pre-specified rules and `sessionInfo.txt` for package versions.
