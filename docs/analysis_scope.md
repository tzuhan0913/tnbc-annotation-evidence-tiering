# Analysis scope and provenance limits

## What the reanalysis uses

Tiers 0 to 5 start from the processed atlas object released with GSE161529 (54,820 cells that
passed the atlas authors' quality control), of which 30,353 tumour epithelial cells enter the
pre-specified inferCNV rule. The patient-level result reported for Tier 2 is the paired edgeR
quasi-likelihood model on the same epithelial cells and labels as Tier 1, with at least 10
cells per patient and condition, across six patients
(`results/epithelial_only_pseudobulk_R/tier2a_paired_edgeR_QL_min10.csv`).

## What it does not use, and why

Two earlier patient-level pseudo-bulk results for this cohort exist in the material handed over
from the original analysis. Neither is combined with, or compared against, the Tier 2 result:

- a comparison of the released atlas tumour clusters with all other cells of the same patients,
  whose comparator includes non-epithelial cells and therefore answers a different question;
- a paired malignant-versus-non-malignant table whose cell labels could not be traced to a
  recorded labelling step.

Both were set aside because of what they compare and where their labels come from, not because
of their results.

## Comparisons with the original analysis

The original analysis is used only through what it reported: per-patient malignant fractions,
the 174- and 869-gene cell-level lists, the sixteen pathway claims and the 27 nominated
compounds. Its per-cell labels were not available, so comparisons of malignant fractions are
patient-level summaries over different cell sets, not cell-by-cell relabelling rates.

## Original processing and Figure 2

The original processing script records a per-sample filter that retained 51,189 of 65,376
cells. The final cell membership of the t-SNE embedding shown in Figure 2 could not be
reconstructed; Figure 2 is historical context, and no quantitative result depends on it.

## Gene symbols

Gene lists are matched to DGIdb by the symbols they carry. The original 869-gene list contains two
spreadsheet date conversions (`3-Sep`, `4-Sep` for SEPT3 and SEPT4), kept as recorded.
`scripts/tier4_symbol_sensitivity.py` repeats Tier 4 after converting symbols to current approved
HGNC symbols where the match is unambiguous (other symbols are kept as recorded), using the HGNC
table distributed with INDRA 1.24.0.
