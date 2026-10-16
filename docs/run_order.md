# Run order

Everything is run from the repository root. Scripts expect the primary data under `data/` and
the external resources under `data_external/`, as described in [`data_sources.md`](data_sources.md).

**Where outputs go.** The R scripts locate the root by its `data/GSE161529` folder and write their
tables to `analysis/<step>/`. The tables shipped in `results/<step>/` are the outputs used for the
manuscript, under the same sub-folder and file names, and the descriptions below refer to them
by their `results/` location. The Python scripts read the shipped `results/` tables by default;
the environment variable that points them at freshly regenerated tables instead is given with
each script (for the two Tier 4 Python scripts, `TNBC_RESULTS=analysis`).

## 0. Setup

```
Rscript scripts/00_install_r_packages.R
Rscript scripts/01_r_environment_check.R
powershell -File scripts/00_download_gse161529.ps1
powershell -File scripts/00_download_figshare_objects.ps1
```

## 1. Primary cohort, Tiers 0 to 2

```
Rscript scripts/run_epithelial_only_pseudobulk.R
```

Consumes the GSE161529 matrices and the processed atlas objects.
Writes `results/epithelial_only_pseudobulk_R/`:

- `infercnv_annotations.txt`, `epithelial_labels_infercnv.csv` — the Tier 0 calls
- `label_agreement_infercnv_vs_python.csv`, `label_agreement_by_patient_infercnv.csv`,
  `label_agreement_crosstab.csv` — the two implementations of the same rule, cell by cell
- `tier1_cell_level_findmarkers.csv` — Tier 1
- `tier2a_paired_edgeR_QL_min10.csv`, `tier2a_loso_key_genes_min10.csv` — Tier 2

**Fixed historical input, not regenerated.** `results/epithelial_only_pseudobulk_R/legacy_gene_lists_in_new_analysis.csv`
holds the original analysis's 174- and 869-gene lists (Tables S5 and S6 of the article) with
each gene's detection rate and Tier 1 result under the new labels. No script in this repository
writes it, so it is shipped as an input: the Python scripts always read it from `results/`, even
when `TNBC_RESULTS` points them at regenerated tables (override with `TNBC_LEGACY`). The optional
`legacy` modes of `run_tier3_pathfindr.R` and `run_tier4_dgidb.R` read the same 869-gene list from
a non-public file of the original analysis; the default runs do not need it.

The Python implementation of the same malignancy rule, and the direct edgeR-versus-voom
comparison, live in `results/epithelial_only_pseudobulk_20260912/` and
`results/epithelial_only_pseudobulk_R/edgeR_vs_python_voom_min10.csv`.

## 2. Tier 2 sensitivity

```
Rscript scripts/run_tier2_sensitivity.R
```

Writes the `_min20` counterparts, the unpaired model, and the G2/M results
(`g2m_camera_fry_min10.csv`, `g2m_per_patient_min10.csv`,
`tier2_threshold_sensitivity.csv`) into the same folder.

Earlier pseudo-bulk comparisons from the original analysis are not part of the reanalysis; see
[`analysis_scope.md`](analysis_scope.md) for the reasons.

## 3. Tier 3 — pathways

```
Rscript scripts/run_tier3_pathfindr.R
```

Writes `results/tier3_pathfindr_reanalysis/`, including
`original_pathway_claims_recheck.csv`, which is the sixteen-versus-nine comparison.

## 4. Tier 4 — drug–gene matching

```
Rscript scripts/run_tier4_dgidb.R
```

Writes `results/tier4_dgidb_reanalysis/`, including `original_27_compounds_recheck.csv` and
`candidate_drugs_for_external_validation.csv` (532 molecules).

```
python scripts/tier4_release_vs_gene_set.py
```

Separates the two things that changed between the original 27 compounds and the re-run: the
original 869-gene list against the current DGIdb release returns 332 multi-target compounds,
the 830-gene patient-level list 532 (the two lists share 114 genes). The original list is read
from `results/epithelial_only_pseudobulk_R/legacy_gene_lists_in_new_analysis.csv` (rows
`legacy_broad869`; the same list as Table S6). Writes `dgidb_release_vs_gene_set.csv`. `run_tier4_dgidb.R` itself writes
`original_27_compounds_recheck.csv`, with separate columns for "in DGIdb at all" (checked against
the full interactions table, not only the matched rows) and "number of query genes hit". Drug
names given as the literal string "NULL" (94 interaction records over 91 genes in this query) are
removed before ranking, so the candidate list has 532 rows and the edge table 4,008.

## 5. Tier 5 — external pharmacogenomics

The cell-line panel first:

```
Rscript scripts/check_tnbc_cell_line_panel.R
```

writes `results/tier5_gdsc_reanalysis/tnbc_panel_curation.csv` and both panel files. Then:

```
Rscript scripts/run_tier5_gdsc.R                 # PANEL <- "full"   (24 lines, curated)
Rscript -e 'PANEL <- "strict"; source("scripts/run_tier5_gdsc.R")'   # 22 lines
```

The strict run writes the same file names with a `_strictpanel` suffix. Outputs include
`gdsc_molecule_level_calls.csv` (the collapse of repeated screens),
`gdsc_axis_vs_background_calibration.csv` (the empirical background, not a 0.5 null) and
`tier5_target_axis_vs_drug_list.csv`.

## 5b. Tier 4 symbol-matching sensitivity (after Tier 5)

```
python scripts/tier4_symbol_sensitivity.py
```

Run after Tier 5: the last part reuses the Tier 5 molecule calls. It repeats the Tier 4 count
after converting each query symbol to its current approved HGNC symbol where a unique previous
symbol or alias allows it (the two spreadsheet dates in the original list are restored first);
symbols without a unique match are kept as recorded (151 in the original 869-gene list, 13 in the
830-gene list). The HGNC table is the one distributed with the INDRA Python package v1.24.0
(`data_external/hgnc/`, see `data_sources.md`), not necessarily the current HGNC release.
Result: 546 instead of 532 compounds for the patient-level list and 335 instead of 332 for the
original list. Thirteen of the 14 added compounds reach the two-gene threshold through one
renamed gene, GARS (now GARS1), and one through HIST3H2A (now H2AC25); only entinostat is an
evaluable GDSC molecule, and adding it leaves the drug list below background under every
collapse rule. Writes `dgidb_symbol_sensitivity.csv`, `dgidb_symbol_sensitivity_changes.csv` and
`dgidb_symbol_sensitivity_tier5.csv`.

## 6. Replication cohort

```
Rscript scripts/run_replication_GSE176078.R
Rscript scripts/run_replication_GSE176078_followup.R
```

The first writes the label sets, the agreement statistics
(`annotation_agreement_summary.csv`, `annotation_agreement_per_patient.csv`,
`A2_agreement_balanced_accuracy.csv`) and the per-label-set tiered analysis; it also writes
`analysis_order_and_thresholds.csv`, which records that the malignancy rule and its thresholds
were fixed before any differential expression was run.

The second fits the six models — two label sets × three patient sets (`all`, `intersect`,
`intersect_min20`) — and writes `E_model_summary.csv`, `E_intersection_gene_comparison.csv`,
`E_key_genes_by_model.csv`, `E_g2m_by_model.csv` and `E_VERDICT.txt`.

## 7. Figures

```
python scripts/make_figure_from_R_outputs.py   # Figure 3
python scripts/make_figure_4_volcano.py        # Figure 4
python scripts/make_figures_4_5.py             # Figures 5 and 6
python scripts/make_figure_6.py                # Figure 7
```

The script names keep the numbering the figures had in an earlier draft; the comments above
give the numbering used in the submitted manuscript. The output files are named after the
old numbering too (`figure_R_final`, `fig4_volcano`, `fig4_new`, `fig5_new`, `fig6_draft`),
in that order.

Figures are written to `figures/`. The results root can be overridden with the `TNBC_RESULTS`
environment variable (`TNBC_REP_RESULTS` for the last script) and the output folder with
`TNBC_FIGOUT`.

### Figures not produced by code in this repository

**Figure 1** is a schematic of the tiered workflow, drawn by hand.

**Figure 2** shows the t-SNE of the primary cohort coloured by scGNN cluster and by merged
cell-type label. It visualises the clustering step described in the Methods, which belongs to
the upstream processing of the cohort rather than to the reanalysis: the embedding was never
written out as a result file, and the Seurat/SingleCellExperiment object that held it is not
redistributed here. The figure as published therefore cannot be re-rendered point for point
from this repository. The final cell membership of the historical embedding is also unavailable. Re-running
clustering would produce a new analysis, not verify the cell set in the historical panel.
The recorded original QC stage retained 51,189 of 65,376 cells; it is not a verified final
cell count for Figure 2. See `analysis_scope.md` for the provenance limitation.

No quantitative claim in the manuscript depends on Figure 2. Quantitative reanalysis results come from the Tier 0 to Tier 5 outputs under `results/`.
Reported historical malignant fractions and legacy gene lists are retained as historical
inputs, with their provenance limits identified in the manuscript and `analysis_scope.md`.

### One panel needs an intermediate that is not redistributed

Panel F of `make_figure_from_R_outputs.py` (panel F of Figure 3) compares the detection rate
of the manuscript's two legacy gene lists against the detection rate of every gene
FindMarkers tested. The two legacy curves are drawn from
`legacy_gene_lists_in_new_analysis.csv`, which is included here. The grey reference curve
covering all tested genes needs the detection rate of each gene across all 54,820 cells,
which is a property of the whole dataset rather than a result of this analysis, so it is not
shipped. The script prints a note and omits that one curve.

To draw it, write a two-column CSV (gene, fraction of cells in which the gene is detected)
and point `TNBC_DETECTION` at it. From the Seurat object used in step 2:

```r
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
det <- Matrix::rowMeans(counts > 0)
write.csv(data.frame(gene = names(det), detection_rate = as.numeric(det)),
          "detection_rates.csv", row.names = FALSE)
```
