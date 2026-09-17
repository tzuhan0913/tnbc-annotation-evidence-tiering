# Run order

Everything is run from the repository root. Scripts expect the primary data under `data/` and
the external resources under `data_external/`, as described in [`data_sources.md`](data_sources.md).

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
- `legacy_gene_lists_in_new_analysis.csv` — where the original gene lists land under the new labels

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
python scripts/make_figure_from_R_outputs.py   # Figures 2 and 3
python scripts/make_figures_4_5.py             # Figures 4 and 5
python scripts/make_figure_6.py                # Figure 6
```

Figures are written to `figures/`. The results root can be overridden with the `TNBC_RESULTS`
environment variable (`TNBC_REP_RESULTS` for Figure 6) and the output folder with `TNBC_FIGOUT`.

Figure 1 is a schematic and is not produced by code.
