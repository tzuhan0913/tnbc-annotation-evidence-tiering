# Annotation as the zeroth evidence tier — analysis code

Analysis code and derived result tables for:

> Yang T-H, Liao J-S, Yang C-Y, Liao G-S, Hong Z-J, Pan C-W, Chen C-I, Wang Y-J, Hsu K-F,
> Feng A-C, Yu J-C, Chu C-M, Chang Y-T.
> **Annotation is the zeroth evidence tier: malignant-cell labelling determines which drug
> targets a single-cell atlas nominates in triple-negative breast cancer.**

The manuscript reanalyses public triple-negative breast cancer (TNBC) single-cell RNA-seq data
through a tiered evidence workflow that begins with malignant-cell annotation (Tier 0) and
continues through cell-level discovery (Tier 1), patient-level pseudo-bulk testing (Tier 2),
program-level testing (Tier 3), drug–gene matching (Tier 4) and external pharmacogenomic
cross-reference (Tier 5), with an external replication cohort for the annotation tier.

This repository contains **every script used to produce the numbers, figures and tables in the
manuscript**, together with the derived result tables those scripts wrote. It does not contain
the primary data, which are redistributed by their original sources under their own terms; see
[`docs/data_sources.md`](docs/data_sources.md) for how to obtain each input.

---

## What is here

```
scripts/     analysis code, in run order
results/     the derived tables the scripts wrote (the numbers quoted in the manuscript)
figures/     written by the figure scripts; not tracked
docs/        run order, data sources, environment
```

### scripts/

| Script | Tier | What it does |
|---|---|---|
| `00_install_r_packages.R` | — | installs the R packages used below |
| `01_r_environment_check.R` | — | prints the environment and checks package versions |
| `00_download_gse161529.ps1` | — | downloads the primary cohort (GSE161529) |
| `00_download_figshare_objects.ps1` | — | downloads the processed atlas objects of the source publication |
| `run_epithelial_only_pseudobulk.R` | 0–2 | inferCNV malignancy calls, the independent implementation of the same rule, cell-level Wilcoxon testing, paired patient-level edgeR quasi-likelihood testing, LOSO |
| `run_tier2_sensitivity.R` | 2 | minimum-cells-per-group sensitivity (10 vs 20), unpaired sensitivity model, G2/M camera and fry |
| `run_tier3_pathfindr.R` | 3 | pathfindR enrichment on the cell-level and the patient-level gene lists; recheck of the sixteen original pathway claims |
| `run_tier4_dgidb.R` | 4 | DGIdb drug–gene matching on the patient-level list; recheck of the twenty-seven original compounds |
| `check_tnbc_cell_line_panel.R` | 5 | curates the TNBC cell-line panel from Cellosaurus and DepMap expression, and writes both the curated and the strict panel |
| `run_tier5_gdsc.R` | 5 | GDSC cross-reference: molecule-level collapse of repeated screens, empirical background calibration, target-axis and drug-list tests, under both panels |
| `run_replication_GSE176078.R` | 0–3 | replication cohort: inferCNV under the same pre-specified rule, agreement with the published labels, and the tiered analysis under each label set |
| `run_replication_GSE176078_followup.R` | 0–2 | the intersection analysis: the same patients under both label sets, and the ≥ 20 cells-per-group model |
| `make_figure_from_R_outputs.py` | — | Figures 2 and 3 |
| `make_figures_4_5.py` | — | Figures 4 and 5 |
| `make_figure_6.py` | — | Figure 6 |

The run order and what each step consumes and produces are in
[`docs/run_order.md`](docs/run_order.md).

### results/

The derived tables, arranged by the analysis that wrote them:

| Folder | Contents |
|---|---|
| `epithelial_only_pseudobulk_R/` | Tier 0–2 on the primary cohort: inferCNV labels, the two-implementation agreement, cell-level and patient-level statistics, LOSO, G2/M, threshold sensitivity |
| `epithelial_only_pseudobulk_20260912/` | the independent Python implementation of the same malignancy rule, and its comparison with the R implementation |
| `tier3_pathfindr_reanalysis/` | pathfindR enrichment for both gene lists and the recheck of the original pathway claims |
| `tier4_dgidb_reanalysis/` | DGIdb edges, per-drug summaries, the candidate compound list, and the recheck of the original twenty-seven |
| `tier5_gdsc_reanalysis/` | GDSC per-compound and molecule-level summaries, background calibration, the target-axis and drug-list results under both cell-line panels, and the panel curation |
| `replication_GSE176078/` | the replication cohort: label sets, agreement, the six models (two label sets × three patient sets), the intersection gene comparison, and the pre-specified analysis order |

Every number quoted in the manuscript, in Tables 1–3 and in Additional file 1 (Tables S1–S11)
can be traced to one of these files.

---

## Reproducing the analysis

1. Obtain the inputs as described in [`docs/data_sources.md`](docs/data_sources.md).
2. Install the R packages: `Rscript scripts/00_install_r_packages.R`.
3. Check the environment: `Rscript scripts/01_r_environment_check.R`.
4. Run the scripts in the order given in [`docs/run_order.md`](docs/run_order.md).

The scripts write into the same `results/` folder names used here, so a full re-run can be
compared file by file against the tables in this repository.

Two notes that matter for exact reproduction:

- **inferCNV is stochastic in its HMM step.** Cell-level calls are stable in our hands, but a
  re-run may move a small number of boundary cells between the resolved and unresolved classes.
  The malignancy rule itself, and the thresholds, were fixed before any differential expression
  was run and are recorded in the analysis-order files.
- **Database versions move.** DGIdb, GDSC, Cellosaurus and DepMap are living resources. The
  version or release used for each is recorded alongside the results
  (`results/tier4_dgidb_reanalysis/dgidb_source_version.txt`,
  `results/replication_GSE176078/package_versions.csv` and
  `results/replication_GSE176078/sessionInfo.txt`). Re-running against a newer release will not
  give identical counts; that instability is itself one of the manuscript's findings.

---

## Environment

R and Python versions, and the exact package versions used for the reported results, are in
`results/replication_GSE176078/sessionInfo.txt` and
`results/replication_GSE176078/package_versions.csv`.

Principal dependencies: R with Seurat, SingleR, inferCNV v1.18.1, edgeR, limma, pathfindR,
data.table; Python with NumPy, SciPy, pandas and Matplotlib.

---

## License

Code is released under the MIT License (see [`LICENSE`](LICENSE)).

The derived result tables in `results/` are made available under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). They are derived from public
resources that carry their own terms; `docs/data_sources.md` lists each source and its terms.

## Citation

See [`CITATION.cff`](CITATION.cff). Please cite the manuscript, and, if you use a specific
archived version of this code, its Zenodo DOI.
