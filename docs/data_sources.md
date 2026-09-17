# Data sources

None of the primary data are redistributed in this repository. Each source below is public and
is redistributed by its own provider under its own terms.

## Primary cohort — GSE161529

Pal B, Chen Y, Vaillant F, et al. *A single-cell RNA expression atlas of normal, preneoplastic
and tumorigenic states in the human breast.* EMBO J. 2021.

- GEO accession: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529>
- Download helper: `scripts/00_download_gse161529.ps1`
- Expected location after download: `data/GSE161529/`

Eight treatment-naive TNBC tumours were used: TN-0106, TN-0114-T2, TN-0126, TN-0135,
TN-B1-0131, TN-B1-0177, TN-B1-0554, TN-B1-4031.

## Processed atlas objects and the original cell-level analysis

Chen Y, Pal B, Visvader JE, Smyth GK. *Differential abundance and expression analysis in
single-cell transcriptomics.* Sci Data. 2022.

- Figshare: <https://doi.org/10.6084/m9.figshare.17058077>
- Descriptor: <https://doi.org/10.1038/s41597-022-01236-2>
- Code: <https://github.com/yunshun/HumanBreast10X>
- Download helper: `scripts/00_download_figshare_objects.ps1`

These objects supply the published cell-level labels that the annotation tier is compared
against.

## Replication cohort — GSE176078

Wu SZ, Al-Eryani G, Roden DL, et al. *A single-cell and spatially resolved atlas of human
breast cancers.* Nat Genet. 2021.

- GEO accession: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078>
- Processed matrix and cell metadata: <https://doi.org/10.5281/zenodo.4739739>
- Expected location after download: `data/GSE176078/Wu_etal_2021_BRCA_scRNASeq/`

The `subtype` column of the supplied metadata selects the TNBC subset; `celltype_major`
(Cancer Epithelial / Normal Epithelial) supplies the published labels.

## Gene position file for inferCNV

GENCODE v19 gene positions, in the three-column format inferCNV expects
(`gencode_v19_gene_pos.txt`). Obtainable from the inferCNV distribution or from
<https://www.gencodegenes.org/>.

## Drug–gene interactions — DGIdb

- <https://www.dgidb.org/> — bulk download `interactions.tsv` (v5).
- The release used, and its download date, are recorded in
  `results/tier4_dgidb_reanalysis/dgidb_source_version.txt`.

Note: the bulk file contains interaction rows whose drug name is the literal string `NULL`.
`run_tier4_dgidb.R` drops these; leaving them in makes `NULL` the highest-degree "drug" in the
network.

## Pharmacogenomics — GDSC

- Genomics of Drug Sensitivity in Cancer, via Cell Model Passports / Sanger:
  <https://cellmodelpassports.sanger.ac.uk/downloads>
- Expected location: `data_external/`

GDSC screens some molecules more than once under different `DRUG_ID` values.
`run_tier5_gdsc.R` collapses repeated screens to the molecule before any test; the three
collapse rules (`meanz`, `any`, `all`) are all reported.

## Cell-line identity and receptor status

- Cellosaurus: <https://www.cellosaurus.org/>
- DepMap, Expression (Short-read) Public 26Q1: <https://depmap.org/portal/>
- Expected location: the DepMap expression CSV anywhere under `data_external/`
  (`check_tnbc_cell_line_panel.R` searches recursively).

`check_tnbc_cell_line_panel.R` calls a line HER2-positive at `ERBB2` log2(TPM+1) ≥ 9 and flags
7–9 as borderline. The commonly quoted threshold of 7 does not separate the classes in this
release: HER2-amplified lines sit at ≥ 9.65 while MCF7, which is HER2-negative, sits at 5.40.
Both the curated (24-line) and the strict (22-line) panels are written, and Tier 5 is reported
under each.

## Gene sets

- G2/M program: the Seurat `cc.genes` G2/M list.
- KEGG and Gene Ontology gene sets as shipped with pathfindR.
