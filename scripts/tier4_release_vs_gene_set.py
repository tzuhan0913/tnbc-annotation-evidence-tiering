"""Tier 4: separate the two things that changed between the original compound list and the re-run.

The original 27 multi-target compounds came from the 869-gene cell-level list queried against
DGIdb v4.4.  The re-run queries the 830-gene patient-level list against the current release.
Querying the original 869-gene list against the current release, with the same multi-target
criterion (a compound interacting with >= 2 query genes, "NULL" drug names removed), separates
the effect of the database release from the effect of the gene list.

(The re-check of the original 27 compounds - "still in DGIdb" versus "hits a query gene" - is
written by run_tier4_dgidb.R.)

Run from the repository root:  python scripts/tier4_release_vs_gene_set.py
Inputs:
  data_external/dgidb/interactions.tsv   (see docs/data_sources.md; found in the repository root
                                          or a folder above it, or set TNBC_DATA_ROOT)
  results/epithelial_only_pseudobulk_R/legacy_gene_lists_in_new_analysis.csv
                                          the original 869-gene list (list == legacy_broad869);
                                          a fixed input, always read from the shipped results/
  results/tier4_dgidb_reanalysis/query_genes_patient.csv   830 patient-level genes
Output:
  results/tier4_dgidb_reanalysis/dgidb_release_vs_gene_set.csv
Set TNBC_RESULTS=analysis to read and write freshly regenerated tables instead of the shipped ones
(the original gene list is still read from results/).
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, '..'))


def find_data_root(marker=os.path.join('data_external', 'dgidb', 'interactions.tsv')):
    """Folder holding data_external/: TNBC_DATA_ROOT, else the working directory, the
    repository root, or a folder above either (when the repository sits in a project folder)."""
    if os.environ.get('TNBC_DATA_ROOT'):
        return os.environ['TNBC_DATA_ROOT']
    for start in (os.getcwd(), REPO):
        p = os.path.abspath(start)
        for _ in range(4):
            if os.path.exists(os.path.join(p, marker)):
                return p
            p = os.path.dirname(p)
    raise SystemExit('data_external/dgidb/interactions.tsv not found - see docs/data_sources.md, '
                     'or set TNBC_DATA_ROOT')


DATA = find_data_root()
RES = os.environ.get('TNBC_RESULTS', os.path.join(REPO, 'results'))
if not os.path.isabs(RES):
    RES = os.path.join(os.getcwd(), RES)
OUT = os.environ.get('TNBC_TIER4_OUT', os.path.join(RES, 'tier4_dgidb_reanalysis'))
# The original gene lists are a fixed historical input: no script regenerates this file, so it is
# always read from the copy shipped in results/, whatever TNBC_RESULTS points at.
LEGACY = os.environ.get('TNBC_LEGACY', os.path.join(REPO, 'results', 'epithelial_only_pseudobulk_R',
                                                    'legacy_gene_lists_in_new_analysis.csv'))
MIN_GENES = 2
BAD = {'NULL', 'NA', 'NONE', ''}

dg = pd.read_csv(os.path.join(DATA, 'data_external', 'dgidb', 'interactions.tsv'), sep='\t',
                 dtype=str, keep_default_na=False, quoting=3)
dg['drug'] = dg['drug_name'].str.strip().str.upper()
named = dg[~dg['drug'].isin(BAD)]

leg = pd.read_csv(LEGACY)
leg = leg[leg['list'] == 'legacy_broad869'].rename(columns={'gene': 'Gene'})
legacy = set(leg['Gene'].astype(str).str.strip()) - {''}
patient = set(pd.read_csv(os.path.join(RES, 'tier4_dgidb_reanalysis', 'query_genes_patient.csv'))['gene'])


def match(genes):
    e = named[named['gene_name'].isin(genes)][['gene_name', 'drug']].drop_duplicates()
    n = e.groupby('drug')['gene_name'].nunique()
    return e, n


rows = []
for label, genes in [('original 869-gene cell-level list', legacy),
                     ('830-gene patient-level list', patient)]:
    e, n = match(genes)
    rows.append({'gene_list': label, 'n_query_genes': len(genes),
                 'genes_with_named_compound': e['gene_name'].nunique(),
                 'compounds_any': len(n), 'compounds_multi_target': int((n >= MIN_GENES).sum()),
                 'dgidb_release': 'current (data_external/dgidb/interactions.tsv)'})
rows.append({'gene_list': 'overlap of the two lists', 'n_query_genes': len(legacy & patient)})
tab = pd.DataFrame(rows)
for c in ('genes_with_named_compound', 'compounds_any', 'compounds_multi_target'):
    tab[c] = tab[c].astype('Int64')
tab.to_csv(os.path.join(OUT, 'dgidb_release_vs_gene_set.csv'), index=False)
print(tab.to_string(index=False))
