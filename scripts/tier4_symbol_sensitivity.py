"""Tier 4 sensitivity: does matching genes to DGIdb by current HGNC symbols change the counts?

The query gene lists carry the symbols of their source data (the atlas object uses pre-2020 HGNC
symbols, and the original 869-gene list contains two spreadsheet date conversions, 3-Sep and
4-Sep for SEPT3 and SEPT4), whereas DGIdb records genes under current approved symbols.  The main
analysis matches symbols as recorded.  This script converts a query symbol to its current
approved HGNC symbol where that is unambiguous (approved symbol kept; otherwise a unique previous
symbol, then a unique alias; spreadsheet dates restored first), keeps every other symbol as
recorded, and repeats the multi-target count.  The HGNC table is the one distributed with the
INDRA Python package v1.24.0, not necessarily the current HGNC release.

Run from the repository root, AFTER Tier 5 (the last part reuses the Tier 5 molecule calls):
  python scripts/tier4_symbol_sensitivity.py
Inputs:
  data_external/dgidb/interactions.tsv
  data_external/hgnc/hgnc_entries_indra-1.24.0.tsv   indra/resources/hgnc_entries.tsv from the
                                                      INDRA 1.24.0 wheel (docs/data_sources.md)
  results/epithelial_only_pseudobulk_R/legacy_gene_lists_in_new_analysis.csv   original 869 genes
                                          (fixed input, always read from the shipped results/)
  results/tier4_dgidb_reanalysis/query_genes_patient.csv                       830 patient-level genes
  results/tier5_gdsc_reanalysis/gdsc_molecule_level_calls.csv, tier5_target_axis_vs_drug_list.csv
Outputs (results/tier4_dgidb_reanalysis/):
  dgidb_symbol_sensitivity.csv          counts as recorded and after symbol mapping
  dgidb_symbol_sensitivity_changes.csv  every query symbol that the mapping changed
  dgidb_symbol_sensitivity_tier5.csv    the Tier 5 drug-list test with the added compounds
Set TNBC_RESULTS=analysis to read and write freshly regenerated tables instead of the shipped ones
(the original gene list is still read from results/).
"""
import os
import re
from collections import defaultdict

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
DATE = re.compile(r'^(\d{1,2})-(Sep|Mar)$')
DATE_PREFIX = {'Sep': 'SEPT', 'Mar': 'MARCH'}

hg = pd.read_csv(os.path.join(DATA, 'data_external', 'hgnc', 'hgnc_entries_indra-1.24.0.tsv'),
                 sep='\t', dtype=str, keep_default_na=False)
hg = hg[hg['Status'] == 'Approved']
approved = set(hg['Approved symbol'])
prev, alias = defaultdict(set), defaultdict(set)
for sym, al, pv in zip(hg['Approved symbol'], hg['Alias symbols'], hg['Previous symbols']):
    for s in filter(None, (x.strip() for x in pv.split(','))):
        prev[s].add(sym)
    for s in filter(None, (x.strip() for x in al.split(','))):
        alias[s].add(sym)


def to_current(sym):
    m = DATE.match(sym)
    if m:
        sym = DATE_PREFIX[m.group(2)] + m.group(1)
    if sym in approved:
        return sym, 'approved'
    if len(prev.get(sym, ())) == 1:
        return next(iter(prev[sym])), 'previous symbol'
    if len(alias.get(sym, ())) == 1:
        return next(iter(alias[sym])), 'alias'
    return sym, 'not mapped'


dg = pd.read_csv(os.path.join(DATA, 'data_external', 'dgidb', 'interactions.tsv'), sep='\t',
                 dtype=str, keep_default_na=False, quoting=3)
dg['drug'] = dg['drug_name'].str.strip().str.upper()
named = dg[~dg['drug'].isin(BAD)]


def count(genes):
    e = named[named['gene_name'].isin(genes)][['gene_name', 'drug']].drop_duplicates()
    n = e.groupby('drug')['gene_name'].nunique()
    return e['gene_name'].nunique(), len(n), int((n >= MIN_GENES).sum())


leg = pd.read_csv(LEGACY)
leg = leg[leg['list'] == 'legacy_broad869'].rename(columns={'gene': 'Gene'})
lists = {'original 869-gene cell-level list': sorted(set(leg['Gene'].astype(str).str.strip()) - {''}),
         '830-gene patient-level list':
             sorted(set(pd.read_csv(os.path.join(RES, 'tier4_dgidb_reanalysis', 'query_genes_patient.csv'))['gene']))}

rows, changes = [], []
for name, genes in lists.items():
    mapped = {}
    for g in genes:
        cur, how = to_current(g)
        mapped[g] = cur
        if cur != g:
            changes.append({'gene_list': name, 'recorded_symbol': g, 'current_symbol': cur,
                            'matched_by': how, 'in_dgidb_after_mapping': cur in set(named['gene_name'])})
    for label, gs in (('as recorded', set(genes)), ('mapped to current HGNC symbols', set(mapped.values()))):
        g_n, c_any, c_multi = count(gs)
        rows.append({'gene_list': name, 'symbols': label, 'n_query_genes': len(gs),
                     'genes_with_named_compound': g_n, 'compounds_any': c_any,
                     'compounds_multi_target': c_multi})
tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(OUT, 'dgidb_symbol_sensitivity.csv'), index=False)
pd.DataFrame(changes).to_csv(os.path.join(OUT, 'dgidb_symbol_sensitivity_changes.csv'), index=False)
print(tab.to_string(index=False))
ch = pd.DataFrame(changes)
print(len(ch), 'symbols changed;', int(ch['in_dgidb_after_mapping'].sum()), 'of them now found in DGIdb')
print(ch[ch['in_dgidb_after_mapping']].to_string(index=False))

# ---- what the extra compounds would do to the Tier 5 drug-list test ------------------------
# Compounds that become multi-target only after symbol mapping and that are evaluable GDSC
# molecules move from the background into the drug list.  The molecule-level calls written by
# run_tier5_gdsc.R (full panel) are reused; nothing is refitted.
from math import comb


def fisher_two_sided(a, b, c, d):
    r1, c1, n = a + b, a + c, a + b + c + d
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)
    p = lambda x: comb(r1, x) * comb(n - r1, c1 - x) / comb(n, c1)
    p0 = p(a)
    return min(1.0, sum(p(x) for x in range(lo, hi + 1) if p(x) <= p0 * (1 + 1e-7)))


T5 = os.path.join(RES, 'tier5_gdsc_reanalysis')
key = lambda s: re.sub(r'[^A-Z0-9]', '', str(s).upper())
pat = lists['830-gene patient-level list']


def multi(gs):
    e = named[named['gene_name'].isin(gs)][['gene_name', 'drug']].drop_duplicates()
    n = e.groupby('drug')['gene_name'].nunique()
    return set(n[n >= MIN_GENES].index)


added = sorted(multi({to_current(g)[0] for g in pat}) - multi(set(pat)))
mol = pd.read_csv(os.path.join(T5, 'gdsc_molecule_level_calls.csv'))
mol['k'] = mol['molkey'].map(key)
hit = mol[mol['k'].isin({key(d) for d in added})]
tvd = pd.read_csv(os.path.join(T5, 'tier5_target_axis_vs_drug_list.csv'))
out = []
for rule in ('meanz', 'any', 'all'):
    r = tvd[(tvd['set'] == 'Tier 4 drug list') & (tvd['collapse_rule'] == rule)].iloc[0]
    k, n = int(r['favourable']), int(r['n_molecules'])
    bk, bn = int(r['background_favourable']), int(r['background_n'])
    f = int(hit[f'fav_{rule}'].astype(str).str.upper().eq('TRUE').sum())
    k2, n2, bk2, bn2 = k + f, n + len(hit), bk - f, bn - len(hit)
    out.append({'collapse_rule': rule, 'added_compounds': len(added), 'added_evaluable_in_gdsc': len(hit),
                'fraction_as_recorded': round(k / n, 4), 'background_as_recorded': round(bk / bn, 4),
                'fisher_as_recorded': round(fisher_two_sided(k, n - k, bk, bn - bk), 4),
                'fraction_mapped': round(k2 / n2, 4), 'background_mapped': round(bk2 / bn2, 4),
                'fisher_mapped': round(fisher_two_sided(k2, n2 - k2, bk2, bn2 - bk2), 4)})
t5 = pd.DataFrame(out)
t5.to_csv(os.path.join(OUT, 'dgidb_symbol_sensitivity_tier5.csv'), index=False)
print('added compounds:', ', '.join(added))
print('evaluable in GDSC:', ', '.join(hit['molkey']))
print(t5.to_string(index=False))
