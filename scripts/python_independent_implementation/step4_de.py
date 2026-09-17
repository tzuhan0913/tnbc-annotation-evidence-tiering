"""Step 4: harmonised Tier 1 (cell-level Wilcoxon) and Tier 2 (pseudo-bulk voom-limma) on the SAME labels.
Tier 2a: paired design ~ patient + condition, patients with >= MINC cells in both groups.
Tier 2b: unpaired design ~ condition, all pseudo-bulks with >= MINC cells.
LOSO for 2a; key-gene table; G2/M programme score test; overlap with legacy DEG lists."""
import sys, json, numpy as np, pandas as pd, scipy.sparse as sp
from scipy import stats
sys.path.insert(0, '/home/claude/cnv'); from common import *; from limma_py import *
OUT = '/home/claude/cnv/out/'
MINC = int(sys.argv[1]) if len(sys.argv) > 1 else 20
KEY = ['PLK1','AURKA','AURKB','TTK','NEK2','MELK','CDK1','CCNB1','CCNB2','GTSE1','TOP2A','MKI67','IGLC2','IGHG1','IGHA1','IGKC','SPP1']
G2M = "HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 PIMREG FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 JPT1 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA".split()

M, obs = load()
lab = pd.read_csv(OUT + 'epithelial_labels.csv', index_col=0)
lab = lab[lab.label.isin(['malignant', 'non_malignant'])]
idx = {n: i for i, n in enumerate(obs.index)}
rows = np.array([idx[n] for n in lab.index])
Me = M[rows]
# gene symbols: keep first occurrence of each symbol
_, first = np.unique(SYM, return_index=True); first = np.sort(first)
Me = Me[:, first]; sym = SYM[first]
mal = (lab.label.values == 'malignant')
summary = {'MINC': MINC, 'n_malignant_cells': int(mal.sum()), 'n_nonmalignant_cells': int((~mal).sum())}

# ---------------- Tier 1: cell-level Wilcoxon (Seurat FindMarkers-like) ----------------
X = lognorm(Me)
pct1 = np.asarray((Me[mal] > 0).mean(0)).ravel(); pct2 = np.asarray((Me[~mal] > 0).mean(0)).ravel()
m1 = np.asarray(np.expm1(X[mal]).mean(0)).ravel() if False else None
Xe = X.copy(); Xe.data = np.expm1(Xe.data)
mean1 = np.asarray(Xe[mal].mean(0)).ravel(); mean2 = np.asarray(Xe[~mal].mean(0)).ravel()
lfc = np.log2(mean1 + 1) - np.log2(mean2 + 1)
test = np.where((np.maximum(pct1, pct2) >= 0.01))[0]
minpct_ok = np.maximum(pct1, pct2) >= 0.1
pv = np.ones(len(sym))
Xc = X.tocsc()
for s in range(0, len(test), 400):
    cols = test[s:s + 400]
    D = Xc[:, cols].toarray()
    pv[cols] = stats.mannwhitneyu(D[mal], D[~mal], axis=0, alternative='two-sided', method='asymptotic').pvalue
t1 = pd.DataFrame({'gene': sym, 'avg_log2FC': lfc, 'pct.1': pct1, 'pct.2': pct2, 'p_val': pv, 'minpct_0.1': minpct_ok}).iloc[test]
t1['p_val_adj_bonferroni'] = np.minimum(t1.p_val * len(sym), 1)
t1['FDR_BH'] = bh(t1.p_val.values)
t1 = t1.sort_values('p_val')
t1['rank'] = np.arange(1, len(t1) + 1)
t1.to_csv(OUT + f'tier1_cell_level_wilcoxon.csv', index=False)
sig05 = t1[(t1['minpct_0.1']) & (t1.p_val_adj_bonferroni < 0.05) & (t1.avg_log2FC.abs() > 0.5)]; sig10 = t1[(t1['minpct_0.1']) & (t1.p_val_adj_bonferroni < 0.05) & (t1.avg_log2FC.abs() > 1)]
summary.update(tier1_sig_abs05=len(sig05), tier1_sig_abs1=len(sig10))
leg174 = pd.read_csv('/home/claude/cnv/data/legacy_stringent174.csv').Gene.str.strip()
leg869 = pd.read_csv('/home/claude/cnv/data/legacy_broad869.csv').Gene.str.strip()
summary.update(overlap_legacy174_with_new_abs1=int(len(set(leg174) & set(sig10.gene))),
               overlap_legacy869_with_new_abs05=int(len(set(leg869) & set(sig05.gene))))

# ---------------- pseudo-bulk ----------------
pb, info = [], []
for (p, l), g in lab.groupby(['patient', 'label']):
    ii = np.where((lab.patient.values == p) & (lab.label.values == l))[0]
    if len(ii) >= MINC:
        pb.append(np.asarray(Me[ii].sum(0)).ravel()); info.append((p, l, len(ii)))
C = np.array(pb).T.astype(float)
info = pd.DataFrame(info, columns=['patient', 'label', 'n_cells'])
info.to_csv(OUT + f'pseudobulk_samples_min{MINC}.csv', index=False)
paired_pts = sorted(set(info[info.label == 'malignant'].patient) & set(info[info.label == 'non_malignant'].patient))
summary['paired_patients'] = paired_pts
summary['pseudobulk_samples'] = info.to_dict('records')

def run(sel, design):
    Cs = C[:, sel]; inf = info.iloc[sel].reset_index(drop=True)
    ncond_min = min((inf.label == 'malignant').sum(), (inf.label == 'non_malignant').sum())
    keep = filter_by_expr(Cs, ncond_min) | np.isin(sym, KEY)
    Cs = Cs[keep]; gs = sym[keep]
    nf = tmm_factors(Cs)
    cond = (inf.label == 'malignant').astype(float).values
    if design == 'paired':
        pts = sorted(inf.patient.unique())
        Xd = np.column_stack([np.ones(len(inf))] + [(inf.patient == p).astype(float).values for p in pts[1:]] + [cond])
    else:
        Xd = np.column_stack([np.ones(len(inf)), cond])
    E, W = voom(Cs, Xd, nf)
    coef, su, sigma, df = lm_weighted(E, Xd, W)
    j = Xd.shape[1] - 1
    t, p, d0, s0 = ebayes_test(coef, su, sigma, df, j)
    tab = pd.DataFrame({'gene': gs, 'logFC': coef[:, j], 't': t, 'P': p, 'FDR': bh(p), 'AveExpr': E.mean(1), 'passes_filterByExpr': filter_by_expr(Cs, ncond_min)})
    return tab.sort_values('P'), E, inf, gs

res = {}
if len(paired_pts) >= 2:
    selp = np.where(info.patient.isin(paired_pts))[0]
    tp, Ep, infp, gsp = run(selp, 'paired'); tp.to_csv(OUT + f'tier2a_paired_pseudobulk_min{MINC}.csv', index=False); res['paired'] = tp
    summary['tier2a_sig_fdr05_abs05'] = int(((tp.FDR < 0.05) & (tp.logFC.abs() > 0.5)).sum())
    # LOSO
    loso = []
    if len(paired_pts) >= 3:
        for drop in paired_pts:
            s2 = np.where(info.patient.isin([p for p in paired_pts if p != drop]))[0]
            tl, *_ = run(s2, 'paired')
            m = tp.merge(tl, on='gene', suffixes=('', '_loso'))
            k = tl.set_index('gene').reindex(KEY)
            loso.append({'dropped': drop, 'spearman_logFC': stats.spearmanr(m.logFC, m.logFC_loso)[0],
                         'top50_overlap': len(set(tp.gene[:50]) & set(tl.gene[:50])),
                         **{f'{g}_logFC': k.loc[g, 'logFC'] for g in KEY}, **{f'{g}_FDR': k.loc[g, 'FDR'] for g in KEY}})
        pd.DataFrame(loso).to_csv(OUT + f'tier2a_loso_min{MINC}.csv', index=False)
    # G2/M programme score: paired differences
    g2 = [g for g in G2M if g in set(gsp)]
    Z = pd.DataFrame(Ep, index=gsp).loc[g2]
    Z = Z.sub(Z.mean(1), axis=0).div(Z.std(1) + 1e-9, axis=0)
    sc = Z.mean(0).values
    diffs = {p: sc[(infp.patient == p) & (infp.label == 'malignant')][0] - sc[(infp.patient == p) & (infp.label == 'non_malignant')][0] for p in paired_pts}
    d = np.array(list(diffs.values()))
    summary['g2m_paired_diffs'] = {k: round(float(v), 3) for k, v in diffs.items()}
    summary['g2m_paired_t_p'] = float(stats.ttest_1samp(d, 0).pvalue) if len(d) > 2 else None
    summary['g2m_positive'] = f'{int((d > 0).sum())}/{len(d)}'
    # competitive test: moderated t of G2M genes vs rest
    tt = tp.set_index('gene').t
    summary['g2m_competitive_ranksum_p'] = float(stats.mannwhitneyu(tt.reindex(g2).dropna(), tt.drop(g2, errors='ignore'), alternative='greater').pvalue)
tu, *_ = run(np.arange(len(info)), 'unpaired'); tu.to_csv(OUT + f'tier2b_unpaired_pseudobulk_min{MINC}.csv', index=False); res['unpaired'] = tu
summary['tier2b_sig_fdr05_abs05'] = int(((tu.FDR < 0.05) & (tu.logFC.abs() > 0.5)).sum())

# ---------------- key genes ----------------
kt = t1.set_index('gene').reindex(KEY)[['avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj_bonferroni', 'FDR_BH', 'rank']]
for nm, tb in res.items():
    kt = kt.join(tb.set_index('gene')[['logFC', 'FDR']].rename(columns={'logFC': f'{nm}_logFC', 'FDR': f'{nm}_FDR'}), how='left')
kt.to_csv(OUT + f'key_genes_tier1_vs_tier2_min{MINC}.csv')
json.dump(summary, open(OUT + f'summary_min{MINC}.json', 'w'), indent=1, default=str)
pd.set_option('display.width', 250); print(json.dumps(summary, indent=1, default=str)); print(kt.round(4))
