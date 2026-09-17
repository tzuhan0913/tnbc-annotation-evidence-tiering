"""Step 2: inferCNV-style copy-number inference on tumour epithelial cells (official clusters 0/2).
Reference = pooled non-epithelial cells (T/NK, myeloid, fibroblast, endothelial; official clusters 1,3,5,8),
capped at 400 cells per patient x cell type.
Steps (mirroring inferCNV defaults for 10x): median-library normalisation; gene filter mean >= 0.1;
log2(x+1); subtract mean of reference cell-type means; clip +/-3; 101-gene pyramidal smoothing within
chromosome; per-cell median centring; subtract smoothed reference mean; noise filter = values within
mean +/- 1.5 SD of the reference matrix set to 0.
Per cell: CNV score = mean squared denoised CNV; CNV correlation = Pearson r between the cell's smoothed
profile and its patient's consensus tumour profile (mean of epithelial cells above the reference 99th
percentile of CNV score)."""
import sys, numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, '/home/claude/cnv'); from common import *
OUT = '/home/claude/cnv/out/'
rng = np.random.default_rng(2026)
M, obs = load()
obs['celltype'] = obs.official_cluster.map({0: 'epithelial', 2: 'epithelial', 1: 'T_NK', 3: 'myeloid', 5: 'fibroblast',
                                            8: 'endothelial', 4: 'plasma', 6: 'B', 7: 'B_myeloid_mixed'})
pos_i = {n: i for i, n in enumerate(obs.index)}
epi = np.where(obs.celltype.values == 'epithelial')[0]
ref = []
for (p, ct), g in obs[obs.celltype.isin(['T_NK', 'myeloid', 'fibroblast', 'endothelial'])].groupby(['patient', 'celltype']):
    ii = np.array([pos_i[n] for n in g.index])
    ref.extend(rng.choice(ii, min(400, len(ii)), replace=False))
ref = np.array(sorted(ref))
print('epithelial', len(epi), 'reference', len(ref), obs.iloc[ref].celltype.value_counts().to_dict())

# ---- genes in genomic order ----
pos = pd.read_csv('/home/claude/cnv/gencode_v19_gene_pos.txt', sep='\t', header=None, names=['sym', 'chr', 'start', 'end'])
pos = pos[pos.chr.isin([f'chr{i}' for i in range(1, 23)])].drop_duplicates('sym')
sym2i = {}
for i, s in enumerate(SYM): sym2i.setdefault(s, i)
pos = pos[pos.sym.isin(sym2i)].copy()
pos['gi'] = pos.sym.map(sym2i); pos['chrn'] = pos.chr.str[3:].astype(int)
pos = pos.sort_values(['chrn', 'start'])
cells = np.concatenate([epi, ref]); isref = np.r_[np.zeros(len(epi), bool), np.ones(len(ref), bool)]
Mc = M[cells]
lib = np.asarray(Mc.sum(1)).ravel(); med = np.median(lib)
Mg = Mc.tocsc()[:, pos.gi.values].tocsr()
Nm = (sp.diags((med / lib).astype(np.float32)) @ Mg).tocsr()
keep = np.asarray(Nm.mean(0)).ravel() >= 0.1
pos = pos[keep].reset_index(drop=True); Nm = Nm.tocsc()[:, keep].tocsr()
print('genes used', len(pos))
chr_bounds = [np.where(pos.chrn.values == c)[0] for c in range(1, 23)]

# ---- reference mean (average of cell-type means) ----
ref_ct = obs.iloc[ref].celltype.values
def logdense(rows):
    X = Nm[rows].toarray(); np.log2(X + 1, out=X); return X
Xr = logdense(np.where(isref)[0])
ref_mean = np.mean([Xr[ref_ct == ct].mean(0) for ct in np.unique(ref_ct)], axis=0)

W = 101; half = W // 2
kern = np.r_[np.arange(1, half + 2), np.arange(half, 0, -1)].astype(np.float32)
def smooth(X):
    out = np.empty_like(X)
    for b in chr_bounds:
        seg = X[:, b]; n = seg.shape[1]
        # convolution along genes with edge renormalisation
        if n >= W: k = kern
        else:
            h = (n - 1) // 2; k = np.r_[np.arange(1, h + 2), np.arange(h, 0, -1)].astype(np.float32)
        from scipy.signal import fftconvolve
        num = fftconvolve(seg, k[None, :], mode='same', axes=1)
        den = np.convolve(np.ones(n, np.float32), k, mode='same')
        out[:, b] = num / den[None, :]
    return out
def process(X):
    X = X - ref_mean
    np.clip(X, -3, 3, out=X)
    X = smooth(X)
    X -= np.median(X, axis=1, keepdims=True)
    return X
Sr = process(Xr)
ref_s_mean = Sr.mean(0)
Sr -= ref_s_mean
mu, sd = Sr.mean(), Sr.std()
lo, hi = mu - 1.5 * sd, mu + 1.5 * sd
def denoise(S):
    D = S.copy(); D[(D > lo) & (D < hi)] = 0; return D
score = np.empty(len(cells), np.float32)
score[isref] = (denoise(Sr) ** 2).mean(1)
# epithelial cells in chunks; keep smoothed profiles for consensus (downcast to float16 to save memory)
epi_rows = np.where(~isref)[0]
prof = np.empty((len(epi_rows), len(pos)), np.float16)
for s in range(0, len(epi_rows), 2000):
    rr = epi_rows[s:s + 2000]
    S = process(logdense(rr)) - ref_s_mean
    score[rr] = (denoise(S) ** 2).mean(1)
    prof[s:s + 2000] = S.astype(np.float16)
    print('chunk', s, flush=True)
res = obs.iloc[cells][['patient', 'celltype', 'official_cluster', 'n_umi', 'n_gene', 'pct_mt']].copy()
res['is_ref'] = isref; res['cnv_score'] = score
s99 = np.quantile(score[isref], 0.99); s95 = np.quantile(score[isref], 0.95)
# per-patient consensus & correlations
res['cnv_corr'] = np.nan
Srf = Sr  # reference profiles
cons_tab = {}
for p in res.patient.unique():
    e = np.where((~isref) & (res.patient.values == p))[0]
    pe = prof[e - 0].astype(np.float32) if False else prof[np.searchsorted(epi_rows, e)].astype(np.float32)
    hi_cells = score[e] > s99
    base = pe[hi_cells] if hi_cells.sum() >= 20 else pe[np.argsort(-score[e])[:max(20, len(e) // 10)]]
    cons = base.mean(0); cons_tab[p] = cons
    def corr(A, v):
        A = A - A.mean(1, keepdims=True); v = v - v.mean()
        return (A @ v) / (np.linalg.norm(A, axis=1) * np.linalg.norm(v) + 1e-9)
    res.iloc[e, res.columns.get_loc('cnv_corr')] = corr(pe, cons)
    rc = corr(Srf, cons)
    res.loc[res.index[isref & (res.patient.values == p)], 'cnv_corr'] = rc[res.patient.values[isref] == p]
    res.loc[res.index[e], 'corr_thr99'] = np.quantile(rc, 0.99)
    res.loc[res.index[e], 'corr_thr95'] = np.quantile(rc, 0.95)
res['score_thr99'] = s99; res['score_thr95'] = s95
res.to_csv(OUT + 'cnv_cell_scores.csv')
np.save(OUT + 'cnv_consensus.npy', np.array([cons_tab[p] for p in sorted(cons_tab)]))
pos.to_csv(OUT + 'cnv_genes.csv', index=False)
# chromosome-level mean per cell (for heatmap/QC): 22 columns
chrmean = np.stack([prof[:, b].astype(np.float32).mean(1) for b in chr_bounds], 1)
pd.DataFrame(chrmean, index=res.index[~isref], columns=[f'chr{c}' for c in range(1, 23)]).to_csv(OUT + 'cnv_chr_means_epithelial.csv')
# save downsampled profiles for heatmap
sel = rng.choice(len(epi_rows), min(3000, len(epi_rows)), replace=False)
np.savez_compressed(OUT + 'cnv_profiles_sample.npz', prof=prof[sel], cells=res.index[epi_rows[sel]].values,
                    refprof=Sr[rng.choice(len(Sr), min(800, len(Sr)), replace=False)].astype(np.float16))
print('done', 's99', s99, 's95', s95)
# ---- chromosome means for reference cells too (diagnostics) ----
refchr = np.stack([Sr[:, b].mean(1) for b in chr_bounds], 1)
pd.DataFrame(refchr, index=res.index[isref], columns=[f'chr{c}' for c in range(1, 23)]).to_csv(OUT + 'cnv_chr_means_reference.csv')
