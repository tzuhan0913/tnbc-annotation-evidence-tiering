import numpy as np, pandas as pd, scipy.sparse as sp, os
D = '/home/claude/cnv/data/'
feat = pd.read_csv(D + 'features.csv')
SYM = feat.symbol.values
PATIENTS = ['TN_0106','TN_0114_T2','TN_0126','TN_0135','TN_B1_0131','TN_B1_0177','TN_B1_0554','TN_B1_4031']
def load(patients=None):
    Ms, obs = [], []
    for p in (patients or PATIENTS):
        if not os.path.exists(D + p + '.npz'): continue
        Ms.append(sp.load_npz(D + p + '.npz')); obs.append(pd.read_csv(D + p + '_obs.csv', index_col=0))
    return sp.vstack(Ms).tocsr(), pd.concat(obs)
def lognorm(M, target=1e4):
    s = np.asarray(M.sum(1)).ravel()
    X = sp.diags((target / np.maximum(s, 1)).astype(np.float32)) @ M
    X = X.tocsr(); X.data = np.log1p(X.data); return X
def gidx(genes):
    d = {}
    for i, s in enumerate(SYM):
        d.setdefault(s, i)
    return [d[g] for g in genes if g in d]
MARKERS = {
 'epithelial': ['EPCAM','KRT8','KRT18','KRT19','CDH1','KRT7','KRT5','KRT14','KRT17'],
 'immune': ['PTPRC'], 'T_NK': ['CD3D','CD3E','CD2','NKG7','GNLY'],
 'myeloid': ['LYZ','CD68','CD14','C1QA','C1QB','AIF1','TYROBP'],
 'B_plasma': ['MS4A1','CD79A','MZB1','JCHAIN','IGKC'],
 'fibroblast': ['COL1A1','COL1A2','DCN','LUM','PDGFRB','ACTA2'],
 'endothelial': ['PECAM1','VWF','CDH5','PLVAP'],
}
def scores(X):
    out = {}
    for k, g in MARKERS.items():
        ix = gidx(g); out[k] = np.asarray(X[:, ix].mean(1)).ravel()
    return pd.DataFrame(out)
