"""Step 1: load GEO GSE161529 TNBC matrices, keep official-QC cells, add QC metrics and official clusters."""
import sys, gzip, time, numpy as np, pandas as pd, scipy.io, scipy.sparse as sp
U = '/mnt/user-data/uploads/投稿文章(JS)/'
import os
RAW = os.environ.get('RAW_DIR', U + 'data/GSE161529/raw/')
BC = U + 'data/GSE161529/raw/'
SAMPLES = {  # GEO prefix -> official sample name
 'GSM4909281_TN-MH0126': 'TN_0126', 'GSM4909282_TN-MH0135': 'TN_0135', 'GSM4909283_TN-SH0106': 'TN_0106',
 'GSM4909284_TN-MH0114-T2': 'TN_0114_T2', 'GSM4909285_TN-B1-MH4031': 'TN_B1_4031', 'GSM4909286_TN-B1-MH0131': 'TN_B1_0131',
 'GSM4909287_TN-B1-Tum0554': 'TN_B1_0554', 'GSM4909288_TN-B1-MH0177': 'TN_B1_0177'}
feat = pd.read_csv(U + 'data/GSE161529/GSE161529_features.tsv.gz', sep='\t', header=None)
feat.columns = ['ensembl', 'symbol', 'type'][:feat.shape[1]]
meta = pd.read_csv(U + 'analysis/reanalysis/tnbc_object_inspection/TNBC_metadata.csv', index_col=0)
mito = feat.symbol.str.startswith('MT-').values
todo = sys.argv[1:] or list(SAMPLES)
for pre in todo:
    name = SAMPLES[pre]
    import os
    f = RAW + pre + '-matrix.mtx.gz'
    if not os.path.exists(f):
        print('missing', name); continue
    t = time.time()
    M = scipy.io.mmread(f).tocsc().astype(np.float32)          # genes x cells
    bc = pd.read_csv(BC + pre + '-barcodes.tsv.gz', header=None)[0].values
    cells = np.array([name + '_' + b for b in bc])
    keep = np.isin(cells, meta.index.values)
    M = M[:, keep].T.tocsr()                                   # cells x genes
    cells = cells[keep]
    n_umi = np.asarray(M.sum(1)).ravel(); n_gene = np.diff(M.indptr)
    pmt = np.asarray(M[:, mito].sum(1)).ravel() / np.maximum(n_umi, 1)
    obs = pd.DataFrame({'patient': name, 'n_umi': n_umi, 'n_gene': n_gene, 'pct_mt': pmt * 100,
                        'official_cluster': meta.loc[cells, 'seurat_clusters'].values}, index=cells)
    sp.save_npz(f'/home/claude/cnv/data/{name}.npz', M)
    obs.to_csv(f'/home/claude/cnv/data/{name}_obs.csv')
    print(name, M.shape, 'official cells in meta', (meta.group == name).sum(), f'{time.time()-t:.0f}s')
feat.to_csv('/home/claude/cnv/data/features.csv', index=False)
