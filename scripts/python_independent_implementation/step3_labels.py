"""Step 3: pre-specified malignancy rule + comparison with legacy (Table S3) and official labels.
malignant      : CNV correlation >= 0.5 and CNV score >= reference median
non-malignant  : CNV correlation < reference 99th pct (per-patient consensus) and CNV score < reference 99th pct
unresolved     : all other epithelial cells (excluded from DE)"""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/home/claude/cnv'); from common import *
r = pd.read_csv('/home/claude/cnv/out/cnv_cell_scores.csv', index_col=0)
ref = r[r.is_ref]; e = r[~r.is_ref].copy()
s50, s99 = ref.cnv_score.quantile(.5), ref.cnv_score.quantile(.99)
e['label'] = 'unresolved'
e.loc[(e.cnv_corr >= 0.5) & (e.cnv_score >= s50), 'label'] = 'malignant'
e.loc[(e.cnv_corr < e.corr_thr99) & (e.cnv_score < s99), 'label'] = 'non_malignant'
e.to_csv('/home/claude/cnv/out/epithelial_labels.csv')
tab = e.groupby(['patient', 'label']).size().unstack(fill_value=0)
tab['epithelial_total'] = tab.sum(1)
tab['malignant_pct_of_resolved'] = (100 * tab.malignant / (tab.malignant + tab.non_malignant)).round(1)
legacy = pd.DataFrame({'patient': ['TN_0106','TN_0114_T2','TN_0126','TN_0135','TN_B1_0131','TN_B1_0177','TN_B1_0554','TN_B1_4031'],
    'legacy_epithelial': [37,153,999,8490,2997,400,411,4797], 'legacy_malignant': [11,90,938,6987,228,184,210,2317]}).set_index('patient')
legacy['legacy_malignant_pct'] = (100 * legacy.legacy_malignant / legacy.legacy_epithelial).round(1)
out = tab.join(legacy, how='left')
out.to_csv('/home/claude/cnv/out/label_comparison_by_patient.csv')
pd.set_option('display.width', 220); print(out)
print('totals', e.label.value_counts().to_dict())
# lineage marker validation
M, obs = load(); X = lognorm(M)
idx = {n: i for i, n in enumerate(obs.index)}
rows = np.array([idx[n] for n in e.index])
genes = ['KRT14','KRT5','ACTA2','TP63','KIT','ELF5','ANKRD30A','ESR1','PGR','AR','AZGP1','MKI67','TOP2A','EPCAM','KRT8','KRT18']
gi = gidx(genes); sub = X[rows][:, gi].toarray()
val = pd.DataFrame(sub, columns=[g for g in genes if g in set(SYM)], index=e.index); val['label'] = e.label.values
print((val.groupby('label').mean()).round(2).T)
val.groupby('label').mean().T.to_csv('/home/claude/cnv/out/label_marker_validation.csv')
