"""Final figure, built from the official inferCNV + edgeR QL outputs (R)."""
import json, os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
_HERE = os.path.dirname(os.path.abspath(__file__))
R = os.environ.get('TNBC_RESULTS',
     os.path.join(_HERE, '..', 'results', 'epithelial_only_pseudobulk_R')) + os.sep
OUT = os.environ.get('TNBC_FIGOUT', os.path.join(_HERE, '..', 'figures')) + os.sep
os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, GREEN, GRAY, INK, MUTED = '#2E6BA8', '#D9822B', '#4C9A6B', '#B5B5B5', '#1a1a1a', '#6b6b6b'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.edgecolor': '#999999',
                     'axes.linewidth': .6, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.labelcolor': INK, 'text.color': INK, 'figure.facecolor': 'white'})
nice = lambda p: p.replace('TN_', 'TN-').replace('_', '-')
KEY_M = ['PLK1', 'AURKA', 'AURKB', 'TTK', 'NEK2', 'MELK', 'CDK1', 'CCNB1', 'CCNB2', 'GTSE1', 'TOP2A', 'MKI67']
KEY_I = ['IGKC', 'IGLC2', 'IGHG1', 'IGHA1']

lab = pd.read_csv(R + 'epithelial_labels_infercnv.csv')
t1 = pd.read_csv(R + 'tier1_cell_level_findmarkers.csv').sort_values('p_val').reset_index(drop=True)
t1['rank'] = np.arange(1, len(t1) + 1)
t2 = pd.read_csv(R + 'tier2a_paired_edgeR_QL_min10.csv')
cf = pd.read_csv(R + 'g2m_camera_fry_min10.csv')
g2 = json.load(open(R + 'g2m_per_patient_scores.json'))
legacy_pct = {'TN_0106': 29.7, 'TN_0114_T2': 58.8, 'TN_0126': 93.9, 'TN_0135': 82.0,
              'TN_B1_0131': 7.6, 'TN_B1_0177': 46.0, 'TN_B1_0554': 51.1, 'TN_B1_4031': 48.3}

ct = pd.crosstab(lab.patient, lab.label)
for c in ['malignant', 'non_malignant', 'unresolved']:
    if c not in ct: ct[c] = 0
ct = ct[['malignant', 'non_malignant', 'unresolved']]
pts = list(ct.index)

fig = plt.figure(figsize=(11, 14.2))
gs = fig.add_gridspec(5, 4, height_ratios=[1.05, .92, .92, 1.2, 1.15], hspace=.72, wspace=.62)

# ---- A: label composition ----
ax = fig.add_subplot(gs[0, :2])
y = np.arange(len(pts))[::-1]
tot = ct.sum(axis=1).values
left = np.zeros(len(pts))
for col, c, name in [('malignant', BLUE, 'malignant'), ('non_malignant', ORANGE, 'non-malignant'), ('unresolved', GRAY, 'unresolved')]:
    v = ct[col].values / tot * 100
    ax.barh(y, v, left=left, color=c, height=.62, label=name)
    for yy, vv, ll, n in zip(y, v, left, ct[col].values):
        if vv > 9: ax.text(ll + vv / 2, yy, f'{n:,}', ha='center', va='center', color='white', fontsize=6.8)
    left += v
ax.set_yticks(y); ax.set_yticklabels([nice(p) for p in pts], fontsize=7.5)
ax.set_xlim(0, 100); ax.set_xlabel('epithelial cells (%)')
ax.set_title('A  Malignancy calls from the inferCNV package (v1.18.1)', loc='left', fontweight='bold', fontsize=9, pad=18)
for s in ['top', 'right']: ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=7, ncol=3, loc='lower center', bbox_to_anchor=(.5, .99))

# ---- B: new vs legacy ----
ax = fig.add_subplot(gs[0, 2:])
new_pct = (100 * ct.malignant / (ct.malignant + ct.non_malignant)).values
old_pct = np.array([legacy_pct[p] for p in pts])
for yy, a, b in zip(y, old_pct, new_pct):
    ax.plot([a, b], [yy, yy], color=GRAY, lw=1.4, zorder=1, solid_capstyle='round')
ax.scatter(old_pct, y, color=GREEN, s=30, zorder=2, label='manuscript (legacy labels)')
ax.scatter(new_pct, y, color=BLUE, s=30, zorder=2, label='this reanalysis (inferCNV 1.18.1)')
for yy, a, b in zip(y, old_pct, new_pct):
    ax.annotate(f'{a:.0f}', (a, yy), textcoords='offset points', xytext=(0, -11), ha='center', fontsize=6.3, color=GREEN)
    ax.annotate(f'{b:.0f}', (b, yy), textcoords='offset points', xytext=(0, 6), ha='center', fontsize=6.3, color=BLUE)
ax.set_yticks([]); ax.set_xlim(-6, 112); ax.set_xlabel('malignant epithelial cells (%)')
ax.set_title('B  Same patients, two label sets', loc='left', fontweight='bold', fontsize=9, pad=18)
for s in ['top', 'right', 'left']: ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=7, ncol=2, loc='lower center', bbox_to_anchor=(.5, .99), handletextpad=.3, columnspacing=1.2)

# ---- C: CNV score vs correlation ----
cax = []
for k, p in enumerate(pts):
    ax = fig.add_subplot(gs[1 + k // 4, k % 4]); cax.append(ax)
    g = lab[lab.patient == p]
    for name, c in [('unresolved', GRAY), ('non_malignant', ORANGE), ('malignant', BLUE)]:
        m = (g.label == name).values
        ax.scatter(g.cnv_score[m] * 1000, g.cnv_corr[m], s=2, color=c, alpha=.5, linewidths=0)
    ax.axhline(.5, color=MUTED, lw=.6, ls='--')
    ax.set_ylim(-.35, 1.02)
    ax.set_title(f'{nice(p)}  (n = {len(g):,})', fontsize=7.5, color=INK, pad=3)
    if k % 4 == 0: ax.set_ylabel('CNV correlation', fontsize=7.5)
    else: ax.set_yticklabels([])
    if k // 4 == 1: ax.set_xlabel('CNV score (×10³)', fontsize=7.5)
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
_p = cax[0].get_position()
fig.text(_p.x0, _p.y1 + .022, 'C  Per-cell CNV score vs correlation with the patient consensus '
         '(dashed line: malignant cut-off r = 0.5)', fontweight='bold', fontsize=9, ha='left')

# ---- D: key genes, patient level (edgeR QL) ----
ax = fig.add_subplot(gs[3, :2])
kg = t2.set_index('gene').reindex(KEY_M + KEY_I[:3])
yy = np.arange(len(kg))[::-1]
sig = (kg.FDR < 0.05).values
ax.hlines(yy, 0, kg.logFC, color=GRAY, lw=1.2)
ax.scatter(kg.logFC, yy, s=36, color=[BLUE if s else 'white' for s in sig], edgecolors=BLUE, zorder=3, linewidths=1.1)
for v, t, f in zip(kg.logFC, yy, kg.FDR):
    off = 9 if v >= 0 else -9
    ax.annotate(f'FDR {f:.3f}', (v, t), textcoords='offset points', xytext=(off, 0), va='center',
                ha='left' if v >= 0 else 'right', fontsize=6.3, color=MUTED)
ax.axvline(0, color='#999999', lw=.6)
ax.set_yticks(yy); ax.set_yticklabels(kg.index, fontsize=7.5)
ax.set_xlabel('patient-level pseudo-bulk log$_2$FC (edgeR QL, paired, 6 patients)')
ax.set_xlim(-4.4, 7.4)
ax.set_title('D  Mitotic and immunoglobulin genes, patient level', loc='left', fontweight='bold', fontsize=9)
for s in ['top', 'right', 'left']: ax.spines[s].set_visible(False)
ax.legend(handles=[Line2D([], [], marker='o', ls='', color=BLUE, label='FDR < 0.05'),
                   Line2D([], [], marker='o', ls='', mfc='white', mec=BLUE, color=BLUE, label='FDR ≥ 0.05')],
          frameon=False, fontsize=7, loc='lower right', bbox_to_anchor=(1.02, -.02))

# ---- E: G2/M programme ----
ax = fig.add_subplot(gs[3, 2:])
d = g2['g2m_paired_diffs']
names = list(d); vals = [d[n] for n in names]
xx = np.arange(len(names))
ax.bar(xx, vals, color=BLUE, width=.55)
for x, v in zip(xx, vals):
    ax.text(x, v + .04, f'{v:+.2f}', ha='center', fontsize=7, color=INK)
ax.axhline(0, color='#999999', lw=.6)
ax.set_xticks(xx); ax.set_xticklabels([nice(n) for n in names], fontsize=6.8, rotation=20, ha='right')
ax.set_ylabel('G2/M score, malignant − non-malignant')
ax.set_ylim(0, max(vals) * 1.32)
camp = cf.PValue[cf.test == 'camera'].iloc[0]; fryp = cf.PValue[cf.test == 'fry'].iloc[0]
ax.set_title(f'E  G2/M programme, {g2["g2m_positive"]} patients positive', loc='left', fontweight='bold', fontsize=9, pad=16)
ax.annotate(f'paired t p = {g2["g2m_paired_t_p"]:.4f}   ·   camera p = {camp:.0e}   ·   fry p = {fryp:.0e}',
            xy=(0, 1.005), xycoords='axes fraction', fontsize=7, color=MUTED)
for s in ['top', 'right']: ax.spines[s].set_visible(False)

# ---- F: detection rate of the legacy gene lists ----
ax = fig.add_subplot(gs[4, :2])
leg = pd.read_csv(R + 'legacy_gene_lists_in_new_analysis.csv')
det = leg.drop_duplicates('gene').set_index('gene')['detection_rate_all_cells']
l174 = leg.loc[leg.list == 'legacy_stringent174', 'gene'].str.strip()
l869 = leg.loc[leg.list == 'legacy_broad869', 'gene'].str.strip()
# The grey reference curve covers every gene FindMarkers tested, so it needs the
# detection rate of each gene across all 54,820 cells. That is a whole-dataset
# intermediate rather than a result, and it is not redistributed here; set
# TNBC_DETECTION to a two-column CSV (gene, detection rate) to draw it. See
# docs/run_order.md for the one Seurat call that produces it.
_dpath = os.environ.get('TNBC_DETECTION', os.path.join(R, 'detection_rates.csv'))
series = [('manuscript 869-gene list', l869, GREEN),
          ('manuscript 174-gene list', l174, ORANGE)]
if os.path.exists(_dpath):
    det = pd.read_csv(_dpath, index_col=0).iloc[:, 0]
    series.insert(0, ('all 15,022 genes tested', pd.Series(t1.gene), GRAY))
else:
    print('note: %s not found; panel F omits the all-tested-genes reference curve'
          % _dpath, file=sys.stderr)
bins = np.logspace(-3.2, 0, 30)
for nm, L, c in series:
    s = det.reindex(L).dropna()
    ax.hist(s, bins=bins, density=True, histtype='step', lw=1.6, color=c, label=f'{nm} (median {s.median()*100:.1f}%)')
ax.axvline(0.01, color=INK, lw=.7, ls='--')
ax.annotate('Seurat min.pct = 1%', xy=(0.0105, .06), xycoords=('data', 'axes fraction'),
            fontsize=6.5, color=INK)
ax.set_xscale('log'); ax.set_xlabel('fraction of the 54,820 cells in which the gene is detected')
ax.set_ylabel('density'); ax.set_yticks([]); ax.set_ylim(0, ax.get_ylim()[1] * 1.22)
ax.set_title('F  How well expressed are the manuscript\'s Tier 1 genes?', loc='left', fontweight='bold', fontsize=9)
for s in ['top', 'right', 'left']: ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=6.8, loc='upper left', bbox_to_anchor=(-.02, 1.02))

# ---- G: cell-level rank vs patient-level rank ----
ax = fig.add_subplot(gs[4, 2:])
mg = t1[['gene', 'p_val', 'rank']].merge(t2[['gene', 'PValue', 'FDR']], on='gene')
mg['r_pat'] = mg.PValue.rank(method='first')
ax.scatter(mg['rank'], mg.r_pat, s=1.5, color=GRAY, alpha=.28, linewidths=0)
for genes, c, nm in [(KEY_M, BLUE, 'mitotic'), (KEY_I, ORANGE, 'immunoglobulin')]:
    s = mg[mg.gene.isin(genes)]
    ax.scatter(s['rank'], s.r_pat, s=26, color=c, zorder=3, linewidths=0, label=nm)
_off = {'PLK1': (7, -9), 'AURKA': (7, 5), 'IGKC': (9, -2), 'IGHG1': (8, 3), 'CDK1': (-4, -11)}
for _, r in mg[mg.gene.isin(_off)].iterrows():
    ax.annotate(r.gene, (r['rank'], r.r_pat), textcoords='offset points', xytext=_off[r.gene],
                fontsize=6.5, color=INK)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('rank at cell level (Seurat Wilcoxon)'); ax.set_ylabel('rank at patient level (edgeR QL)')
rho = mg['rank'].corr(mg.r_pat, method='spearman')
ax.set_title(f'G  The two tiers rank genes differently (Spearman ρ = {rho:.2f})', loc='left', fontweight='bold', fontsize=9)
for s in ['top', 'right']: ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=7, loc='lower right', bbox_to_anchor=(1.02, -.02))

fig.suptitle('Epithelial-only malignancy calling with the inferCNV package and harmonised Tier 1 / Tier 2 analysis\n'
             '8 TNBC tumours, 30,353 epithelial cells; 6 patients paired at ≥10 cells per group',
             fontsize=10.5, fontweight='bold', x=.09, ha='left', y=.982)
fig.savefig(OUT + 'figure_R_final.png', dpi=200, bbox_inches='tight')
fig.savefig(OUT + 'figure_R_final.pdf', bbox_inches='tight')
print('saved', OUT + 'figure_R_final.png')
