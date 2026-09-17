import sys, json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT = '/home/claude/cnv/out/'
BLUE, ORANGE, GREEN, GRAY, INK, MUTED = '#2E6BA8', '#D9822B', '#4C9A6B', '#B5B5B5', '#1a1a1a', '#6b6b6b'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.edgecolor': '#999999',
                     'axes.linewidth': .6, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.labelcolor': INK, 'text.color': INK, 'figure.facecolor': 'white'})
MINC = int(sys.argv[1]) if len(sys.argv) > 1 else 10
lab = pd.read_csv(OUT + 'epithelial_labels.csv', index_col=0)
cmp_ = pd.read_csv(OUT + 'label_comparison_by_patient.csv', index_col=0)
key = pd.read_csv(OUT + f'key_genes_tier1_vs_tier2_min{MINC}.csv', index_col=0)
summ = json.load(open(OUT + f'summary_min{MINC}.json'))
pts = list(cmp_.index)
nice = lambda p: p.replace('TN_', 'TN-').replace('_', '-')

fig = plt.figure(figsize=(11, 11.4))
gs = fig.add_gridspec(4, 4, height_ratios=[1.15, .95, .95, 1.3], hspace=.85, wspace=.65)

# ---- A: label composition per patient ----
ax = fig.add_subplot(gs[0, :2])
y = np.arange(len(pts))[::-1]
tot = cmp_[['malignant', 'non_malignant', 'unresolved']].sum(axis=1).values
left = np.zeros(len(pts))
for col, c, name in [('malignant', BLUE, 'malignant'), ('non_malignant', ORANGE, 'non-malignant'), ('unresolved', GRAY, 'unresolved')]:
    v = cmp_[col].values / tot * 100
    ax.barh(y, v, left=left, color=c, height=.62, label=name)
    for yy, vv, ll, n in zip(y, v, left, cmp_[col].values):
        if vv > 9: ax.text(ll + vv / 2, yy, f'{n:,}', ha='center', va='center', color='white', fontsize=6.8)
    left += v
ax.set_yticks(y); ax.set_yticklabels([nice(p) for p in pts], fontsize=7.5)
ax.set_xlim(0, 100); ax.set_xlabel('epithelial cells (%)')
ax.set_title('A  Harmonised CNV labels (this reanalysis)', loc='left', fontweight='bold', fontsize=9, pad=18)
for s in ['top', 'right']: ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=7, ncol=3, loc='lower center', bbox_to_anchor=(.5, .99))

# ---- B: new vs legacy malignant fraction ----
ax = fig.add_subplot(gs[0, 2:])
new_pct = (100 * cmp_.malignant / (cmp_.malignant + cmp_.non_malignant)).values
old_pct = cmp_.legacy_malignant_pct.values
for yy, a, b in zip(y, old_pct, new_pct):
    ax.plot([a, b], [yy, yy], color=GRAY, lw=1.4, zorder=1, solid_capstyle='round')
ax.scatter(old_pct, y, color=GREEN, s=30, zorder=2, label='manuscript (legacy inferCNV)')
ax.scatter(new_pct, y, color=BLUE, s=30, zorder=2, label='this reanalysis')
for yy, a, b in zip(y, old_pct, new_pct):
    ax.annotate(f'{a:.0f}', (a, yy), textcoords='offset points', xytext=(0, -11), ha='center', fontsize=6.3, color=GREEN)
    ax.annotate(f'{b:.0f}', (b, yy), textcoords='offset points', xytext=(0, 6), ha='center', fontsize=6.3, color=BLUE)
ax.set_yticks([]); ax.set_xlim(-6, 112); ax.set_xlabel('malignant epithelial cells (%)')
ax.set_title('B  Same patients, two label sets', loc='left', fontweight='bold', fontsize=9, pad=18)
for s in ['top', 'right', 'left']: ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=7, loc='lower left', bbox_to_anchor=(-.02, -.03), handletextpad=.3)

# ---- C: CNV score vs correlation (2 x 4) ----
sc = pd.read_csv(OUT + 'cnv_cell_scores.csv', index_col=0)
e = sc[~sc.is_ref]
for k, p in enumerate(pts):
    ax = fig.add_subplot(gs[1 + k // 4, k % 4])
    g = e[e.patient == p]; lb = lab.loc[g.index, 'label']
    for name, c in [('unresolved', GRAY), ('non_malignant', ORANGE), ('malignant', BLUE)]:
        m = (lb == name).values
        ax.scatter(g.cnv_score[m] * 1000, g.cnv_corr[m], s=2, color=c, alpha=.5, linewidths=0)
    ax.axhline(.5, color=MUTED, lw=.6, ls='--')
    ax.set_xlim(0, 35); ax.set_ylim(-.35, 1.02)
    ax.set_title(f'{nice(p)}  (n = {len(g):,})', fontsize=7.5, color=INK, pad=3)
    if k % 4 == 0: ax.set_ylabel('CNV correlation', fontsize=7.5)
    else: ax.set_yticklabels([])
    if k // 4 == 1: ax.set_xlabel('CNV score (×10³)', fontsize=7.5)
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
fig.text(.09, .672, 'C  Per-cell CNV score vs correlation with the patient consensus (dashed line: malignant cut-off r = 0.5)',
         fontweight='bold', fontsize=9, ha='left')

# ---- D: key genes at the patient level ----
ax = fig.add_subplot(gs[3, :2])
kg = key.loc[[g for g in ['PLK1','AURKA','AURKB','TTK','NEK2','MELK','CDK1','CCNB1','CCNB2','GTSE1','TOP2A','MKI67',
                          'IGKC','IGLC2','IGHG1'] if g in key.index]]
yy = np.arange(len(kg))[::-1]
sig = kg.paired_FDR < 0.05
ax.hlines(yy, 0, kg.paired_logFC, color=GRAY, lw=1.2)
ax.scatter(kg.paired_logFC, yy, s=36, color=[BLUE if s else 'white' for s in sig],
           edgecolors=BLUE, zorder=3, linewidths=1.1)
for v, t, f in zip(kg.paired_logFC, yy, kg.paired_FDR):
    off = 9 if v >= 0 else -9
    ax.annotate(f'FDR {f:.3f}', (v, t), textcoords='offset points', xytext=(off, 0), va='center',
                ha='left' if v >= 0 else 'right', fontsize=6.3, color=MUTED)
ax.axvline(0, color='#999999', lw=.6)
ax.set_yticks(yy); ax.set_yticklabels(kg.index, fontsize=7.5)
ax.set_xlabel('patient-level pseudo-bulk log$_2$FC (paired, malignant vs non-malignant)')
ax.set_xlim(-4.2, 6.4)
ax.set_title('D  Mitotic and immunoglobulin genes, patient level',
             loc='left', fontweight='bold', fontsize=9)
for s in ['top', 'right', 'left']: ax.spines[s].set_visible(False)
ax.legend(handles=[Line2D([], [], marker='o', ls='', color=BLUE, label='FDR < 0.05'),
                   Line2D([], [], marker='o', ls='', mfc='white', mec=BLUE, color=BLUE, label='FDR ≥ 0.05')],
          frameon=False, fontsize=7, loc='lower right')

# ---- E: G2/M programme score ----
ax = fig.add_subplot(gs[3, 2:])
d = summ['g2m_paired_diffs']
names = list(d); vals = [d[n] for n in names]
xx = np.arange(len(names))
ax.bar(xx, vals, color=BLUE, width=.55)
for x, v in zip(xx, vals):
    ax.text(x, v + .04, f'{v:+.2f}', ha='center', fontsize=7, color=INK)
ax.axhline(0, color='#999999', lw=.6)
ax.set_xticks(xx); ax.set_xticklabels([nice(n) for n in names], fontsize=6.8, rotation=20, ha='right')
ax.set_ylabel('G2/M score, malignant − non-malignant')
ax.set_ylim(0, max(vals) * 1.3)
sub = f"paired t p = {summ['g2m_paired_t_p']:.4f}" if summ.get('g2m_paired_t_p') else ''
ax.set_title(f'E  Patient-level G2/M programme ({summ["g2m_positive"]} patients positive; {sub})',
             loc='left', fontweight='bold', fontsize=9)
for s in ['top', 'right']: ax.spines[s].set_visible(False)

npat = len(summ['paired_patients'])
fig.suptitle(f'Epithelial-only CNV labels and harmonised Tier 1 / Tier 2 analysis (8 patients; {npat} paired, ≥{MINC} cells per group)',
             fontsize=10.5, fontweight='bold', x=.09, ha='left', y=.975)
fig.savefig(OUT + f'figure_full8_min{MINC}.png', dpi=200, bbox_inches='tight')
fig.savefig(OUT + f'figure_full8_min{MINC}.pdf', bbox_inches='tight')
print('saved', OUT + f'figure_full8_min{MINC}.png')
