"""Figure 6: what an annotation actually changes is which patients are comparable.
Cohort 1 = GSE161529 (primary); cohort 2 = GSE176078 (Wu et al.)."""
import csv, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BLUE, ORANGE, GREEN, GREY = '#2E6BA8', '#D9822B', '#4C9A6B', '#8A8A8A'
R = os.environ.get('TNBC_REP_RESULTS',
     os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'replication_GSE176078'))
OUT = os.environ.get('TNBC_FIGOUT',
       os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures'))
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.linewidth': 0.7, 'xtick.major.width': 0.7, 'ytick.major.width': 0.7})

C1 = [('TN-0106', 30.2, 29.7), ('TN-0114-T2', 50.3, 58.8), ('TN-0126', 99.5, 93.9),
      ('TN-0135', 99.8, 82.3), ('TN-B1-0131', 100.0, 7.6), ('TN-B1-0177', 65.0, 46.0),
      ('TN-B1-0554', 86.3, 51.1), ('TN-B1-4031', 99.3, 48.3)]
C2 = [(r['patient'], 100 * float(r['infercnv_malignant_fraction']),
       100 * float(r['author_malignant_fraction']))
      for r in csv.DictReader(open(f'{R}/headline_summary.csv'))
      if r['row_type'] == 'annotation_patient']

fig = plt.figure(figsize=(7.2, 6.8))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.2], width_ratios=[1.28, 1.0],
                      hspace=0.52, wspace=0.32, left=0.095, right=0.985, top=0.925, bottom=0.085)

# ---------------- A: per-patient malignant fraction, both cohorts ----------------
ax = fig.add_subplot(gs[0, :])
off, ticks, labels = 0, [], []
for name, rows, col in [('Cohort 1  GSE161529', C1, BLUE), ('Cohort 2  GSE176078', C2, GREEN)]:
    for k, (p, ours, theirs) in enumerate(rows):
        x = off + k
        ax.plot([x, x], [min(ours, theirs), max(ours, theirs)], color=GREY, linewidth=1.0, zorder=1)
        ax.scatter(x, ours, s=26, color=col, zorder=3)
        ax.scatter(x, theirs, s=26, color='white', edgecolor=col, linewidth=1.1, zorder=3)
        ticks.append(x); labels.append(p)
    gap = max(abs(a - b) for _, a, b in rows)
    ax.text(off + (len(rows) - 1) / 2, 104, f'{name}   —   largest gap {gap:.0f} pp',
            ha='center', va='bottom', fontsize=7.6, fontweight='bold', color=col)
    off += len(rows) + 1.6
ax.set_xticks(ticks); ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=6.6)
ax.set_ylabel('malignant fraction of\nepithelial cells (%)')
ax.set_ylim(-6, 136); ax.set_yticks([0, 25, 50, 75, 100])
ax.axvline(len(C1) + 0.3, color='#DDDDDD', linewidth=0.8)
ax.legend(handles=[plt.Line2D([], [], marker='o', linestyle='', color='#555555',
                              label='our pre-specified inferCNV rule'),
                   plt.Line2D([], [], marker='o', linestyle='', markerfacecolor='white',
                              markeredgecolor='#555555', label='the published labels')],
          loc='upper left', bbox_to_anchor=(0.0, 1.0), frameon=False, fontsize=7,
          handletextpad=0.3, ncol=2, columnspacing=1.4)
ax.set_title('A   Per-patient malignant fraction under two defensible annotations', loc='left',
             fontsize=8.5, fontweight='bold')

# ---------------- B: discoveries by eligible patient set ----------------
S = {r['model']: r for r in csv.DictReader(open(f'{R}/E_model_summary.csv'))}
steps = [('each annotation\nchooses its own\npatients (4 + 4)',
          int(S['authors|all']['n_FDR05']), int(S['infercnv|all']['n_FDR05'])),
         ('the same\n3 patients',
          int(S['authors|intersect']['n_FDR05']), int(S['infercnv|intersect']['n_FDR05'])),
         ('the same 2 patients,\n≥ 20 cells per group',
          int(S['authors|intersect_min20']['n_FDR05']), int(S['infercnv|intersect_min20']['n_FDR05']))]
ax = fig.add_subplot(gs[1, 0])
xs = np.arange(len(steps)); w = 0.34
ax.bar(xs - w/2, [s[1] for s in steps], width=w, color=BLUE, label='published labels')
ax.bar(xs + w/2, [s[2] for s in steps], width=w, color=ORANGE, label='our inferCNV rule')
for x, (_, a, b) in zip(xs, steps):
    ax.text(x - w/2, a + 16, str(a), ha='center', fontsize=7.6, fontweight='bold')
    ax.text(x + w/2, b + 16, str(b), ha='center', fontsize=7.6, fontweight='bold')
    ax.text(x, max(a, b) + 62, f'{max(a, b) / max(min(a, b), 1):.1f}×',
            ha='center', fontsize=7.2, color='#555555')
ax.set_xticks(xs); ax.set_xticklabels([s[0] for s in steps], fontsize=6.9)
ax.set_ylabel('genes at patient-level FDR < 0.05')
ax.set_ylim(0, 860)
ax.legend(frameon=False, fontsize=7, loc='upper right', handlelength=1.0)
ax.annotate('the extra patient is\nCID4465 (10 non-malignant cells)',
            xy=(-w/2, 690), xytext=(0.44, 610), fontsize=6.4, color=BLUE, ha='left',
            arrowprops=dict(arrowstyle='-|>', color=BLUE, linewidth=0.7,
                            connectionstyle='arc3,rad=0.25', shrinkB=2))
ax.annotate('the extra patient is\nCID4513 (614 malignant cells)',
            xy=(w/2, 62), xytext=(0.44, 430), fontsize=6.4, color=ORANGE, ha='left',
            arrowprops=dict(arrowstyle='-|>', color=ORANGE, linewidth=0.7,
                            connectionstyle='arc3,rad=0.20', shrinkB=2))
ax.set_title('B   Patient-level discoveries by eligible patient set',
             loc='left', fontsize=8.5, fontweight='bold')

# ---------------- C: effect-size concordance on identical patients ----------------
rows = list(csv.DictReader(open(f'{R}/E_intersection_gene_comparison.csv')))
xa = np.array([float(r['logFC_authors']) for r in rows])
xi = np.array([float(r['logFC_infercnv']) for r in rows])
sa = np.array([float(r['FDR_authors']) < 0.05 for r in rows])
si = np.array([float(r['FDR_infercnv']) < 0.05 for r in rows])
ax = fig.add_subplot(gs[1, 1])
both = sa & si; disc = sa ^ si; neither = ~(sa | si)
ax.scatter(xa[neither], xi[neither], s=1.5, color='#D8D8D8', linewidths=0, rasterized=True)
ax.scatter(xa[disc], xi[disc], s=5, color=ORANGE, linewidths=0, rasterized=True)
ax.scatter(xa[both], xi[both], s=5, color=BLUE, linewidths=0, rasterized=True)
lim = [min(xa.min(), xi.min()) - 0.4, max(xa.max(), xi.max()) + 0.4]
ax.plot(lim, lim, color='#666666', linewidth=0.7, linestyle='--')
ax.set_xlim(lim); ax.set_ylim(lim)
r_ = np.corrcoef(xa, xi)[0, 1]
ax.set_xlabel('log₂FC, published labels', fontsize=7.4)
ax.set_ylabel('log₂FC, our inferCNV rule', fontsize=7.4)
ax.text(0.04, 0.965, f'r = {r_:.3f}\nall {int(disc.sum())} discordant calls\nshare the sign of log₂FC',
        transform=ax.transAxes, va='top', fontsize=6.9)
ax.legend(handles=[Patch(color=BLUE, label='significant under both'),
                   Patch(color=ORANGE, label='significant under one only'),
                   Patch(color='#D8D8D8', label='neither')],
          loc='lower right', frameon=False, fontsize=6.6, handlelength=1.0)
ax.set_title('C   On identical patients, the two annotations\n      yield concordant effect sizes',
             loc='left', fontsize=8.5, fontweight='bold')

fig.savefig(f'{OUT}/fig6_draft.png', dpi=400)
plt.close(fig)
print('fig6 written | cohort 2 patients', len(C2), '| genes', len(rows),
      '| r = %.3f' % r_, '| discordant', int(disc.sum()))
