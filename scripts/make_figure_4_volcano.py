"""Figure 4: cell-level discovery volcano plots (Tier 1).

Built from `tier1_cell_level_findmarkers.csv`, the Seurat FindMarkers output shipped in
this repository, so the figure is reproducible from repository contents alone.

The y axis is the Bonferroni-adjusted P value that Seurat returns in `p_val_adj`, which is
the quantity the significance counts in the manuscript are based on. Forty-one genes have
an adjusted P that underflows to zero in double precision; they are drawn at the axis cap.

Palette matches Figures 3, 5 and 6.
"""
import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
R = os.environ.get('TNBC_RESULTS',
                   os.path.join(_HERE, '..', 'results', 'epithelial_only_pseudobulk_R'))
OUT = os.environ.get('TNBC_FIGOUT', os.path.join(_HERE, '..', 'figures'))
os.makedirs(OUT, exist_ok=True)

BLUE, RED, GREY, MUTED = '#2E6BA8', '#C0392B', '#C8C8C8', '#6b6b6b'
CAP = 300.0                      # axis cap for -log10 adjusted P
ALPHA = 0.05

KEY_I = ['IGKC', 'IGLC2', 'IGHG1', 'IGHA1', 'SPP1']
KEY_M = ['PLK1', 'AURKA', 'AURKB', 'TTK', 'NEK2', 'MELK',
         'CDK1', 'CCNB1', 'CCNB2', 'GTSE1', 'TOP2A', 'MKI67']

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.linewidth': 0.7, 'xtick.major.width': 0.7,
                     'ytick.major.width': 0.7, 'figure.facecolor': 'white'})


def load():
    t = pd.read_csv(os.path.join(R, 'tier1_cell_level_findmarkers.csv'))
    t = t.sort_values('p_val').reset_index(drop=True)
    t['y'] = np.where(t.p_val_adj > 0, -np.log10(t.p_val_adj.clip(lower=1e-320)), CAP)
    t['y'] = t.y.clip(upper=CAP)
    return t


def leaders(ax, tab, genes, colour, side, lo, hi, xpad):
    """Circle each gene and run a leader line out to a label stacked at the panel edge."""
    sub = tab[tab.gene.isin(genes)].sort_values('p_val')
    xlo, xhi = ax.get_xlim()
    xtext = xlo + xpad if side == 'left' else xhi - xpad
    ha = 'left' if side == 'left' else 'right'
    for ytext, (_, r) in zip(np.linspace(hi, lo, len(sub)), sub.iterrows()):
        ax.plot([r.avg_log2FC], [r.y], marker='o', ms=5, mfc='none',
                mec=colour, mew=0.9, zorder=5)
        ax.annotate(r.gene, xy=(r.avg_log2FC, r.y), xytext=(xtext, ytext),
                    ha=ha, va='center', color=colour, fontsize=7.5, zorder=6,
                    annotation_clip=False,
                    arrowprops=dict(arrowstyle='-', color=colour, lw=0.5,
                                    shrinkA=2, shrinkB=4))


def panel(ax, tab, lfc, letter):
    up = (tab.p_val_adj < ALPHA) & (tab.avg_log2FC > lfc)
    dn = (tab.p_val_adj < ALPHA) & (tab.avg_log2FC < -lfc)
    ns = ~(up | dn)
    for mask, colour, z in ((ns, GREY, 1), (dn, BLUE, 2), (up, RED, 2)):
        ax.scatter(tab.avg_log2FC[mask], tab.y[mask], s=2.0, c=colour,
                   linewidths=0, zorder=z)

    ax.axhline(-np.log10(ALPHA), ls='--', lw=0.7, color='#555555', zorder=3)
    for v in (-lfc, lfc):
        ax.axvline(v, ls='--', lw=0.7, color='#555555', zorder=3)

    ax.set_xlabel('average log$_2$ fold change (malignant / non-malignant)')
    ax.set_ylabel('$-$log$_{10}$ adjusted $P$')
    ax.set_ylim(-12, CAP + 12)
    ax.set_title(f'{letter}   |log$_2$FC| > {lfc:g} and adjusted $P$ < {ALPHA}: '
                 f'{int(up.sum() + dn.sum()):,} genes '
                 f'({int(up.sum()):,} up, {int(dn.sum()):,} down)',
                 loc='left', fontweight='bold', fontsize=9.5, pad=10)
    return int(up.sum()), int(dn.sum())


def main():
    tab = load()
    n_cap = int((tab.p_val_adj == 0).sum())
    note = (f'{n_cap} of {len(tab):,} genes have an adjusted $P$ below the '
            f'double-precision limit and are drawn at the axis cap')

    fig, axes = plt.subplots(2, 1, figsize=(6.6, 8.6))
    fig.subplots_adjust(hspace=0.55)

    for ax, lfc, letter in ((axes[0], 0.5, 'A'), (axes[1], 1.0, 'B')):
        n_up, n_dn = panel(ax, tab, lfc, letter)
        print(f'panel {letter}: {n_up + n_dn} significant ({n_up} up, {n_dn} down)')
        leaders(ax, tab, KEY_I, BLUE, 'left', 140, CAP, 0.3)
        if letter == 'B':
            leaders(ax, tab, KEY_M, RED, 'right', 107, CAP, 0.3)
        ax.annotate(note, xy=(0.5, -0.17), xycoords='axes fraction',
                    ha='center', va='top', fontsize=7, color=MUTED)

    png = os.path.join(OUT, 'fig4_volcano.png')
    fig.savefig(png, dpi=400, bbox_inches='tight', pad_inches=0.04)
    fig.savefig(os.path.join(OUT, 'fig4_volcano.pdf'),
                bbox_inches='tight', pad_inches=0.04)
    print('wrote', png)


if __name__ == '__main__':
    main()
