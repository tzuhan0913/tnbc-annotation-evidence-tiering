"""Rebuild Figure 4 (Tier 3 pathway reproducibility) and Figure 5 (Tier 5 GDSC calibration)
from the re-run outputs.  Palette matches Figure 2/3."""
import csv, os, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BLUE, ORANGE, GREEN, GREY = '#2E6BA8', '#D9822B', '#4C9A6B', '#8A8A8A'
A = os.environ.get('TNBC_RESULTS',
     os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results'))
OUT = os.environ.get('TNBC_FIGOUT',
       os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures'))
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.linewidth': 0.7, 'xtick.major.width': 0.7, 'ytick.major.width': 0.7})

# ----------------------------------------------------------------- Figure 4 --
chk = list(csv.DictReader(open(f'{A}/tier3_pathfindr_reanalysis/original_pathway_claims_recheck.csv')))
terms = list(csv.DictReader(open(f'{A}/tier3_pathfindr_reanalysis/pathfindr_all_enriched_terms_reanalysis.csv')))

claims, seen = [], set()
for r in chk:
    k = (r['claim'], r['database'])
    if k not in seen:
        seen.add(k); claims.append(k)
state = {(r['claim'], r['database'], r['input_set']): r['found'] == 'TRUE' for r in chk}

fig = plt.figure(figsize=(7.2, 5.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.90,
                      left=0.205, right=0.985, top=0.885, bottom=0.115)

# -- A: which of the 16 original claims each input recovers
ax = fig.add_subplot(gs[0, 0])
order = sorted(range(len(claims)), key=lambda i: (['KEGG', 'GO-BP', 'GO-CC', 'GO-MF'].index(claims[i][1]), i))
ys = np.arange(len(order))[::-1]
for y, i in zip(ys, order):
    cl, db = claims[i]
    for x, iset, col in ((0, 'patient', BLUE), (1, 'cell', ORANGE)):
        ok = state.get((cl, db, iset), False)
        ax.scatter(x, y, s=62, marker='o' if ok else 'X',
                   color=col if ok else 'white', edgecolor=col, linewidth=1.1, zorder=3)
ax.set_yticks(ys)
SHORT = {'protein serine/threonine kinase activity': 'Ser/Thr kinase activity'}
ax.set_yticklabels([SHORT.get(claims[i][0], claims[i][0]) for i in order])
ax.set_xticks([0, 1]); ax.set_xticklabels(['Patient-\nlevel', 'Cell-\nlevel'])
ax.set_xlim(-0.55, 1.55); ax.set_ylim(-0.8, len(order) - 0.2)
ax.grid(axis='y', color='#E4E4E4', linewidth=0.6, zorder=0)
ax.set_axisbelow(True)
# database brackets
prev, start = None, 0
for n, i in enumerate(order):
    db = claims[i][1]
    if prev is not None and db != prev:
        ax.axhline(ys[n] + 0.5, color='#BBBBBB', linewidth=0.6, linestyle=':')
    prev = db
ax.set_title('A   Original Figure 4 pathway claims recovered\n      (patient-level 16/16, cell-level 9/16)', loc='left',
             fontsize=8.2, fontweight='bold')
ax.legend(handles=[plt.Line2D([], [], marker='o', color=BLUE, linestyle='', label='recovered'),
                   plt.Line2D([], [], marker='X', color='white', markeredgecolor=ORANGE,
                              linestyle='', label='not recovered')],
          loc='lower center', bbox_to_anchor=(0.5, -0.155), ncol=2, frameon=False,
          handletextpad=0.2, columnspacing=1.0, fontsize=7.2)

# -- B: fold enrichment of the mitotic terms, patient-level input
def best(term_pat, db, iset='patient'):
    hits = [t for t in terms if t['input_set'] == iset and t['database'] == db
            and re.search(term_pat, t['Term_Description'], re.I)]
    if not hits: return None
    return min(hits, key=lambda t: float(t['lowest_p']))

show = [('kinetochore', 'GO-CC'), ('spindle', 'GO-CC'), ('midbody', 'GO-CC'), ('centrosome', 'GO-CC'),
        ('chromosome segregation', 'GO-BP'), ('mitotic cytokinesis', 'GO-BP'),
        ('mitotic cell cycle', 'GO-BP'), ('G2/M transition', 'GO-BP'),
        ('microtubule binding', 'GO-MF'), ('protein kinase binding', 'GO-MF'),
        ('^Cell cycle$', 'KEGG'), ('p53 signaling', 'KEGG')]
ax = fig.add_subplot(gs[0, 1])
labs, fe, lp, cols = [], [], [], []
DBC = {'KEGG': BLUE, 'GO-BP': GREEN, 'GO-CC': ORANGE, 'GO-MF': '#7B5EA7'}
for pat, db in show:
    h = best(pat, db)
    if h is None: continue
    d0 = h['Term_Description']
    labs.append(d0 if len(d0) <= 26 else d0[:25] + '\u2026'); fe.append(float(h['Fold_Enrichment']))
    lp.append(-np.log10(max(float(h['lowest_p']), 1e-300))); cols.append(DBC[db])
y = np.arange(len(labs))[::-1]
ax.barh(y, fe, color=cols, height=0.68)
for yy, f, p in zip(y, fe, lp):
    ax.text(f + 0.25, yy, f'{f:.1f}  (−log₁₀P {p:.0f})', va='center', fontsize=6.6, color='#444444')
ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=7.2)
ax.set_xlabel('Fold enrichment (patient-level input)')
ax.set_xlim(0, max(fe) * 1.72)
ax.set_title('B   Enrichment, patient-level input', loc='left',
             fontsize=8.2, fontweight='bold')
ax.legend(handles=[Patch(color=DBC[k], label=k) for k in ('KEGG', 'GO-BP', 'GO-CC', 'GO-MF')],
          loc='lower center', bbox_to_anchor=(0.5, -0.155), ncol=4, frameon=False,
          fontsize=7, handlelength=0.9, columnspacing=0.9, handletextpad=0.35)
fig.savefig(f'{OUT}/fig4_new.png', dpi=400, bbox_inches='tight', pad_inches=0.04)
plt.close(fig)
print('fig4 written; panel A rows', len(order), '| panel B bars', len(labs))

# ----------------------------------------------------------------- Figure 5 --
per = list(csv.DictReader(open(f'{A}/tier5_gdsc_reanalysis/gdsc_per_compound_summary.csv')))
cand = [r for r in csv.DictReader(open(f'{A}/tier5_gdsc_reanalysis/gdsc_exact_candidate_drug_match.csv'))]
ing = [r for r in cand if r['in_gdsc'] == 'TRUE']

fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.25], width_ratios=[1.05, 1.0],
                      hspace=0.46, wspace=0.30, left=0.145, right=0.985, top=0.925, bottom=0.075)

# -- A: molecule-level calibration under three collapse rules
ax = fig.add_subplot(gs[0, 0])
# Counts, background rates and two-sided Fisher P values are read from the Tier 5 output rather
# than typed in.  The two sets have DIFFERENT backgrounds (the background is "GDSC molecules
# outside the set"), so each bar gets its own dashed background line.
tvd = list(csv.DictReader(open(f'{A}/tier5_gdsc_reanalysis/tier5_target_axis_vs_drug_list.csv')))
def _row(setname, rule):
    return next(r for r in tvd if r['set'] == setname and r['collapse_rule'] == rule)
def _p(x):
    from decimal import Decimal, ROUND_HALF_UP
    d = Decimal(x)
    return str(d.quantize(Decimal('0.001') if d < Decimal('0.1') else Decimal('0.01'), ROUND_HALF_UP))
RULE_LABELS = [('meanz', 'mean z\nacross screens'), ('any', 'any screen\nfavorable'),
               ('all', 'every screen\nfavorable')]
RULES = []
for rule, lab in RULE_LABELS:
    a_, b_ = _row('mitotic target axis', rule), _row('Tier 4 drug list', rule)
    RULES.append((lab, (int(a_['favourable']), int(a_['n_molecules'])),
                  (int(b_['favourable']), int(b_['n_molecules'])),
                  100 * float(a_['background_fraction']), 100 * float(b_['background_fraction']),
                  _p(a_['fisher_two_sided']), _p(b_['fisher_two_sided'])))
w = 0.33
xs = np.arange(len(RULES))
ax.bar(xs - w/2, [100*a[0]/a[1] for _, a, _, _, _, _, _ in RULES], width=w, color=BLUE,
       label='mitotic target axis')
ax.bar(xs + w/2, [100*b[0]/b[1] for _, _, b, _, _, _, _ in RULES], width=w, color=ORANGE,
       label='Tier 4 drug list')
for x, (lab, a, b, bga, bgb, pa, pb) in zip(xs, RULES):
    for cx, bg in ((x - w/2, bga), (x + w/2, bgb)):
        ax.plot([cx - w/2 - 0.03, cx + w/2 + 0.03], [bg, bg], color='#555555', linestyle='--',
                linewidth=1.0, zorder=3)
    ax.text(x - w/2, 100*a[0]/a[1] + 1.8, f'{a[0]}/{a[1]}\nP={pa}', ha='center', fontsize=6.4)
    ax.text(x + w/2, max(100*b[0]/b[1], bgb) + 2.4, f'{b[0]}/{b[1]}\nP={pb}', ha='center',
            fontsize=6.4)
ax.set_xticks(xs); ax.set_xticklabels([r[0] for r in RULES], fontsize=6.9)
ax.set_xlabel('rule for collapsing repeated GDSC screens', fontsize=7.2)
ax.set_ylabel('TNBC-favorable molecules (%)')
ax.set_ylim(0, 84)
from matplotlib.lines import Line2D
_h, _l = ax.get_legend_handles_labels()
_h.append(Line2D([0], [0], color='#555555', linestyle='--', linewidth=1.0))
_l.append('background for that set and rule')
ax.legend(_h, _l, frameon=False, fontsize=6.6, loc='upper right', bbox_to_anchor=(1.0, 1.02),
          handlelength=1.6)
ax.set_title('A   Molecule-level calibration, three collapse rules',
             loc='left', fontsize=8.5, fontweight='bold')

# -- B: what survives from the Tier 4 candidate list
ax = fig.add_subplot(gs[0, 1])
steps = [('Candidate\ncompounds', 532, BLUE), ('Present in\nGDSC', 59, BLUE),
         ('TNBC-\nfavorable', 17, ORANGE), ('Survive BH\ncorrection', 1, '#B23A48')]  # molecule level
xs = np.arange(len(steps))
ax.bar(xs, [s_[1] for s_ in steps], color=[s_[2] for s_ in steps], width=0.62)
for x, (lab, v, _) in zip(xs, steps):
    ax.text(x, v * 1.45, f'{v}', ha='center', fontsize=7.6, fontweight='bold')
ax.set_yscale('log'); ax.set_ylim(0.6, 20000)
ax.set_xticks(xs); ax.set_xticklabels([s_[0] for s_ in steps], fontsize=6.9)
ax.set_ylabel('Compounds (log scale)')
ax.annotate('bortezomib — reaches the gene list\nthrough proteasome subunits,\nnot the mitotic axis',
            xy=(3, 3.2), xytext=(2.0, 1500), fontsize=6.4, color='#B23A48', ha='center',
            arrowprops=dict(arrowstyle='-|>', color='#B23A48', linewidth=0.7,
                            connectionstyle='arc3,rad=-0.25', shrinkB=3))
ax.set_title('B   Attrition of the Tier 4 drug list', loc='left', fontsize=8.5, fontweight='bold')

# -- C: molecule-level evidence within the target axis
# Panels A and B count MOLECULES, so this panel does too: one bar per molecule, showing the mean of
# the per-screen (median TNBC z - median other z).  The three molecules GDSC screened twice carry a
# marker per screen, which is how the alisertib disagreement stays visible without inflating n.
ax = fig.add_subplot(gs[1, :])
NM = {'CD532': 'CD532', 'TOZASERTIB': 'Tozasertib', 'BI-2536': 'BI-2536', 'GSK1070916': 'GSK1070916',
      'AT-7519': 'AT-7519', 'ZM-447439': 'ZM-447439', 'SL 0101-1': 'SL 0101-1',
      'DINACICLIB': 'Dinaciclib', 'RO-3306': 'RO-3306', 'GENENTECH CPD 10': 'Genentech Cpd 10',
      'GW843682X': 'GW843682X', 'CGP-60474': 'CGP-60474', 'NPK76-II-72-1': 'NPK76-II-72-1',
      'ALISERTIB': 'Alisertib', 'BX-795': 'BX-795', 'MITOXANTRONE': 'Mitoxantrone',
      'MK-1775': 'MK-1775', 'ETOPOSIDE': 'Etoposide', 'HG-5-113-01': 'HG-5-113-01',
      'MPS-1-IN-1': 'MPS-1-IN-1'}
key = lambda s_: re.sub(r'[^A-Z0-9]', '', s_.upper())
NM = {key(k): v for k, v in NM.items()}
mol = {}
for r in per:
    if r['is_axis'] != 'TRUE' or r['p_z_tnbc_lower'] in ('', 'NA'): continue
    mol.setdefault(key(r['DRUG_NAME']), []).append(r)
rows = []
for k, v in mol.items():
    dz = [float(x['med_z_tnbc']) - float(x['med_z_other']) for x in v]
    rows.append((NM.get(k, v[0]['DRUG_NAME'].title()), -sum(dz) / len(dz), [-d for d in dz]))
rows.sort(key=lambda t: t[1])
y = np.arange(len(rows))
TOP2 = {'Etoposide', 'Mitoxantrone'}
cols = [GREEN if nm in TOP2 else (BLUE if v > 0 else ORANGE) for nm, v, _ in rows]
ax.barh(y, [v for _, v, _ in rows], color=cols, height=0.66)
for yy, (_, _, screens) in zip(y, rows):
    if len(screens) > 1:
        ax.scatter(screens, [yy] * len(screens), s=16, facecolor='white',
                   edgecolor='#333333', linewidth=0.8, zorder=4)
ax.axvline(0, color='#333333', linewidth=0.8)
ax.set_yticks(y); ax.set_yticklabels([nm for nm, _, _ in rows], fontsize=7.0)
ax.set_xlabel('TNBC sensitivity advantage (\u2212\u0394 median published z-score);  '
              'right = TNBC panel more sensitive', fontsize=7.6)
ax.set_title('C   The 20 evaluable target-axis molecules', loc='left', fontsize=8.5, fontweight='bold')
ax.legend(handles=[Patch(color=BLUE, label='TNBC-favorable'),
                   Patch(color=ORANGE, label='TNBC-unfavorable'),
                   Patch(color=GREEN, label='TOP2 cytotoxic (negative control)'),
                   plt.Line2D([], [], marker='o', linestyle='', markerfacecolor='white',
                              markeredgecolor='#333333', markersize=4.5,
                              label='individual screens (molecules screened twice)')],
          loc='lower right', frameon=False, fontsize=6.9, handlelength=1.1)
lo = min(min(s_) for _, _, s_ in rows); hi = max(max(s_) for _, _, s_ in rows)
pad = (hi - lo) * 0.10
ax.set_xlim(lo - pad, hi + pad * 2.2)
ax.set_ylim(-0.8, len(rows) - 0.2)
fig.savefig(f'{OUT}/fig5_new.png', dpi=400, bbox_inches='tight', pad_inches=0.04)
plt.close(fig)
print('fig5 written; panel C molecules', len(rows),
      '| favourable', sum(1 for _, v, _ in rows if v > 0))
