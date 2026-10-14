"""
Figures for the independent validation (Scientific Data revision, round 2).

  python figures_validation.py field      # GBOV RM-7 plot PAI vs retrieved product, by biome
  python figures_validation.py gapfill    # gap-fill accuracy by QA code (field + hold-out)
  python figures_validation.py lp3        # 30 m product vs GBOV Sentinel-2 20 m maps
  python figures_validation.py modis      # MODIS consistency by vegetation type, common cell set

Figure files are named by their descriptive title (titles are set in the manuscript, not on
the figure).
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm

VAL = r"K:\Hangkai\CONUS_LAI_validation"
FIG = r"C:\Users\hyou34\OneDrive - UW-Madison\manuscripts\paper\CONUS LAI\SD R2\figures"
INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#898781', '#e1e0d9'
BLUE = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b']
SEQ = LinearSegmentedColormap.from_list('seqblue', BLUE)
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.edgecolor': '#c3c2b7',
                     'axes.labelcolor': INK2, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.linewidth': 0.6, 'savefig.dpi': 400})
BIOMES = ['Evergreen Needleleaf Forest', 'Deciduous Broadleaf Forest', 'Mixed Forest',
          'Closed Shrublands', 'Open Shrublands', 'Grassland', 'Croplands']
LETTERS = 'abcdefghijkl'
SHORT = {'Evergreen Needleleaf Forest': 'Evergreen needleleaf forest', 'Deciduous Broadleaf Forest':
         'Deciduous broadleaf forest', 'Mixed Forest': 'Mixed forest', 'Closed Shrublands':
         'Closed shrubland', 'Open Shrublands': 'Open shrubland', 'Croplands': 'Cropland'}


def metrics(r, e):
    m = np.isfinite(r) & np.isfinite(e); r, e = r[m], e[m]
    return dict(n=len(r), r2=np.corrcoef(r, e)[0, 1] ** 2, rmse=np.sqrt(np.mean((e - r) ** 2)),
                bias=np.mean(e - r))


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    p = os.path.join(FIG, name + '.png')
    fig.savefig(p, bbox_inches='tight', facecolor='white')
    print('saved', p)


def density_panel(ax, ref, est, title, lim=9, norm=None):
    m = np.isfinite(ref) & np.isfinite(est)
    hb = ax.hexbin(ref[m], est[m], gridsize=36, extent=(0, lim, 0, lim), cmap=SEQ,
                   norm=norm or LogNorm(1, 300), mincnt=1, linewidths=0.2, edgecolors='white')
    ax.plot([0, lim], [0, lim], color=MUTED, lw=1, ls='--', zorder=3)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect('equal')
    ax.set_xticks(range(0, lim + 1, 3)); ax.set_yticks(range(0, lim + 1, 3))
    ax.grid(color=GRID, lw=0.4); ax.set_axisbelow(True)
    s = metrics(ref, est)
    ax.text(0.04, 0.96, f"N = {s['n']:,}\nR² = {s['r2']:.2f}\nRMSE = {s['rmse']:.2f}\n"
            f"Bias = {s['bias']:+.2f}", transform=ax.transAxes, va='top', ha='left',
            color=INK, fontsize=7, linespacing=1.3)
    ax.set_title(title, color=INK, fontsize=7.5, loc='left', pad=4)
    return hb


def fig_field():
    d = pd.read_csv(os.path.join(VAL, 'gbov_rm7_matchups.csv'))
    d = d[d['ret'].notna()]
    fig, axs = plt.subplots(2, 4, figsize=(7.2, 3.9), constrained_layout=True)
    groups = [('All plots', d)] + [(b, d[d['IGBP_class'] == b]) for b in BIOMES]
    for k, (ax, (name, g)) in enumerate(zip(axs.flat, groups)):
        hb = density_panel(ax, g['pai'].values, g['ret'].values, f"({LETTERS[k]}) {SHORT.get(name, name)}")
        if k % 4: ax.set_yticklabels([])
        if k < 4: ax.set_xticklabels([])
    fig.supxlabel('GBOV field PAI (m$^2$ m$^{-2}$)', color=INK2, fontsize=8)
    fig.supylabel('CONUS-30m-LAI retrieved (m$^2$ m$^{-2}$)', color=INK2, fontsize=8)
    cb = fig.colorbar(hb, ax=axs, shrink=0.6, pad=0.01)
    cb.set_label('Plot-dates per bin', color=INK2); cb.outline.set_visible(False)
    save(fig, 'Validation of retrieved CONUS-30m-LAI against GBOV field PAI at 27 NEON sites by biome')

    # subsets: Kang et al. (2021) records vs records not used by Kang
    fig, axs = plt.subplots(1, 2, figsize=(4.4, 2.3), constrained_layout=True)
    for k, (ax, (sub, name)) in enumerate(zip(axs, [('kang', 'Records used by Kang et al. (2021)'),
                                                   ('new', 'Records not used by Kang et al.')])):
        g = d[d['subset'] == sub]
        hb = density_panel(ax, g['pai'].values, g['ret'].values, f"({LETTERS[k]}) {name}")
        if k: ax.set_yticklabels([])
    axs[0].set_ylabel('CONUS-30m-LAI retrieved (m$^2$ m$^{-2}$)')
    fig.supxlabel('GBOV field PAI (m$^2$ m$^{-2}$)', color=INK2, fontsize=8)
    cb = fig.colorbar(hb, ax=axs, shrink=0.8, pad=0.02)
    cb.set_label('Plot-dates per bin', color=INK2); cb.outline.set_visible(False)
    save(fig, 'Field validation for records used and not used by Kang et al 2021')


def fig_gapfill():
    d = pd.read_csv(os.path.join(VAL, 'gbov_rm7_matchups.csv'))
    h = pd.read_csv(os.path.join(VAL, 'gapfill_holdout_samples.csv'))
    rows = [('Field PAI: retrieved (QA 0)', metrics(d['pai'].values, d['ret'].values)),
            ('Field PAI: filled (QA 1)', metrics(d.loc[d['ret'].isna() & (d['qa'] == 1), 'pai'].values,
                                             d.loc[d['ret'].isna() & (d['qa'] == 1), 'gf'].values))]
    cases = [('Withholding: QA 1, short gaps', (h['qa'] == 1) & (h['scenario'] == 'random25')),
             ('Withholding: QA 1, long gaps', (h['qa'] == 1) & (h['scenario'] == 'keep2')),
             ('Withholding: QA 2', h['qa'] == 2), ('Withholding: QA 3', h['qa'] == 3),
             ('Withholding: QA 4', h['qa'] == 4)]
    for name, m in cases:
        if m.sum() >= 50:
            rows.append((name, metrics(h.loc[m, 'truth'].values, h.loc[m, 'filled'].values)))
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.8), constrained_layout=True)
    names = [r[0] for r in rows]; y = np.arange(len(rows))[::-1]
    for ax, key, xl, letter in [(axs[0], 'r2', 'R$^2$', 'a'), (axs[1], 'rmse', 'RMSE (m$^2$ m$^{-2}$)', 'b'),
                                (axs[2], 'bias', 'Mean bias (m$^2$ m$^{-2}$)', 'c')]:
        v = np.array([r[1][key] for r in rows])
        ax.barh(y, v, height=0.6, color=BLUE[4], edgecolor='white', linewidth=1)
        for yi, vi, r in zip(y, v, rows):
            lab = ('0.00' if abs(vi) < 0.005 else f"{vi:+.2f}") if key == 'bias' else f"{vi:.2f}"
            ax.text(vi + (0.02 if vi >= 0 else -0.02), yi, lab,
                    va='center', ha='left' if vi >= 0 else 'right', color=INK2, fontsize=7)
        ax.set_yticks(y); ax.set_yticklabels(names if letter == 'a' else [], color=INK2)
        ax.axvline(0, color='#c3c2b7', lw=0.8)
        ax.set_xlabel(xl); ax.grid(axis='x', color=GRID, lw=0.4); ax.set_axisbelow(True)
        ax.set_title(f"({letter})", loc='left', color=INK, fontsize=8)
        for s in ['top', 'right']: ax.spines[s].set_visible(False)
        lo, hi = ax.get_xlim(); ax.set_xlim(min(lo, -0.05 if key == 'bias' else 0) * (1.6 if key == 'bias' else 1.25), hi * (1.6 if key == 'bias' else 1.2))
        if key == 'r2': ax.set_xlim(0, 1.15)
        ax.axhline(y[1] - 0.5, color='#c3c2b7', lw=0.8, ls=':')
    for yi, r in zip(y, rows):
        axs[2].text(1.03, yi, f"N = {r[1]['n']:,}", transform=axs[2].get_yaxis_transform(),
                    va='center', color=MUTED, fontsize=7)
    save(fig, 'Accuracy of retrieved and gap-filled values by QA code from field data and withholding tests')


def fig_lp3():
    px = pd.read_csv(os.path.join(VAL, 'gbov_lp3_pixels.csv'))
    mp = pd.read_csv(os.path.join(VAL, 'gbov_lp3_maps.csv'))
    mp['sd_ratio'] = mp['sd_prod'] / mp['sd_ref']
    fig = plt.figure(figsize=(7.2, 3.0), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.1, 1.1])
    ax = fig.add_subplot(gs[0])
    density_panel(ax, px['ref'].values, px['prod'].values, '(a) 30 m pixels, all maps',
                  lim=9, norm=LogNorm(10, 30000))
    ax.set_xlabel('GBOV LP-3 LAI (m$^2$ m$^{-2}$)'); ax.set_ylabel('CONUS-30m-LAI (m$^2$ m$^{-2}$)')
    order = mp.groupby('site')['r'].median().sort_values().index
    for k, (col, xl, ref) in enumerate([('r', 'Pearson r (30 m pixels within each map)', None),
                                        ('sd_ratio', 'SD ratio, product / GBOV (log scale)', 1.0)]):
        ax = fig.add_subplot(gs[k + 1])
        data = [mp.loc[mp['site'] == s, col].dropna().values for s in order]
        bp = ax.boxplot(data, vert=False, widths=0.6, showfliers=False, patch_artist=True,
                        medianprops=dict(color=INK, lw=1), whiskerprops=dict(color=MUTED, lw=0.6),
                        capprops=dict(color=MUTED, lw=0.6))
        for b in bp['boxes']: b.set(facecolor=BLUE[2], edgecolor=BLUE[5], lw=0.6)
        ax.set_yticks(range(1, len(order) + 1)); ax.set_yticklabels(order if k == 0 else [], fontsize=6.5, color=INK2)
        if ref is not None:
            ax.axvline(ref, color=MUTED, lw=1, ls='--'); ax.set_xscale('log')
            ax.set_xticks([0.25, 0.5, 1, 2, 4, 8]); ax.set_xticklabels(['0.25', '0.5', '1', '2', '4', '8'])
            ax.set_xlim(0.2, 14)
        ax.set_xlabel(xl); ax.grid(axis='x', color=GRID, lw=0.4); ax.set_axisbelow(True)
        ax.set_title(f"({LETTERS[k + 1]}) " + ['Within-map correlation', 'Spatial variability'][k],
                     loc='left', color=INK, fontsize=7.5, pad=4)
        for s in ['top', 'right']: ax.spines[s].set_visible(False)
    save(fig, 'Comparison of CONUS-30m-LAI with GBOV Sentinel-2 20 m reference LAI maps')


PFT_NAMES = {'shrub': 'Shrubland', 'grass': 'Grassland and pasture', 'decid': 'Deciduous forest',
             'crop': 'Cropland', 'mixed': 'Mixed forest', 'herbwet': 'Herbaceous wetland',
             'evergreen': 'Evergreen forest', 'woodywet': 'Woody wetland', 'ALL': 'All classes'}
SERIES = {'Rc': ('Retrieved', '#2a78d6'), 'Fc': ('Gap-filled', '#eb6834')}   # categorical slots 1-2


def modis_common():
    """MODIS consistency results on the common cell set, R1 settings (agg_validate_r2.py)."""
    d = pd.read_csv(os.path.join(VAL, 'agg_validate_r2_results.csv'))
    return d[(d['threshold'] == 'px50_R1') & (d['period'] == 'ALL')]


def fig_modis():
    """MODIS consistency by vegetation type on the common cell set. Every bar is labelled with
       its value, and each class label carries its sample size: N common cells (all cells with
       a gap-filled value in parentheses)."""
    d = modis_common()
    c = d[d['product'].isin(['Rc', 'Fc'])].set_index(['product', 'pft'])
    n_gf = d[d['product'] == 'F'].set_index('pft')['N']             # gap-filled product, own cells
    order = (c.loc['Rc'].drop('ALL')['R2'].sort_values(ascending=False).index.tolist()) + ['ALL']
    y = np.arange(len(order))[::-1].astype(float)
    y[-1] -= 0.5                                                   # separate the "All" row
    ylabels = [PFT_NAMES[p] + "\n" + f"N = {int(c.loc[('Rc', p), 'N']):,} ({int(n_gf[p]):,})"
               for p in order]
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 4.4), constrained_layout=True)
    specs = [('R2', 'R$^2$', 1.0, [0, 0.25, 0.5, 0.75, 1.0], '{:.2f}'),
             ('RMSE', 'RMSE (m$^2$ m$^{-2}$)', 1.25, [0, 0.5, 1.0], '{:.2f}'),
             ('bias', 'Mean bias (m$^2$ m$^{-2}$)', 0.65, [0, 0.2, 0.4, 0.6], '{:+.2f}')]
    for k, (ax, (key, xl, lim, ticks, fmt)) in enumerate(zip(axs, specs)):
        for j, (prod, (lab, col)) in enumerate(SERIES.items()):
            off = 0.19 if j == 0 else -0.19
            v = c.loc[prod].loc[order, key].values
            ax.barh(y + off, v, height=0.36, color=col, label=lab, edgecolor='white', linewidth=0.8)
            for yi, vi in zip(y, v):
                ax.text(vi + lim * 0.02, yi + off, fmt.format(vi), va='center', ha='left',
                        color=INK2, fontsize=6)
        ax.axhline(y[-2] - 0.75, color='#c3c2b7', lw=0.8, ls=':')
        ax.set_xlim(0, lim * 1.3); ax.set_xticks(ticks); ax.spines['bottom'].set_bounds(0, ticks[-1])
        ax.set_xlabel(xl)
        ax.set_yticks(y)
        ax.set_yticklabels(ylabels if k == 0 else [], color=INK2, fontsize=6.5, linespacing=1.15)
        ax.grid(axis='x', color=GRID, lw=0.4); ax.set_axisbelow(True)
        for s_ in ['top', 'right']: ax.spines[s_].set_visible(False)
        ax.set_title(f"({LETTERS[k]})", loc='left', color=INK, fontsize=8)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=7.5)
    save(fig, 'Consistency of CONUS-30m-LAI with MODIS LAI by vegetation type on a common sample set')


if __name__ == '__main__':
    for what in sys.argv[1:] or ['field']:
        {'field': fig_field, 'gapfill': fig_gapfill, 'lp3': fig_lp3, 'modis': fig_modis}[what]()
