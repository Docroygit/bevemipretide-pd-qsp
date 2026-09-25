"""Corrected figures for the BJP manuscript, drawn only from saved results in ../results."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from model import ROOT

R, F = ROOT / 'results', ROOT / 'figures'
F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.titlesize': 11, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'savefig.dpi': 300,
                     'svg.fonttype': 'none', 'legend.frameon': False})
ARM = {'placebo': ('Placebo', '#8C8C8C'), '10': ('10 mg', '#9CC3E4'), '30': ('30 mg', '#2F7FC1'), '60': ('60 mg', '#16426B')}
PATH = {'mitochondrial': '#C0504D', 'inflammation': '#E0A030', 'proteotoxic': '#6A5A9E'}
SCEN = {'healthy': ('Healthy', '#4E8A4A', ':'), 'disease': ('Untreated disease', '#A54452', '--'),
        'A_30mg': ('30 mg, scenario A', '#2F7FC1', '-'), 'B_30mg': ('30 mg, scenario B', '#2F7FC1', (0, (4, 2)))}


def load(name): return json.loads((R / name).read_text(encoding='utf-8'))


def save(fig, name):
    for ext in ('png', 'svg'): fig.savefig(F / f'{name}.{ext}', bbox_inches='tight', facecolor='white')
    plt.close(fig)


def fig5():
    pk = pd.read_csv(R / 'pk_accumulation.csv'); audit = load('pk.json')
    traj = pd.read_csv(R / 'human_20y_trajectories.csv'); hn = load('human_nominal.json')
    st4 = pd.read_csv(R / 'ST4_phased.csv')
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw={'width_ratios': [1, 1.25, 1]})
    ax = axs[0]
    for sp, col in [('mouse', '#B78629'), ('human', '#2F7FC1')]:
        d = pk[pk.species == sp]
        ax.plot(d.day, d.fraction_of_SS_trough, color=col,
                label=f"{sp.capitalize()} (slowest t½ {audit[sp]['slowest_mode_half_life_days']:.0f} d)")
    ax.axvline(120, color='#555', ls=':', lw=1); ax.text(124, .08, 'day 120', fontsize=8, color='#555')
    ax.set(xlim=(0, 730), ylim=(0, 1.02), xlabel='Days of once-daily dosing', ylabel='Mitochondrial trough / steady-state trough',
           title='A  Accumulation to steady state'); ax.legend(loc='lower right', fontsize=8)
    ax = axs[1]
    for key, (lab, col, ls) in SCEN.items():
        d = traj[traj.scenario == key]; ax.plot(d.year, d.DA, color=col, ls=ls, label=lab, lw=1.8)
    for thr, lab in [(.7, 'symptom onset'), (.5, 'clinical PD'), (.4, 'Bernheimer')]:
        ax.axhline(thr, color='#999', lw=.8, ls='-'); ax.text(20.2, thr, lab, va='center', fontsize=8, color='#666')
    ax.set(xlim=(0, 20), ylim=(0, 1.02), xlabel='Years', ylabel='DA neuron survival fraction',
           title='B  Human disease course, treatment from year 0'); ax.legend(loc='lower left', fontsize=8)
    ax = axs[2]
    sub = st4[st4.scenario == 'A']; delays = sorted(sub.delay_months.unique()); w = .25
    for i, mg in enumerate([10, 30, 60]):
        vals = [sub[(sub.delay_months == dl) & (sub.dose_mg == mg)].delta_DA.iloc[0] for dl in delays]
        ax.bar(np.arange(len(delays)) + (i - 1) * w, vals, w, color=ARM[str(mg)][1], label=ARM[str(mg)][0])
    ax.set_xticks(range(len(delays)), [f'{d} mo' for d in delays])
    ax.set(xlabel='Diagnostic delay after crossing 70% survival', ylabel='ΔDA vs placebo after 18 months',
           title='C  Benefit by dose and diagnostic delay'); ax.legend(fontsize=8)
    fig.tight_layout(); save(fig, 'Fig5_human_translation_corrected')


def fig6():
    e = pd.read_csv(R / 'vct_endpoints.csv'); s = load('vct_summary.json')
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw={'width_ratios': [1.1, 1.2, .9]})
    ax = axs[0]; arms = ['placebo', 'A_10', 'A_30', 'A_60']
    data = [e[e.arm == a].DA_18mo.values for a in arms]
    parts = ax.violinplot(data, showextrema=False, widths=.8)
    for body, a in zip(parts['bodies'], arms):
        body.set_facecolor(ARM[a.split('_')[-1]][1]); body.set_alpha(.45); body.set_edgecolor('none')
    ax.boxplot(data, widths=.18, showfliers=False, medianprops={'color': 'k'})
    ax.axhline(.7, color='#999', lw=.8, ls=':'); ax.text(4.45, .7, 'enrolment\nthreshold', fontsize=7, va='center', color='#666')
    ax.set_xticks(range(1, 5), [ARM[a.split('_')[-1]][0] for a in arms])
    ax.set(ylabel='DA neuron survival at 18 months', title=f"A  Survival at 18 months (n = {s['flow']['included_analysed']})")
    ax = axs[1]
    pl = e[e.arm == 'placebo'].set_index('id').DA_18mo; tr = e[e.arm == 'A_30'].set_index('id').DA_18mo
    ben = np.sort((tr - pl).values)[::-1] * 100
    ax.bar(np.arange(len(ben)), ben, width=1, color=np.where(ben > 5, '#2F7FC1', '#9CC3E4'))
    ax.axhline(5, color='#C0504D', ls='--', lw=1); ax.text(len(ben) * .6, 5.3, '5-point simulation threshold', color='#C0504D', fontsize=8)
    ax.set(xlabel='Virtual patient (ranked)', ylabel='ΔDA vs own placebo (percentage points)',
           title='B  Individual benefit, 30 mg (scenario A)', xlim=(-2, len(ben) + 2))
    ax = axs[2]
    rows = ['mitochondrial', 'inflammation', 'proteotoxic', 'total']
    mat = np.array([[s['arms'][f'A_{mg}']['rescue_fraction_cumulative_loss'][r] * 100 for mg in (10, 30, 60)] for r in rows])
    im = ax.imshow(mat, cmap='Greens', vmin=0, vmax=max(40, mat.max()), aspect='auto')
    for i in range(len(rows)):
        for j in range(3): ax.text(j, i, f'{mat[i, j]:.0f}%', ha='center', va='center', fontsize=10,
                                   color='white' if mat[i, j] > .6 * mat.max() else 'black')
    ax.set_xticks(range(3), ['10 mg', '30 mg', '60 mg']); ax.set_yticks(range(len(rows)), ['Mitochondrial', 'Neuroinflammation', 'Proteotoxicity', 'Total'])
    ax.set_title('C  Reduction in cumulative neuron loss'); ax.spines[:].set_visible(False)
    fig.colorbar(im, ax=ax, fraction=.05, label='% of placebo loss prevented')
    fig.tight_layout(); save(fig, 'Fig6_virtual_trial_corrected')


def fig7():
    env = pd.read_csv(R / 'envelope.csv'); envs = load('envelope.json'); s = load('vct_summary.json')
    e = pd.read_csv(R / 'vct_endpoints.csv')
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
    ax = axs[0]; ok = env[env.status == 'included']
    for scen, col in [('A', '#2F7FC1'), ('B', '#9CC3E4')]:
        ax.hist(ok[f'DA_saved_{scen}'], bins=20, color=col, alpha=.7, label=f'Scenario {scen}')
    ax.set(xlabel='DA saved at 18 months, 30 mg vs placebo', ylabel='Grid combinations',
           title=f"A  Translational parameter grid ({len(ok)}/{envs['grid_points']})"); ax.legend(fontsize=8)
    ax = axs[1]; pl = e[e.arm == 'placebo'].set_index('id')
    data = [(e[e.arm == f'A_{mg}'].set_index('id').DA_18mo - pl.DA_18mo).values * 100 for mg in (10, 30, 60)]
    ax.boxplot(data, widths=.4, showfliers=True, flierprops={'markersize': 2, 'alpha': .4}, medianprops={'color': 'k'})
    ax.set_xticks(range(1, 4), ['10 mg', '30 mg', '60 mg'])
    ax.set(ylabel='Matched ΔDA (percentage points)', title='B  Patient-level benefit (scenario A)')
    ax = axs[2]; y = 0; labels = []
    for ep, name in [('DA', 'DA'), ('motor', 'Motor')]:
        for scen, mk in [('A', 'o'), ('B', 's')]:
            for mg in (10, 30, 60):
                st = s['arms'][f'{scen}_{mg}'][ep]; lo, hi = st['d_unpaired_95']
                ax.errorbar(st['d_unpaired'], y, xerr=[[st['d_unpaired'] - lo], [hi - st['d_unpaired']]], fmt=mk,
                            color=ARM[str(mg)][1], ms=5, capsize=2)
                labels.append(f'{name} {mg} mg ({scen})'); y += 1
            y += .5; labels.append('')
    pos = []; k = 0
    for l in labels:
        pos.append(k); k += 1 if l else .5
    ax.set_yticks([p for p, l in zip(pos, labels) if l], [l for l in labels if l], fontsize=8)
    ax.axvline(0, color='#999', lw=.8); ax.invert_yaxis()
    ax.set(xlabel="Unpaired Cohen's d (95% bootstrap interval)", title='C  Standardised effect sizes')
    fig.tight_layout(); save(fig, 'Fig7_uncertainty_corrected')


def figS1():
    so = load('sobol_pm50.json'); names = so['parameters']
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
    for j, (ax, out) in enumerate(zip(axs, ['DA neuron survival, day 35', 'Motor score, day 35'])):
        st = np.array(so['ST'])[:, j]; s1 = np.array(so['S1'])[:, j]
        stci = np.array(so['ST_95'])[:, :, j]; order = np.argsort(st)
        yy = np.arange(len(names))
        ax.barh(yy, st[order], color='#2F7FC1', alpha=.35, label='Total order $S_T$')
        ax.errorbar(st[order], yy, xerr=[st[order] - stci[0][order], stci[1][order] - st[order]], fmt='none', ecolor='#2F7FC1', lw=1)
        ax.scatter(s1[order], yy, color='#16426B', s=18, zorder=3, label='First order $S_1$')
        ax.set_yticks(yy, [names[i] for i in order]); ax.axvline(0, color='#999', lw=.8)
        ax.set(xlabel='Sobol index', title=f"{'AB'[j]}  {out}")
    axs[0].legend(loc='lower right', fontsize=8)
    fig.suptitle(f"Sobol indices, ±50% uniform ranges, N = {so['n_base']} ({so['valid_rows']} valid rows)", fontsize=10, x=.01, ha='left')
    fig.tight_layout(); save(fig, 'FigS1_sobol_pm50_corrected')


def figS2():
    ts = pd.read_csv(R / 'mouse_timeseries.csv')
    fig, axs = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, (key, title) in zip(axs, [('disease', 'A  Untreated disease'), ('A_5mgkg', 'B  5 mg kg⁻¹ daily')]):
        d = ts[ts.scenario == key]
        ax.stackplot(d.day, d.D_mito, d.D_inflam, d.D_proteo, colors=[PATH['mitochondrial'], PATH['inflammation'], PATH['proteotoxic']],
                     labels=['Mitochondrial', 'Neuroinflammation', 'Proteotoxicity'], alpha=.85)
        ax.set(xlabel='Day', title=title, xlim=(0, 50))
    axs[0].set_ylabel('Cumulative neuron loss (fraction)'); axs[0].legend(loc='upper left', fontsize=8)
    fig.tight_layout(); save(fig, 'FigS2_cumulative_death_corrected')


def figS3():
    dr = pd.read_csv(R / 'mouse_dose_response_day50.csv')
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for scen, ls in [('A', '-'), ('B', '--')]:
        d = dr[dr.scenario == scen]
        axs[0].plot(d.dose_mg_per_kg, d.DA, 'o', ls=ls, color='#2F7FC1', label=f'Scenario {scen}', ms=4)
        axs[1].plot(d.dose_mg_per_kg, d.CL_ratio, 'o', ls=ls, color='#4E8A4A', label=f'Scenario {scen}', ms=4)
    axs[0].set(xlabel='Dose (mg kg⁻¹ day⁻¹)', ylabel='DA neuron survival, day 50', title='A  Neuron survival')
    axs[1].set(xlabel='Dose (mg kg⁻¹ day⁻¹)', ylabel='Functional CL ratio, day 50', title='B  Cardiolipin ratio')
    for ax in axs: ax.legend(fontsize=8)
    fig.tight_layout(); save(fig, 'FigS3_dose_response_corrected')


if __name__ == '__main__':
    for f in (fig5, fig6, fig7, figS1, figS2, figS3): f(); print('drew', f.__name__)
