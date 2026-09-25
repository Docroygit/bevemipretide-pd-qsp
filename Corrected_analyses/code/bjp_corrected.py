"""Corrected versions of the BJP manuscript analyses, run on the independent Python rebuild.

model.py, fast_kernel.py and parameters.csv are unmodified copies from
'Bevemipretide QSP/Revisions_2026-09-18' (hashes in README.md). This script only
changes how analyses are set up, following the corrections listed in README.md:

1. Drug exposure is normalised by the exact steady-state mean mitochondrial amount at the
   reference dose of each species (5 mg/kg mouse, 30 mg human), not by C_ref = 4.43 / 2.15.
   Scenario A re-anchors k_drug so the published mouse drug response is reproduced exactly
   (k_drug * A/4.43 == k_drug_A * A/A_ss). Scenario B keeps the literal k_drug = 0.1.
2. Virtual patients are enrolled only after untreated DA survival falls below 0.70, followed
   by a 12-month diagnostic delay; every candidate and exclusion is recorded.
3. Pathway attribution uses cumulative loss (sums to total loss), not instantaneous shares.
4. Response thresholds are simulation thresholds, not validated MCIDs.
"""
import os
for _name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[_name] = '1'
import argparse, csv, json, time, platform
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.stats import nct, qmc, rankdata, t as tdist
from model import *
from fast_kernel import make_rhs

OUT = ROOT / 'results'
YEAR_H = 365.25 * 24
C_REF_MOUSE_BJP = 4.43
K_DRUG_BJP = 0.1
HUMAN_DOSES_MG = [10, 30, 60]
REF_DOSE_MG = 30.
DA_THRESHOLD = 0.70
SCREEN_YEARS = 30.
SEED = 20260923


def k_drug_anchored():
    return K_DRUG_BJP * parameters()['exposure_reference'] / C_REF_MOUSE_BJP


SCENARIOS = {'A': k_drug_anchored(), 'B': K_DRUG_BJP}


def save_json(name, obj):
    (OUT / name).write_text(json.dumps(obj, indent=2, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else float(x)),
                            encoding='utf-8')


def write_csv(name, fields, rows):
    with (OUT / name).open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f); w.writerow(fields); w.writerows(rows)


def cl_ratio(y):
    return y[0] / (y[0] + y[1])


def end(p, days, dose=0., y0=None):
    """Endpoint state including the three cumulative-loss integrals from the start of this run."""
    return simulate(p, days=days, dose=dose, y0=y0, points=2).y[:, -1]


def crossing(p, y0, threshold=DA_THRESHOLD, years=SCREEN_YEARS, dose=0.):
    pk = PK(p, dose) if dose else None
    event = lambda t, y: y[18] - threshold
    event.terminal = True; event.direction = -1
    sol = solve_ivp(make_rhs(p, pk, p['exposure_reference']), (0, years * YEAR_H), np.asarray(y0)[:23],
                    method='LSODA', events=event, rtol=1e-8, atol=1e-10, max_step=168.)
    if not sol.success: raise RuntimeError(sol.message)
    if not len(sol.t_events[0]): return None, None
    return sol.t_events[0][0] / YEAR_H, sol.y_events[0][0]


def with_kdrug(p, scenario):
    q = dict(p); q['k_drug'] = SCENARIOS[scenario]; return q


# ---------------------------------------------------------------- mouse (Table 2, Tiers 1-3)
def mouse():
    y0 = initial_state()
    healthy, disease = parameters(healthy=True), parameters()
    A = with_kdrug(disease, 'A')
    H = simulate(healthy, days=365, y0=y0, points=731)
    D = simulate(disease, days=365, y0=y0, points=731)
    day = lambda s, d: s.y[:, int(round(d * 2))]
    h35, d35, h50, d50 = day(H, 35), day(D, 35), day(H, 50), day(D, 50)

    # Exact equivalence of scenario A with the manuscript normalisation in the mouse.
    src = simulate(disease, days=50, dose=1, y0=y0, points=101, reference_mode='manuscript')
    anc = simulate(A, days=50, dose=1, y0=y0, points=101)
    equivalence = float(np.max(np.abs(src.y - anc.y)))

    drug = {}
    for scen in 'AB':
        q = with_kdrug(disease, scen)
        s = simulate(q, days=365, dose=1, y0=y0, points=731)
        x35, x50, x365 = day(s, 35), day(s, 50), s.y[:, -1]
        drug[scen] = {
            'DA_gap_rescued_day35': (x35[18] - d35[18]) / (h35[18] - d35[18]),
            'DA_relative_gain_day35': x35[18] / d35[18] - 1,
            'CL_ratio_gap_restored_day35': (cl_ratio(x35) - cl_ratio(d35)) / (cl_ratio(h35) - cl_ratio(d35)),
            'aSyn_olig_reduction_day35': 1 - x35[6] / d35[6],
            'motor_improvement_day35': 1 - x35[22] / d35[22],
            'directions_day35': {'DA_up': x35[18] > d35[18], 'CI_up': x35[8] > d35[8], 'aSyn_down': x35[6] < d35[6],
                                 'mROS_down': x35[13] < d35[13], 'motor_down': x35[22] < d35[22]},
            'day35': dict(zip(STATES, x35[:23])), 'day50': dict(zip(STATES, x50[:23])),
            'DA_day365': x365[18],
            'cumulative_loss_day35': x35[23:26], 'cumulative_loss_day50': x50[23:26]}
        doses_mgkg = [0.5, 1, 2, 5]
        vals = [day(simulate(q, days=35, dose=m / 5, y0=y0, points=71), 35) for m in doses_mgkg]
        drug[scen]['tier2_dose_ordering'] = {'dose_mg_per_kg': doses_mgkg, 'DA': [v[18] for v in vals],
                                            'CL_ratio': [cl_ratio(v) for v in vals], 'aSyn_olig': [v[6] for v in vals]}
        drug[scen]['tier2_monotonic'] = bool(np.all(np.diff([v[18] for v in vals]) > 0) and np.all(np.diff([cl_ratio(v) for v in vals]) > 0)
                                             and np.all(np.diff([v[6] for v in vals]) < 0))
    ablation = end(dict(disease, k_drug=0.), 35, dose=1, y0=y0)
    tlr2 = day(simulate(dict(disease, k_act_MG=0.), days=35, y0=y0, points=71), 35)
    impair = day(simulate(dict(disease, k_impair=0.), days=35, y0=y0, points=71), 35)
    tier1_keys = ['CL_ratio', 'aSyn_olig', 'CI_activity', 'mROS', 'mPTP_open', 'dpsi', 'ATP', 'MG_active', 'TNF', 'PINK1_act', 'mito_damage', 'DA_neuron']
    val = lambda y, k: cl_ratio(y) if k == 'CL_ratio' else y[INDEX[k]]
    tier1 = {k: [val(h35, k), val(d35, k)] for k in tier1_keys}
    shares = lambda c: (np.asarray(c) / np.sum(c)).tolist()
    result = {
        'k_drug': {'BJP': K_DRUG_BJP, 'A_mouse_anchored': SCENARIOS['A'], 'B_literal': SCENARIOS['B'],
                   'mouse_exposure_reference_true_SS_mean': disease['exposure_reference'], 'BJP_C_ref': C_REF_MOUSE_BJP,
                   'BJP_C_ref_as_fraction_of_true_SS_mean': C_REF_MOUSE_BJP / disease['exposure_reference'],
                   'max_abs_difference_A_vs_manuscript_normalisation_50d': equivalence},
        'table2': {
            'CI_activity_absolute': d35[8], 'CI_activity_relative_to_healthy': d35[8] / h35[8],
            'CL_native_drop': 1 - d35[0] / h35[0], 'mROS_ratio': d35[13] / h35[13],
            'DA_loss_day35': 1 - d35[18], 'motor_day35': d35[22],
            'TLR2_input_ablation_aSyn_reduction': 1 - tlr2[6] / d35[6],
            'clearance_impairment_ablation_aSyn_reduction': 1 - impair[6] / d35[6],
            'healthy_DA_day50': h50[18], 'healthy_DA_day35': h35[18],
            'k_drug_zero_ablation_DA_change': ablation[18] - d35[18]},
        'tier1_healthy_vs_disease_day35': tier1,
        'tier3': {'DA_checkpoints': {str(d): day(D, d)[18] for d in [30, 90, 180, 365]},
                  'motor_day365': D.y[22, -1], 'healthy_DA_day365': H.y[18, -1],
                  'aSyn_olig_disease_to_healthy_ratio_day35': d35[6] / h35[6]},
        'drug': drug,
        'death': {'disease_cumulative_day35': d35[23:26], 'disease_cumulative_day50': d50[23:26],
                  'disease_share_cumulative_day35': shares(d35[23:26]), 'disease_share_cumulative_day50': shares(d50[23:26]),
                  'disease_share_instantaneous_day35': shares(hazards(d35, disease))},
    }
    # Boundedness under the corrected state definitions.
    checks = {}
    for label, q, dose in [('healthy', healthy, 0), ('disease', disease, 0), ('low', A, 0.2), ('high', A, 2)]:
        s = simulate(q, days=365, dose=dose, y0=y0, points=366)
        checks[label] = {'min_state': float(s.y[:23].min()), 'max_fraction': float(s.y[FRACTIONS].max()),
                         'max_pool_index': float(s.y[[0, 1, 2, 5, 6, 7, 13]].max()), 'monotone_DA': bool(np.all(np.diff(s.y[18]) <= 1e-12))}
    result['boundedness'] = checks
    # Time courses for Figure 3 and phase projection (Figure 4), dose-response for Figure S3A.
    series = [('healthy', healthy, 0), ('disease', disease, 0), ('A_1mgkg', A, 0.2), ('A_5mgkg', A, 1)]
    rows = []
    for label, q, dose in series:
        s = simulate(q, days=50, dose=dose, y0=y0, points=101)
        rows += [[label, t / 24, *yy] for t, yy in zip(s.t, s.y.T)]
    write_csv('mouse_timeseries.csv', ['scenario', 'day', *STATES, 'D_mito', 'D_inflam', 'D_proteo'], rows)
    dr = []
    for scen in 'AB':
        q = with_kdrug(disease, scen)
        for m in [0, 0.25, 0.5, 1, 2, 3, 5, 7.5, 10]:
            y = end(q, 50, dose=m / 5, y0=y0)
            dr.append([scen, m, y[18], cl_ratio(y), y[6], y[22]])
    write_csv('mouse_dose_response_day50.csv', ['scenario', 'dose_mg_per_kg', 'DA', 'CL_ratio', 'aSyn_olig', 'motor'], dr)
    save_json('mouse.json', result)


# ---------------------------------------------------------------- PK (Table 4)
def pk_table():
    out = {}; curves = []
    for species in ['mouse', 'human']:
        p = parameters(species); A = pk_matrix(p); E = expm(A * 24); b = np.array([1., 0, 0, 0, 0])
        trough = np.linalg.solve(np.eye(5) - E, b); trough = E @ trough
        x = np.zeros(5); reached = {}
        for d in range(1, 1461):
            x = E @ (x + b); frac = x[4] / trough[4]; curves.append([species, d, frac])
            for level in (0.5, 0.9, 0.95):
                if frac >= level and level not in reached: reached[level] = d
            if d == 120: day120 = frac
        eig = np.linalg.eigvals(A[1:, 1:]).real
        out[species] = {'k_a': p['k_a'], 'ke_plasma': p['ke_plasma'], 'k_12': p['k_12'], 'k_21': p['k_21'],
                        'k_brain_in': p['k_brain_in'], 'k_brain_out': p['k_brain_out'], 'k_mito_in': p['k_mito_in'], 'k_mito_out': p['k_mito_out'],
                        'central_elimination_half_time_h': np.log(2) / p['ke_plasma'],
                        'slowest_mode_half_life_days': np.log(2) / -eig.max() / 24,
                        'mean_SS_mito_amount_per_unit_daily_dose': p['exposure_reference'],
                        'SS_trough_mito_amount': trough[4], 'BJP_C_ref': p['C_ref'],
                        'day120_fraction_of_SS_trough': day120, 'days_to_fraction_of_SS': {str(k): v for k, v in reached.items()}}
    out['allometric_k_brain_out_from_mouse'] = 0.016 * (70 / 0.025) ** -0.25
    out['HED_mg_70kg_FDA_Km'] = 5 * 3 / 37 * 70
    save_json('pk.json', out)
    write_csv('pk_accumulation.csv', ['species', 'day', 'fraction_of_SS_trough'], curves)


# ---------------------------------------------------------------- human nominal (Fig 5, ST4)
def human_nominal():
    p = parameters('human'); y0 = initial_state(p)
    ph = parameters('human', healthy=True); yh = initial_state(ph)
    out = {'untreated_crossing_years': {str(t): crossing(p, y0, t, 40)[0] for t in (0.9, 0.7, 0.5, 0.4)},
           'healthy_DA_20y': end(ph, 20 * 365.25, y0=yh)[18]}
    years = [1, 2, 5, 10, 15, 20]
    traj_rows = []
    grid = np.linspace(0, 20 * 365.25 * 24, 20 * 52 + 1)
    for label, q, dose in [('healthy', ph, 0), ('disease', p, 0), ('A_30mg', with_kdrug(p, 'A'), 1), ('B_30mg', with_kdrug(p, 'B'), 1)]:
        s = simulate(q, days=20 * 365.25, dose=dose, y0=yh if label == 'healthy' else y0, points=len(grid))
        traj_rows += [[label, t / YEAR_H, yy[18], yy[22]] for t, yy in zip(s.t, s.y.T)]
        if label != 'healthy':
            out[label] = {'DA_at_years': {str(yr): float(np.interp(yr * YEAR_H, s.t, s.y[18])) for yr in years}}
    write_csv('human_20y_trajectories.csv', ['scenario', 'year', 'DA', 'motor'], traj_rows)
    for scen in 'AB':
        q = with_kdrug(p, scen)
        onset = {str(t): crossing(q, y0, t, 45, dose=1)[0] for t in (0.7, 0.5, 0.4)}
        out[f'{scen}_30mg']['treated_from_time0_crossing_years'] = onset
        out[f'{scen}_30mg']['delay_years'] = {k: (v - out['untreated_crossing_years'][k]) if v else None for k, v in onset.items()}
        out[f'{scen}_30mg']['DA_benefit_at_years'] = {k: v - out['disease']['DA_at_years'][k] for k, v in out[f'{scen}_30mg']['DA_at_years'].items()}
    # ST4: prodromal onset -> diagnostic delay -> 18 months treatment, nominal patient.
    _, onset_state = crossing(p, y0)
    st4 = []
    for delay_months in [0, 6, 12, 24]:
        start = end(p, delay_months * 365.25 / 12, y0=onset_state)[:23] if delay_months else onset_state
        placebo = end(p, 547, y0=start)
        for scen in 'AB':
            for mg in HUMAN_DOSES_MG:
                treated = end(with_kdrug(p, scen), 547, dose=mg / REF_DOSE_MG, y0=start)
                st4.append([scen, delay_months, mg, start[18], treated[18], placebo[18], treated[18] - placebo[18],
                            bool(treated[18] <= start[18] + 1e-12)])
    for row in st4:
        ref = [r for r in st4 if r[0] == row[0] and r[2] == row[2] and r[1] == 0][0][6]
        row.append(row[6] / ref - 1)
    write_csv('ST4_phased.csv', ['scenario', 'delay_months', 'dose_mg', 'DA_at_start', 'DA_after_18mo_treated', 'DA_after_18mo_placebo',
                                'delta_DA', 'treated_not_above_start', 'relative_change_vs_0mo_delay'], st4)
    save_json('human_nominal.json', out)


# ---------------------------------------------------------------- virtual trial (Table 5, ST2, ST3, Fig 6/7)
LOGNORMAL = {'k_ROS_basal': .25, 'K_mPTP': .25, 'k_SOD2': .25, 'k_agg': .30, 'k_shunt': .20, 'k_clear': .25,
             'k_death': .30, 'K_death': .25, 'k_death_inflam': .50, 'k_death_aSyn': .50,
             'k_a': .30, 'ke_plasma': .30, 'k_brain_in': .40}
TRUNC_NORMAL = {'CI_max': (.90, .03, .80, .96), 'alpha_clear': (.35, .05, .20, .55)}


def sample_population(n, rng):
    nominal = parameters('human'); rows = []
    for i in range(n):
        row = {}
        for k, cv in LOGNORMAL.items():
            s = np.sqrt(np.log(1 + cv ** 2)); row[k] = float(np.exp(rng.normal(np.log(nominal[k]) - s * s / 2, s)))
        for k, (m, sd, lo, hi) in TRUNC_NORMAL.items():
            while True:
                v = rng.normal(m, sd)
                if lo <= v <= hi: break
            row[k] = float(v)
        rows.append((i, row))
    return rows


def screen_job(task):
    ident, overrides = task
    try:
        p = parameters('human', overrides=overrides); y0 = initial_state(p)
        years, state = crossing(p, y0)
        if years is None: return {'id': ident, 'status': 'no_onset_within_30y', 'parameters': overrides}
        return {'id': ident, 'status': 'eligible', 'onset_years': years, 'onset_state': state.tolist(), 'parameters': overrides}
    except Exception as e:
        return {'id': ident, 'status': 'numerical_failure_screening', 'reason': str(e), 'parameters': overrides}


def trial_job(record):
    try:
        p = parameters('human', overrides=record['parameters'])
        diag = end(p, 365.25, y0=np.array(record['onset_state']))[:23]
        arms = {'placebo': end(p, 547, y0=diag)}
        for scen in 'AB':
            for mg in HUMAN_DOSES_MG:
                arms[f'{scen}_{mg}'] = end(with_kdrug(p, scen), 547, dose=mg / REF_DOSE_MG, y0=diag)
        return {**{k: v for k, v in record.items() if k != 'onset_state'}, 'status': 'included', 'diagnosis_DA': diag[18],
                'DA': {k: v[18] for k, v in arms.items()}, 'motor': {k: v[22] for k, v in arms.items()},
                'cumulative_loss': {k: v[23:26].tolist() for k, v in arms.items()}}
    except Exception as e:
        return {**{k: v for k, v in record.items() if k != 'onset_state'}, 'status': 'numerical_failure_treatment', 'reason': str(e)}


def cliffs_delta(x, y):
    """P(x > y) - P(x < y) using mid-ranks."""
    r = rankdata(np.r_[x, y]); n1, n2 = len(x), len(y)
    u = r[:n1].sum() - n1 * (n1 + 1) / 2
    return 2 * u / (n1 * n2) - 1


def power(n, d):
    df = 2 * n - 2; crit = tdist.ppf(.975, df); nc = abs(d) * np.sqrt(n / 2)
    return nct.sf(crit, df, nc) + nct.cdf(-crit, df, nc)


def n_for_power(d, target=.8):
    if d <= 0: return None
    lo, hi = 2, 2
    while power(hi, d) < target: hi *= 2
    while lo < hi:
        mid = (lo + hi) // 2
        if power(mid, d) >= target: hi = mid
        else: lo = mid + 1
    return lo


def endpoint_stats(placebo, treated, sign, threshold, rng, boots):
    """sign=+1 when larger is better (DA), -1 when smaller is better (motor). Paired virtual patients."""
    benefit = sign * (treated - placebo)
    pooled = np.sqrt((np.var(placebo, ddof=1) + np.var(treated, ddof=1)) / 2)
    def stats(ix):
        a, b = placebo[ix], treated[ix]; ben = sign * (b - a)
        sd = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2)
        return [ben.mean() / sd, cliffs_delta(sign * b, sign * a), np.mean(ben > threshold), ben.mean()]
    n = len(placebo); full = stats(np.arange(n))
    bs = np.array([stats(rng.integers(0, n, n)) for _ in range(boots)])
    lo, hi = np.quantile(bs, [.025, .975], axis=0)
    return {'mean_difference': float(np.mean(treated - placebo)), 'mean_benefit': float(benefit.mean()),
            'd_unpaired': full[0], 'cliffs_delta': full[1], 'response_rate': full[2],
            'd_unpaired_95': [lo[0], hi[0]], 'cliffs_delta_95': [lo[1], hi[1]], 'response_rate_95': [lo[2], hi[2]],
            'mean_benefit_95': [lo[3], hi[3]], 'n_per_arm_80pct_power': n_for_power(full[0]),
            'pooled_sd': pooled, 'median_benefit': float(np.median(benefit)),
            'benefit_quartile_share_top25': float(np.sort(benefit)[-n // 4:].sum() / benefit.sum()) if benefit.sum() > 0 else None}


def trial(pool, n_target, n_candidates, boots):
    rng = np.random.default_rng(SEED)
    candidates = sample_population(n_candidates, rng)
    screened = sorted(pool.map(screen_job, candidates, chunksize=8), key=lambda r: r['id'])
    eligible = [r for r in screened if r['status'] == 'eligible']
    if len(eligible) < n_target: raise RuntimeError(f'Only {len(eligible)} eligible of {n_candidates}; increase candidates')
    last_id = eligible[n_target - 1]['id']
    screened_used = [r for r in screened if r['id'] <= last_id]
    enrolled = eligible[:n_target]
    results = sorted(pool.map(trial_job, enrolled, chunksize=2), key=lambda r: r['id'])
    included = [r for r in results if r['status'] == 'included']
    records = [r for r in screened_used if r['status'] != 'eligible'] + results
    save_json('vct_records.json', sorted([{k: v for k, v in r.items() if k != 'onset_state'} for r in records], key=lambda r: r['id']))
    flow = {'candidates_screened': len(screened_used), 'eligible_enrolled': len(enrolled), 'included_analysed': len(included),
            'status_counts': {s: sum(r['status'] == s for r in records) for s in sorted(set(r['status'] for r in records))}}
    DA = {k: np.array([r['DA'][k] for r in included]) for k in included[0]['DA']}
    MO = {k: np.array([r['motor'][k] for r in included]) for k in included[0]['motor']}
    CUM = {k: np.array([r['cumulative_loss'][k] for r in included]) for k in included[0]['cumulative_loss']}
    brng = np.random.default_rng(SEED + 1)
    arms = {}
    for scen in 'AB':
        for mg in HUMAN_DOSES_MG:
            key = f'{scen}_{mg}'
            rescue = 1 - CUM[key].sum(axis=0) / CUM['placebo'].sum(axis=0)
            arms[key] = {'scenario': scen, 'dose_mg': mg,
                         'DA': endpoint_stats(DA['placebo'], DA[key], +1, .05, brng, boots),
                         'motor': endpoint_stats(MO['placebo'], MO[key], -1, .03, brng, boots),
                         'rescue_fraction_cumulative_loss': {'mitochondrial': rescue[0], 'inflammation': rescue[1], 'proteotoxic': rescue[2],
                                                             'total': 1 - CUM[key].sum() / CUM['placebo'].sum()}}
    diag = np.array([r['diagnosis_DA'] for r in included]); onset = np.array([r['onset_years'] for r in included])
    summary = {'flow': flow, 'bootstrap_resamples': boots,
               'diagnosis_DA_range': [diag.min(), diag.max()], 'diagnosis_DA_median': float(np.median(diag)),
               'onset_years_median_IQR': np.quantile(onset, [.5, .25, .75]),
               'placebo_DA_18mo_median_range': [float(np.median(DA['placebo'])), DA['placebo'].min(), DA['placebo'].max()],
               'placebo_cumulative_share': (CUM['placebo'].sum(axis=0) / CUM['placebo'].sum()).tolist(),
               'arms': arms}
    save_json('vct_summary.json', summary)
    rows = [[r['id'], r['onset_years'], r['diagnosis_DA'], arm, r['DA'][arm], r['motor'][arm], *r['cumulative_loss'][arm]]
            for r in included for arm in r['DA']]
    write_csv('vct_endpoints.csv', ['id', 'onset_years', 'diagnosis_DA', 'arm', 'DA_18mo', 'motor_18mo', 'loss_mito', 'loss_inflam', 'loss_proteo'], rows)
    return summary


# ---------------------------------------------------------------- translational envelope (ST3B, Fig 7A)
def envelope_job(task):
    ci, alpha, scale = task
    base = base_parameters()
    try:
        p = parameters('human', overrides={'CI_max': ci, 'alpha_clear': alpha, 'k_death': base['k_death'] * scale,
                                           'k_death_inflam': base['k_death_inflam'] * scale, 'k_death_aSyn': base['k_death_aSyn'] * scale})
        years, state = crossing(p, initial_state(p))
        if years is None: return {'CI_max': ci, 'alpha_clear': alpha, 'death_scale': scale, 'status': 'no_onset_within_30y'}
        diag = end(p, 365.25, y0=state)[:23]; placebo = end(p, 547, y0=diag)[18]
        saved = {s: end(with_kdrug(p, s), 547, dose=1, y0=diag)[18] - placebo for s in 'AB'}
        return {'CI_max': ci, 'alpha_clear': alpha, 'death_scale': scale, 'status': 'included', 'onset_years': years,
                'DA_saved_A': saved['A'], 'DA_saved_B': saved['B']}
    except Exception as e:
        return {'CI_max': ci, 'alpha_clear': alpha, 'death_scale': scale, 'status': 'numerical_failure', 'reason': str(e)}


def envelope(pool, sd_A, sd_B):
    grid = [(c, a, s) for c in np.linspace(.86, .94, 5) for a in np.linspace(.28, .42, 5) for s in np.linspace(.04, .08, 5)]
    res = list(pool.map(envelope_job, grid, chunksize=2))
    ok = [r for r in res if r['status'] == 'included']
    out = {'grid_points': len(grid), 'included': len(ok), 'status_counts': {s: sum(r['status'] == s for r in res) for s in set(r['status'] for r in res)}}
    for scen, sd in [('A', sd_A), ('B', sd_B)]:
        v = np.array([r[f'DA_saved_{scen}'] for r in ok])
        out[scen] = {'DA_saved_min_max': [v.min(), v.max()], 'DA_saved_2.5_97.5': np.quantile(v, [.025, .975]),
                     'approx_d_min_max': [v.min() / sd, v.max() / sd], 'pooled_sd_used': sd}
        for key in ['CI_max', 'alpha_clear', 'death_scale']:
            means = {}
            for r in ok: means.setdefault(r[key], []).append(r[f'DA_saved_{scen}'])
            m = [np.mean(x) for x in means.values()]
            out[scen][f'fold_range_of_mean_DA_saved_over_{key}'] = max(m) / min(m)
    save_json('envelope.json', out)
    write_csv('envelope.csv', ['CI_max', 'alpha_clear', 'death_scale', 'status', 'onset_years', 'DA_saved_A', 'DA_saved_B'],
              [[r['CI_max'], r['alpha_clear'], r['death_scale'], r['status'], r.get('onset_years'), r.get('DA_saved_A'), r.get('DA_saved_B')] for r in res])


# ---------------------------------------------------------------- Sobol at the manuscript's +/-50% range (Table 3)
SOBOL = ['k_ROS_basal', 'K_mPTP', 'CI_max', 'k_SOD2', 'k_agg', 'k_shunt', 'alpha_clear', 'k_CL_loss', 'k_CI_repair', 'k_red',
         'k_death', 'K_motor', 'k_clear']


def sobol_job(multipliers):
    p = parameters(); q = {k: p[k] * m for k, m in zip(SOBOL, multipliers)}
    q['CI_max'] = min(q['CI_max'], 1.)
    try:
        y = end(parameters(overrides=q), 35)
        return [y[18], y[22]]
    except Exception:
        return [np.nan, np.nan]


def sobol(pool, nbase):
    k = len(SOBOL); u = qmc.Sobol(d=2 * k, scramble=True, seed=SEED).random_base2(int(np.log2(nbase)))
    A, B = .5 + u[:, :k], .5 + u[:, k:]
    mats = [A, B] + [np.where(np.arange(k) == j, B, A) for j in range(k)]
    out = np.array(list(pool.map(sobol_job, [row for m in mats for row in m], chunksize=16))).reshape(len(mats), nbase, 2)
    good = ~np.isnan(out).any(axis=(0, 2)); out = out[:, good]; n = good.sum()
    ya, yb = out[0], out[1]
    def indices(ix):
        a, b = ya[ix], yb[ix]; var = np.var(np.r_[a, b], axis=0, ddof=1)
        st = np.array([np.mean((a - m[ix]) ** 2, axis=0) / (2 * var) for m in out[2:]])
        s1 = np.array([1 - np.mean((b - m[ix]) ** 2, axis=0) / (2 * var) for m in out[2:]])
        return s1, st
    s1, st = indices(np.arange(n)); rng = np.random.default_rng(SEED + 2)
    boot = [indices(rng.integers(0, n, n)) for _ in range(500)]
    s1ci = np.quantile([b[0] for b in boot], [.025, .975], axis=0); stci = np.quantile([b[1] for b in boot], [.025, .975], axis=0)
    save_json('sobol_pm50.json', {'parameters': SOBOL, 'outputs': ['DA_day35', 'motor_day35'], 'n_base': nbase, 'valid_rows': int(n),
                                   'evaluations': len(mats) * nbase, 'bounds': 'Independent uniform 50-150% of nominal (manuscript OAT range); CI_max capped at 1',
                                   'S1': s1, 'ST': st, 'S1_95': s1ci, 'ST_95': stci})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument('--patients', type=int, default=250)
    ap.add_argument('--candidates', type=int, default=1200)
    ap.add_argument('--boots', type=int, default=10000)
    ap.add_argument('--sobol-base', type=int, default=256)
    ap.add_argument('--skip-sobol', action='store_true')
    a = ap.parse_args()
    OUT.mkdir(exist_ok=True); t0 = time.perf_counter(); timings = {}
    def step(name, fn, *args):
        s = time.perf_counter(); print(name, flush=True); r = fn(*args); timings[name] = time.perf_counter() - s; return r
    step('mouse', mouse); step('pk', pk_table); step('human_nominal', human_nominal)
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        summary = step('virtual_trial', trial, pool, a.patients, a.candidates, a.boots)
        step('envelope', envelope, pool, summary['arms']['A_30']['DA']['pooled_sd'], summary['arms']['B_30']['DA']['pooled_sd'])
        if not a.skip_sobol: step('sobol', sobol, pool, a.sobol_base)
    save_json('run_metadata.json', {'seed': SEED, 'workers': a.workers, 'patients': a.patients, 'candidates_generated': a.candidates,
                                    'bootstrap': a.boots, 'sobol_base': a.sobol_base, 'elapsed_seconds': time.perf_counter() - t0,
                                    'step_seconds': timings, 'python': platform.python_version(), 'numpy': np.__version__,
                                    'scipy': scipy.__version__, 'solver': 'LSODA rtol 1e-8 atol 1e-10', 'k_drug': SCENARIOS})
    print('done', time.perf_counter() - t0, flush=True)


if __name__ == '__main__':
    main()
