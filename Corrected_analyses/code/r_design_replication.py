"""Replicate the GitHub R virtual-trial design with the Python engine, then change one design choice at a time.

Designs (all 250 virtual patients sampled as in run_arm.R: set of 13 log-normal parameters, CI_max and
alpha_clear from clipped normals):
  R      : start from the R integrated_init() healthy state, dose from day 0 for 547 days, human exposure
           normalised by C_ref_human = day-120 trough of the 30 mg regimen (as in run_arm.R).
  R+enrol: same normalisation, but each patient first progresses untreated to DA < 0.70 plus 12 months.
Results are compared with scenario A/B of bjp_corrected.py (steady-state normalisation plus enrolment).
"""
import os
for _n in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[_n] = '1'
import json
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from model import ROOT, PK, parameters, simulate
from bjp_corrected import crossing, endpoint_stats, SCENARIOS

R_INIT = np.array([.912, .051, .002, .100, .900, .950, .051, .003, .950, .970, .960, .954, .002, .100, .833, .667,
                   .200, .046, 1.0, .134, .211, .167, .057])
DOSE_NORM_30MG = (30 / 70) / 5          # run_arm.R: dose_mg_kg / dose_ref
LOGN = {'k_ROS_basal': .25, 'K_mPTP': .25, 'k_SOD2': .25, 'k_agg': .30, 'k_shunt': .20, 'k_clear': .25, 'k_death': .30,
        'K_death': .25, 'k_death_inflam': .50, 'k_death_aSyn': .50, 'k_a': .30, 'ke_plasma': .30, 'k_brain_in': .40}


def r_c_ref():
    """Trough mitochondrial amount on day 120 of 30 mg daily with nominal human PK, in R dose units."""
    pk = PK(parameters('human'), DOSE_NORM_30MG)
    return min(pk.mito(t) for t in np.arange(119 * 24, 120 * 24 + 1, 1.0))


def sample(n, seed=42):
    rng = np.random.default_rng(seed); nom = parameters('human'); rows = []
    for i in range(n):
        row = {}
        for k, cv in LOGN.items():
            s = np.sqrt(np.log(1 + cv ** 2)); row[k] = float(np.exp(rng.normal(np.log(nom[k]) - s * s / 2, s)))
        row['CI_max'] = float(np.clip(rng.normal(.90, .03), .80, .96))
        row['alpha_clear'] = float(np.clip(rng.normal(.35, .05), .20, .55))
        rows.append((i, row))
    return rows


def job(task):
    ident, over = task
    out = {'id': ident}
    try:
        p = parameters('human', overrides=over)
        r = dict(p, exposure_reference=r_c_ref() / DOSE_NORM_30MG, k_drug=0.1)   # relative-dose units
        for label, y0 in [('R', R_INIT)]:
            out[label] = {str(d): simulate(r, days=547, dose=d, y0=y0, points=2).y[[18, 22], -1].tolist() for d in (0, 1 / 3, 1, 2)}
        from model import initial_state
        years, state = crossing(p, initial_state(p))
        if years is None:
            out['R_enrol'] = None
        else:
            diag = simulate(p, days=365.25, y0=state, points=2).y[:23, -1]
            out['R_enrol'] = {str(d): simulate(r, days=547, dose=d, y0=diag, points=2).y[[18, 22], -1].tolist() for d in (0, 1 / 3, 1, 2)}
        out['placebo_R_DA'] = out['R']['0'][0]
    except Exception as e:
        out['error'] = str(e)
    return out


def summarise(records, key):
    rows = [r[key] for r in records if r.get(key)]
    rng = np.random.default_rng(7); res = {'n': len(rows)}
    pl = np.array([r['0'] for r in rows])
    res['placebo_DA_median_min_max'] = [float(np.median(pl[:, 0])), float(pl[:, 0].min()), float(pl[:, 0].max())]
    for d, mg in [('0.3333333333333333', 10), ('1', 30), ('2', 60)]:
        tr = np.array([r[d] for r in rows])
        res[f'{mg}mg'] = {'DA': endpoint_stats(pl[:, 0], tr[:, 0], +1, .05, rng, 2000),
                          'motor': endpoint_stats(pl[:, 1], tr[:, 1], -1, .03, rng, 2000)}
    return res


if __name__ == '__main__':
    with ProcessPoolExecutor(max_workers=22) as pool:
        records = list(pool.map(job, sample(250), chunksize=2))
    out = {'C_ref_human_R_units': r_c_ref(), 'effective_kdrug_at_SS': 0.1 * parameters('human')['exposure_reference'] * DOSE_NORM_30MG / r_c_ref(),
           'scenario_A_kdrug': SCENARIOS['A'], 'errors': sum('error' in r for r in records),
           'R_design': summarise(records, 'R'), 'R_normalisation_with_enrolment': summarise(records, 'R_enrol')}
    (ROOT / 'results' / 'r_design_replication.json').write_text(json.dumps(out, indent=2, default=float), encoding='utf-8')
    for key in ['R_design', 'R_normalisation_with_enrolment']:
        s = out[key]
        print(key, 'n', s['n'], 'placebo DA median/min/max', [round(x, 3) for x in s['placebo_DA_median_min_max']])
        for mg in ('10mg', '30mg', '60mg'):
            a, m = s[mg]['DA'], s[mg]['motor']
            print(f"  {mg}: dDA={a['mean_difference']:+.4f} d={a['d_unpaired']:.3f} resp={100*a['response_rate']:.1f}% N80={a['n_per_arm_80pct_power']} | motor d={m['d_unpaired']:.3f} resp={100*m['response_rate']:.1f}%")
    print('C_ref_human', out['C_ref_human_R_units'], 'effective k_drug*u at SS', out['effective_kdrug_at_SS'], 'errors', out['errors'])
