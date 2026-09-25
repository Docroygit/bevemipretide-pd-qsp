"""Checks that the corrected analyses satisfy the invariants the BJP tables violated."""
import csv, json
import numpy as np
from model import ROOT, parameters, initial_state, simulate
from bjp_corrected import SCENARIOS, C_REF_MOUSE_BJP, DA_THRESHOLD

R = ROOT / 'results'


def load(name): return json.loads((R / name).read_text(encoding='utf-8'))


def test_scenario_A_reproduces_manuscript_mouse_normalisation():
    p = parameters(); y0 = initial_state()
    src = simulate(p, days=35, dose=1, y0=y0, points=71, reference_mode='manuscript')
    anc = simulate(dict(p, k_drug=SCENARIOS['A']), days=35, dose=1, y0=y0, points=71)
    assert np.max(np.abs(src.y - anc.y)) < 1e-7
    assert abs(SCENARIOS['A'] - 0.1 * p['exposure_reference'] / C_REF_MOUSE_BJP) < 1e-12


def test_ST4_treatment_never_raises_survival():
    with (R / 'ST4_phased.csv').open(encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    assert rows and all(float(r['DA_after_18mo_treated']) <= float(r['DA_at_start']) + 1e-12 for r in rows)
    assert all(float(r['DA_after_18mo_treated']) >= float(r['DA_after_18mo_placebo']) for r in rows)


def test_trial_enrols_only_symptomatic_patients_and_accounts_for_all():
    s = load('vct_summary.json'); records = load('vct_records.json')
    assert s['diagnosis_DA_range'][1] < DA_THRESHOLD
    assert s['placebo_DA_18mo_median_range'][2] < DA_THRESHOLD
    assert len(records) == s['flow']['candidates_screened']
    assert s['flow']['included_analysed'] == s['flow']['eligible_enrolled']


def test_cumulative_loss_accounting_is_complete():
    p = parameters(); y0 = initial_state()
    s = simulate(dict(p, k_drug=SCENARIOS['A']), days=50, dose=1, y0=y0, points=11)
    assert np.max(np.abs(s.y[18] + s.y[23:26].sum(axis=0) - 1)) < 1e-7
    e = list(csv.DictReader((R / 'vct_endpoints.csv').open(encoding='utf-8')))
    assert all(float(r['DA_18mo']) <= float(r['diagnosis_DA']) + 1e-12 for r in e)
