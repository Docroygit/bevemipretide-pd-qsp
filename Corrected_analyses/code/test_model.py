"""Numerical verification, not biological validation."""
import numpy as np
from scipy.linalg import expm
from model import *

def test_parameter_count():
    assert len(base_parameters())==89 # 90 table records, one derived dose expression.

def test_pk_mass_balance():
    p=parameters(); a=np.arange(1,6.)
    assert abs(np.sum(pk_matrix(p)@a)+p['ke_plasma']*a[1])<1e-12

def test_exact_pk_against_dose_recurrence():
    p=parameters(); pk=PK(p); a=np.zeros(5); b=np.array([1,0,0,0,0]); E=expm(pk_matrix(p)*24)
    for day in range(120):
        a+=b
        assert np.allclose(pk.amount(day*24+3),expm(pk_matrix(p)*3)@a,rtol=1e-9,atol=1e-10)
        a=E@a

def test_pk_periodic_reference():
    p=parameters(); A=pk_matrix(p); E=expm(A*24); b=np.array([1,0,0,0,0]); post=np.linalg.solve(np.eye(5)-E,b)
    mean=np.linalg.solve(A,(E-np.eye(5))@post)/24
    assert abs(mean[4]-pk_reference(p))<1e-9

def test_burnin_and_conservation():
    p=parameters(healthy=True); y=initial_state(p); d=rhs(0,y,p,freeze_neurons=True)
    assert np.max(abs(d))<1e-7
    c,o,e=y[:3]
    assert abs(sum(d[:3])-(p['k_syn_CL']*(1-c/p['K_CL'])-p['k_degrad_ox']*o-p['k_clear_ext_CL']*e))<1e-12

def test_neuron_survival_and_cumulative_hazards():
    s=simulate(days=50,dose=1)
    assert np.max(np.diff(s.y[18]))<=1e-9
    assert np.max(abs(s.y[18]+s.y[23:26].sum(axis=0)-1))<1e-7

def test_solver_and_tolerance_agreement():
    a=simulate(days=10,dose=1,points=21)
    b=simulate(days=10,dose=1,points=21,method='Radau',rtol=1e-9,atol=1e-11,max_step=6)
    assert np.max(abs(a.y-b.y))<2e-6

def test_drug_ablation_and_no_dose_equivalence():
    p=parameters(overrides={'k_drug':0.})
    a=simulate(p,days=10,dose=0,points=2); b=simulate(p,days=10,dose=1,points=2)
    assert np.max(abs(a.y-b.y))<1e-10

def test_initial_condition_continuity():
    p=parameters('human'); first=simulate(p,days=400,points=2)
    second=simulate(p,days=100,dose=1,y0=first.y[:,-1],points=2)
    assert np.allclose(first.y[:23,-1],second.y[:23,0],rtol=0,atol=1e-14)
    assert second.y[18,-1]<=first.y[18,-1]

def test_exposure_reference_not_subject_normalized():
    p=parameters('human');q=parameters('human',overrides={'k_brain_in':.03})
    assert p['exposure_reference']==q['exposure_reference']
    assert pk_reference(q)>pk_reference(p)

def test_compiled_kernel_matches_reference():
    from fast_kernel import make_rhs
    p=parameters();rng=np.random.default_rng(22)
    for mechanism in ['repair','protection']:
        pk=PK(p);fun=make_rhs(p,pk,p['exposure_reference'],False,mechanism)
        for i in range(10):
            y=rng.uniform(.01,.95,size=26);t=13.5+24*i
            assert np.allclose(fun(t,y),rhs(t,y,p,pk,p['exposure_reference'],False,mechanism),rtol=1e-13,atol=1e-13)
