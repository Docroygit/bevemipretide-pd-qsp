"""Independent reconstruction of the manuscript's 28-state QSP equations.

Five linear PK amount-equivalent states are solved analytically and force 23 PD
states. This is mathematically equivalent to the triangular 28-state system.
No concentrations, physical dose conversion, clinical calibration, or original
code recovery is claimed. All time units are hours. No state clamping is used.
"""
from pathlib import Path
from functools import lru_cache
import csv
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
STATES = ['CL_n','CL_ox','CL_ext','ALCAT1','TAZ','aSyn_mono','aSyn_olig','aSyn_ext',
          'CI_activity','SC_integrity','dpsi','ATP','mPTP_open','mROS','SOD2_act',
          'GPx4_act','PINK1_act','mito_damage','DA_neuron','MG_active','TNF','IL1b','Motor_score']
INDEX = {s:i for i,s in enumerate(STATES)}
PK_STATES = ['Depot','A_plasma','A_periph','A_brain','A_mito']
FRACTIONS = [3,4,8,9,10,11,12,14,15,16,17,18,19,20,21,22]

@lru_cache(maxsize=1)
def base_parameters():
    with (ROOT/'data'/'parameters.csv').open(encoding='utf-8-sig') as f:
        return {r['parameter']:float(r['value']) for r in csv.DictReader(f) if r['value']}

def parameters(species='mouse', healthy=False, overrides=None):
    p = base_parameters().copy()
    if species == 'human':
        p.update(k_a=.4,ke_plasma=.032,k_12=.041,k_21=.027,k_brain_in=.015,k_brain_out=.010,
                 C_ref=2.15,CI_max=.90,alpha_clear=.35)
        for k in ['k_death','k_death_inflam','k_death_aSyn']: p[k] *= .06
    elif species != 'mouse': raise ValueError(species)
    if healthy: p.update(CI_max=1.,alpha_clear=1.)
    # Fixed species reference set BEFORE subject-specific perturbations.
    # Otherwise normalizing each subject by its own exposure erases PK variability.
    p['exposure_reference']=pk_reference(p)
    if overrides: p.update(overrides)
    return p

def pk_matrix(p):
    ka,ke,k12,k21,bi,bo,mi,mo = [p[k] for k in ['k_a','ke_plasma','k_12','k_21','k_brain_in','k_brain_out','k_mito_in','k_mito_out']]
    return np.array([[-ka,0,0,0,0],[ka,-ke-k12-bi,k21,bo,0],[0,k12,-k21,0,0],
                     [0,bi,0,-bo-mi,mo],[0,0,0,mi,-mo]],dtype=float)

def pk_reference(p):
    """Exact mean A_mito at periodic steady state for one unit every 24 h."""
    b=np.zeros(5); b[0]=1/24
    return float(np.linalg.solve(-pk_matrix(p),b)[4])

class PK:
    def __init__(self,p,dose=1.,interval=24.):
        self.matrix=pk_matrix(p); self.dose=float(dose); self.interval=float(interval)
        self.lam,self.vec=np.linalg.eig(self.matrix)
        self.coef=np.linalg.solve(self.vec,np.array([dose,0,0,0,0]))
        if np.max(np.abs(self.lam.imag))>1e-10: raise ValueError('Unexpected complex PK modes')
        self.lam=self.lam.real; self.vec=self.vec.real; self.coef=self.coef.real
    def amount(self,t):
        if t<0 or self.dose==0: return np.zeros(5)
        n=int(np.floor(t/self.interval+1e-12))+1
        tau=max(0.,t-(n-1)*self.interval)
        geo=np.expm1(self.lam*n*self.interval)/np.expm1(self.lam*self.interval)
        return self.vec @ (self.coef*np.exp(self.lam*tau)*geo)
    def mito(self,t):
        if self.dose==0 or t<0: return 0.
        n=int(np.floor(t/self.interval+1e-12))+1
        tau=max(0.,t-(n-1)*self.interval)
        geo=np.expm1(self.lam*n*self.interval)/np.expm1(self.lam*self.interval)
        return float(np.dot(self.vec[4],self.coef*np.exp(self.lam*tau)*geo))

def hill(x,k,n):
    # Negative states are errors, not hidden by clipping. Tolerate roundoff near 0.
    # Permit tiny negative Jacobian probes without projecting any state.
    if x< -1e-5: raise FloatingPointError(f'Negative Hill input {x}')
    if x<=0: return 0.
    z=(x/k)**n
    return z/(1+z)

def hazards(y,p):
    return np.array([p['k_death']*hill(y[17],p['K_death'],p['n_death']),
        p['k_death_inflam']*hill((y[20]+y[21])/2,p['K_inflam_death'],p['n_inflam_death']),
        p['k_death_aSyn']*hill(y[6],p['K_aSyn_death'],p['n_aSyn_death'])])

def rhs(t,y,p,pk=None,reference=1.,freeze_neurons=False,mechanism='repair'):
    c,o,e,a,z,m,q,x,ci,sc,v,atp,pore,r,s,g,pink,damage,neur,mg,tnf,il,motor=y[:23]
    ratio=c/(c+o) if c+o>0 else 0.
    drug=pk.mito(t)/reference if pk else 0.
    oxid=p['k_ox']*r*c*(1+p['k_ALCAT']*a)+p['k_aSyn_ox']*q*c
    recycle=p['k_red']*z*o
    if mechanism=='repair': recycle += p['k_drug']*drug*o
    elif mechanism=='protection': oxid /= (1+drug) # New dimensionless potency assumption, not equipotent to repair.
    else: raise ValueError(mechanism)
    external=p['k_ext']*c*r
    stress=1+p['k_ROS_agg']*r+p['k_CL_loss']*(1-ratio)+p['k_seed_CL']*e
    aggregation=p['k_agg']*stress*m
    refold=p['k_refold']*ratio*q
    clearance=p['k_clear']*p['alpha_clear']*atp*(1-p['k_impair']*mg)*q
    cyt=(tnf+il)/2
    h=hazards(y,p)
    target=1-(1-hill(1-neur,p['K_motor'],p['n_motor']))*(1-p['k_inflam']*cyt)
    dy=np.array([
        p['k_syn_CL']*(1-c/p['K_CL'])-oxid+recycle-external,
        oxid-recycle-p['k_degrad_ox']*o,
        external-p['k_clear_ext_CL']*e,
        p['k_act_A']*r*(1-a)-p['k_deact_A']*a,
        p['k_act_T']*(1-z)-p['k_deact_T']*r*z,
        p['k_syn_aSyn']*(1-m)-aggregation+refold-p['k_turn']*m,
        aggregation-refold-clearance-p['k_secrete']*q,
        p['k_secrete']*q-p['k_clear_ext_aSyn']*x,
        p['k_CI_repair']*(p['CI_max']-ci)-p['k_CI_ROS']*r*ci-p['k_CI_aSyn']*q*q*ci,
        p['k_SC_form']*ratio*(1-sc)-p['k_SC_dissoc']*(1-ratio)*sc,
        p['k_resp']*ci*sc*(1-v)-p['k_ATPase']*v-p['k_leak']*v-p['k_mPTP_depol']*pore*v,
        p['k_ATP_syn']*v*ci*sc*(1-atp)-p['k_ATP_consume']*atp-p['k_ATP_mPTP']*pore*atp,
        p['k_mPTP_open']*hill(r,p['K_mPTP'],p['n_Hill'])*(1-pore)-p['k_mPTP_close']*pore*v,
        p['k_ROS_basal']+p['k_shunt']*(1-ci)*sc-(p['k_SOD2']*s+p['k_GPx4']*g+p['k_ROS_other'])*r,
        p['k_SOD2_on']*(1-s)-p['k_SOD2_off']*r*s,
        p['k_GPx4_on']*(1-g)-p['k_GPx4_off']*r*g,
        (p['k_PINK1_on']*(1-v)+p['k_PINK1_CL']*e)*(1-pink)-p['k_PINK1_off']*v*pink,
        (p['k_damage_mPTP']*pore+p['k_damage_energy']*(1-atp)+p['k_damage_dpsi']*(1-v))*(1-damage)-(p['k_repair']*pink+p['k_biogen']*atp)*damage,
        0. if freeze_neurons else -np.sum(h)*neur,
        (p['k_act_MG']*x+p['k_auto']*cyt)*(1-mg)-p['k_deact_MG']*mg,
        p['k_TNF_rel']*mg*(1-tnf)-p['k_TNF_deg']*tnf,
        p['k_IL1b_rel']*mg*(1-il)-p['k_IL1b_deg']*il,
        p['k_motor_rate']*(target-motor)])
    if len(y)==26:
        dy=np.r_[dy, np.zeros(3) if freeze_neurons else h*neur]
    return dy

def initial_state(p=None):
    """Healthy biochemical equilibrium; DA fixed at 1 during burn-in only."""
    p=parameters(healthy=True) if p is None else {**p,'CI_max':1.,'alpha_clear':1.}
    y=np.array([.9,.05,.001,.1,.9,.95,.05,.0025,.95,.95,.95,.94,.02,.1,.83,.67,.2,.06,1.,.13,.2,.16,.05])
    from fast_kernel import make_rhs
    sol=solve_ivp(make_rhs(p,freeze=True),(0,5000),y,method='LSODA',rtol=1e-9,atol=1e-11)
    if not sol.success: raise RuntimeError(sol.message)
    residual=np.max(np.abs(rhs(0,sol.y[:,-1],p,freeze_neurons=True)))
    if residual>1e-7: raise RuntimeError(f'Healthy biochemical burn-in not converged: {residual}')
    return sol.y[:,-1]

def simulate(p=None,days=50,dose=0.,y0=None,species='mouse',method='LSODA',
             rtol=1e-8,atol=1e-10,max_step=24.,points=None,reference_mode='periodic_mean',
             mechanism='repair',cumulative=True):
    p=parameters(species) if p is None else p.copy()
    y0=initial_state(p) if y0 is None else np.asarray(y0)[:23].copy()
    if cumulative: y0=np.r_[y0,np.zeros(3)]
    reference=p['exposure_reference'] if reference_mode=='periodic_mean' else p['C_ref']
    if reference_mode not in ['periodic_mean','manuscript']: raise ValueError(reference_mode)
    pk=PK(p,dose) if dose else None
    times=np.linspace(0,24*days,points or (int(days*2)+1))
    from fast_kernel import make_rhs
    sol=solve_ivp(make_rhs(p,pk,reference,False,mechanism),(0,24*days),y0,
                  method=method,t_eval=times,rtol=rtol,atol=atol,max_step=max_step)
    if not sol.success or not np.all(np.isfinite(sol.y)): raise RuntimeError(sol.message)
    if np.min(sol.y)<-1e-7: raise FloatingPointError('Negative state')
    if np.max(sol.y[FRACTIONS])>1+1e-7: raise FloatingPointError('Fraction exceeds one')
    return sol

def outputs(y):
    y=np.asarray(y)
    return np.array([y[0]/(y[0]+y[1]),y[8],y[13],y[6],y[11],y[18],y[22]])

def normalized_endpoints(p,days=35):
    """CI relative to healthy, ROS relative to healthy, native CL depletion."""
    init=initial_state(p)
    disease=simulate(p,days=days,y0=init,points=2).y[:,-1]
    return np.array([disease[8]/init[8],disease[13]/init[13],1-disease[0]/init[0]])
