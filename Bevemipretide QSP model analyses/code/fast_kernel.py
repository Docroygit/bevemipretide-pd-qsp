"""Numba-compiled PD kernel, algebraically identical to model.rhs.

The Python implementation remains the independent reference for tests.
No fast-math approximations or reduced precision are enabled.
"""
import numpy as np
from numba import njit,types
from numba.typed import Dict

def typed_parameters(p):
    out=Dict.empty(key_type=types.unicode_type,value_type=types.float64)
    for k,v in p.items():out[k]=float(v)
    return out

@njit(cache=True)
def H(x,k,n):
    if x < -1e-5:raise ValueError('Negative Hill input')
    if x<=0:return 0.
    z=(x/k)**n
    return z/(1+z)

@njit(cache=True)
def derivative(y,p,drug,freeze,protection):
    c,o,e,a,z,m,q,x,ci,sc,v,atp,pore,r,s,g,pink,damage,neur,mg,tnf,il,motor=y[:23]
    ratio=c/(c+o) if c+o>0 else 0.
    oxid=p['k_ox']*r*c*(1+p['k_ALCAT']*a)+p['k_aSyn_ox']*q*c
    recycle=p['k_red']*z*o
    if protection:oxid/=(1+drug)
    else:recycle+=p['k_drug']*drug*o
    external=p['k_ext']*c*r
    stress=1+p['k_ROS_agg']*r+p['k_CL_loss']*(1-ratio)+p['k_seed_CL']*e
    aggregation=p['k_agg']*stress*m;refold=p['k_refold']*ratio*q
    clearance=p['k_clear']*p['alpha_clear']*atp*(1-p['k_impair']*mg)*q
    cyt=(tnf+il)/2
    hm=p['k_death']*H(damage,p['K_death'],p['n_death'])
    hi=p['k_death_inflam']*H(cyt,p['K_inflam_death'],p['n_inflam_death'])
    hq=p['k_death_aSyn']*H(q,p['K_aSyn_death'],p['n_aSyn_death'])
    target=1-(1-H(1-neur,p['K_motor'],p['n_motor']))*(1-p['k_inflam']*cyt)
    dy=np.empty(len(y))
    dy[0]=p['k_syn_CL']*(1-c/p['K_CL'])-oxid+recycle-external
    dy[1]=oxid-recycle-p['k_degrad_ox']*o
    dy[2]=external-p['k_clear_ext_CL']*e
    dy[3]=p['k_act_A']*r*(1-a)-p['k_deact_A']*a
    dy[4]=p['k_act_T']*(1-z)-p['k_deact_T']*r*z
    dy[5]=p['k_syn_aSyn']*(1-m)-aggregation+refold-p['k_turn']*m
    dy[6]=aggregation-refold-clearance-p['k_secrete']*q
    dy[7]=p['k_secrete']*q-p['k_clear_ext_aSyn']*x
    dy[8]=p['k_CI_repair']*(p['CI_max']-ci)-p['k_CI_ROS']*r*ci-p['k_CI_aSyn']*q*q*ci
    dy[9]=p['k_SC_form']*ratio*(1-sc)-p['k_SC_dissoc']*(1-ratio)*sc
    dy[10]=p['k_resp']*ci*sc*(1-v)-p['k_ATPase']*v-p['k_leak']*v-p['k_mPTP_depol']*pore*v
    dy[11]=p['k_ATP_syn']*v*ci*sc*(1-atp)-p['k_ATP_consume']*atp-p['k_ATP_mPTP']*pore*atp
    dy[12]=p['k_mPTP_open']*H(r,p['K_mPTP'],p['n_Hill'])*(1-pore)-p['k_mPTP_close']*pore*v
    dy[13]=p['k_ROS_basal']+p['k_shunt']*(1-ci)*sc-(p['k_SOD2']*s+p['k_GPx4']*g+p['k_ROS_other'])*r
    dy[14]=p['k_SOD2_on']*(1-s)-p['k_SOD2_off']*r*s
    dy[15]=p['k_GPx4_on']*(1-g)-p['k_GPx4_off']*r*g
    dy[16]=(p['k_PINK1_on']*(1-v)+p['k_PINK1_CL']*e)*(1-pink)-p['k_PINK1_off']*v*pink
    dy[17]=(p['k_damage_mPTP']*pore+p['k_damage_energy']*(1-atp)+p['k_damage_dpsi']*(1-v))*(1-damage)-(p['k_repair']*pink+p['k_biogen']*atp)*damage
    dy[18]=0. if freeze else -(hm+hi+hq)*neur
    dy[19]=(p['k_act_MG']*x+p['k_auto']*cyt)*(1-mg)-p['k_deact_MG']*mg
    dy[20]=p['k_TNF_rel']*mg*(1-tnf)-p['k_TNF_deg']*tnf
    dy[21]=p['k_IL1b_rel']*mg*(1-il)-p['k_IL1b_deg']*il
    dy[22]=p['k_motor_rate']*(target-motor)
    if len(y)==26:
        dy[23]=0. if freeze else hm*neur
        dy[24]=0. if freeze else hi*neur
        dy[25]=0. if freeze else hq*neur
    return dy

@njit(cache=True)
def forced_derivative(t,y,p,lam,coef,weights,interval,reference,freeze,protection):
    n=int(np.floor(t/interval+1e-12))+1
    tau=max(0.,t-(n-1)*interval)
    amount=0.
    for j in range(5):
        geo=np.expm1(lam[j]*n*interval)/np.expm1(lam[j]*interval)
        amount+=weights[j]*coef[j]*np.exp(lam[j]*tau)*geo
    return derivative(y,p,amount/reference,freeze,protection)

def make_rhs(p,pk=None,reference=1.,freeze=False,mechanism='repair'):
    typed=typed_parameters(p);protection=(mechanism=='protection')
    if pk is None:return lambda t,y:derivative(y,typed,0.,freeze,protection)
    return lambda t,y:forced_derivative(t,y,typed,pk.lam,pk.coef,pk.vec[4],pk.interval,reference,freeze,protection)
