"""lc3: each element emits once, at T=3000 K, at its own recombination crossing (T_elem falls through 3000 K).
Emitter per direction = event on our exact-GR past light cone where the element there has T_elem = 3000 K, with t_e >= t*.
T_elem = K (J_obs,today/J)^(1/3); ASSUMPTION K = 2.725 K = our own element's temperature today (fixes t*, no sky rescale).
T_obs = 3000/(1+z). Ray reaching t* without meeting a crossing -> 'no_crossing_after_t*' (no light, not scored)."""
import numpy as np; from lc2 import *
def emitters_rec(L,tn,ts,ells,K=2.725,budget=40000):
    X0=np.array([0.]); lJo=lnJ(L,tn,X0)[0]; thr=3*np.log(3000./K); out=[]
    def g(t,y): return (lJo-lnJ(L,t,np.array([y[0]]))[0])-thr
    g.terminal=True
    for e in ells:
        r=ray(L,tn,np.radians(e),tmin=ts,rtol=1e-9,budget=budget,evx=g); d=dict(ell=float(e),stop=r['stop'],pts=r)
        if r['stop']=='circle_hit':
            te=r['t'][-1]; xe=np.array([r['x'][-1]]); o=obs(L,tn); lr=np.log(L.rho_rel(te,xe)[0]/o['rho'])
            d.update(status='beyond_CoE' if r['nc'].max()>0 else 'emitter',t_e=float(te),x_e=float(xe[0]),R_e=float(L.R(te,xe)[0]),phi=float(r['phi'][-1]),
                     z=float(np.expm1(r['lnz1'][-1])),T_obs=float(3000*np.exp(-r['lnz1'][-1])),lnrho_rel=float(lr),no_matter=bool(not np.isfinite(lr) or lr<np.log(1e-12)))
        elif r['stop']=='reached_tmin': d['status']='no_crossing_after_t*'
        else: d['status']='stopped:'+r['stop']
        out.append(d)
    return out
