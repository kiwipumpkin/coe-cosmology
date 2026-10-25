"""lc2: continuous circle construction (Chris's confirmed model), exact GR (LTB + central mass, null geodesics).
 - Shells: common start t=0 at R = R_i*r0 (R_i small, default 1e-3; R_i->0 = all shells from one point). tB(x) = -T_kepler(x, R_i r0).
 - Observer O = comoving element r0=1. Element temperature (cools with its own volume expansion): T = K*(J_obs(tn)/J)^(1/3), J=R^2 R'/r0^2.
   K = observer's element temperature today. ASSUMPTION: K = 2.725 K (not tuned).  t* = time O's element is at 3000 K.
 - Continuum of O_n sees C_(n-1): the light reaching O at t was emitted from the emission circle of radius c(t-t*) about
   O(t*) ... in GR: the intersection of O's past null cone with the slice t=t*.  For today's sky: emitter = ray point at t_e = t*.
 - T_obs = T_emit/(1+z), T_emit = emitter element temperature at t* (no matter-edge fallback, no K tuning).
 - Flags: 'beyond_CoE' = ray passed through the centre (r0 -> 0) before t*; 'stopped' = integrator stopped (never scored);
   'no_matter' = emitter density < 1e-12 of ours (reported, not hidden)."""
import numpy as np; from lc import *
class LC2(LC):
    def __init__(s,P,Ri=1e-3): super().__init__(P); s.Ri=Ri
    def tB(s,x):
        x=np.atleast_1d(np.asarray(x,float)); return -s._T(x,s.Ri*np.exp(x))
def make2(law,z,Ri=1e-3):
    L,tn=make(law,z); L2=LC2(L.P,Ri); return L2,tn
def tstar(L,tn,K=2.725):
    X0=np.array([0.]); lJo=lnJ(L,tn,X0)[0]
    f=lambda t: (lJo-lnJ(L,t,X0)[0])/3+np.log(K/3000.)
    tl=1e-9*tn
    return brentq(f,tl,tn) if f(tl)>0>f(tn) else None
def emitters(L,tn,ts,ells,K=2.725,budget=30000):
    o=obs(L,tn); X0=np.array([0.]); lJo=lnJ(L,tn,X0)[0]; out=[]
    for e in ells:
        r=ray(L,tn,np.radians(e),tmin=ts,rtol=1e-9,budget=budget); d=dict(ell=float(e),stop=r['stop'],pts=r)
        if r['stop']=='reached_tmin' and len(r['t']) and abs(r['t'][-1]-ts)<1e-6*tn:
            xe=np.array([r['x'][-1]]); Re=L.R(ts,xe)[0]
            lr=np.log(L.rho_rel(ts,xe)[0]/o['rho']); Te=K*np.exp((lJo-lnJ(L,ts,xe)[0])/3)
            d.update(status='beyond_CoE' if r['nc'].max()>0 else 'emitter',x_e=float(xe[0]),R_e=float(Re),phi=float(r['phi'][-1]),z=float(np.expm1(r['lnz1'][-1])),
                     T_emit=float(Te),T_obs=float(Te*np.exp(-r['lnz1'][-1])),lnrho_rel=float(lr),no_matter=bool(not np.isfinite(lr) or lr<np.log(1e-12)))
        else: d['status']='stopped:'+r['stop']
        out.append(d)
    return out
# ---- dt_step construction (Chris approved ../sketch/dt_step.png) ----
# Emission circle at time t_e: centred on OUR position O(t_e) (areal R_o(t_e) on the CoE->us axis), radius c(t_e - t*).
# Emitter for a direction: first point back along today's past null geodesic where |x_ray(t) - x_O(t)| = c (t - t*)
# (Euclidean distance in areal coordinates (R cos phi, R sin phi) -- ASSUMPTION for the GR form of 'radius'); t_e >= t* by construction.
def emitters_circle(L,tn,ts,ells,K=2.725,budget=30000):
    o=obs(L,tn); X0=np.array([0.]); lJo=lnJ(L,tn,X0)[0]; out=[]
    def g(t,y):
        Rr=L.R(t,np.array([y[0]]))[0]; Ro=L.R(t,X0)[0]
        d=np.sqrt(Rr**2+Ro**2-2*Rr*Ro*np.cos(y[1])); return d-(t-ts)
    g.terminal=True
    for e in ells:
        r=ray(L,tn,np.radians(e),tmin=ts,rtol=1e-9,budget=budget,evx=g); d=dict(ell=float(e),stop=r['stop'],pts=r)
        if r['stop']=='circle_hit':
            te=r['t'][-1]; xe=np.array([r['x'][-1]]); Re=L.R(te,xe)[0]
            lr=np.log(L.rho_rel(te,xe)[0]/o['rho']); Te=K*np.exp((lJo-lnJ(L,te,xe)[0])/3)
            d.update(status='beyond_CoE' if r['nc'].max()>0 else 'emitter',t_e=float(te),x_e=float(xe[0]),R_e=float(Re),phi=float(r['phi'][-1]),R_o_te=float(L.R(te,X0)[0]),
                     z=float(np.expm1(r['lnz1'][-1])),T_emit=float(Te),T_obs=float(Te*np.exp(-r['lnz1'][-1])),lnrho_rel=float(lr),
                     plasma=bool(Te>=3000),no_matter=bool(not np.isfinite(lr) or lr<np.log(1e-12)))
        else: d['status']='stopped:'+r['stop']
        out.append(d)
    return out
