"""Single shared gravitational frequency factor (all CMB codes).
nu_obs/nu_emit = exp[(Phi_emit - Phi_obs)/c^2]  (Phi<0; photon climbing out of a deeper well at the emitter is redshifted).
Weak-field/Schwarzschild equivalent used by cmb-cancel.py: sqrt[(1-2mu/r_e)/(1-2mu/r_o)] (same sign)."""
import numpy as np
def grav_factor(phi_emit,phi_obs,c2): return np.exp((np.asarray(phi_emit)-np.asarray(phi_obs))/c2)
import sys; sys.path.insert(0,'/workspace/coe_grz'); import grav_redshift as GRZ
def gr_S(rs,xe,ye,Ue,xo,yo,Uo,kx,ky):
    """nu_obs/nu_emit = 1/(1+z), exact point-mass Schwarzschild (Cosmology's grav_redshift.py), radial movers (ell=0)
    with proper radial speed U=dr/dtau in units of c -> eps = sqrt(U^2 + 1 - rs/r). Photon: straight-line direction k
    (no lensing; ASSUMPTION), impact L=|x x k|, radial sign = sign(x.k). Returns nan in forbidden regions."""
    out=np.full(np.shape(xe),np.nan)
    for i in np.ndindex(np.shape(xe)):
        try:
            re=np.hypot(xe[i],ye[i]); ro=np.hypot(xo,yo)
            Ae=1-rs/re; Ao=1-rs/ro
            ee=np.sqrt(Ue[i]**2+Ae); eo=np.sqrt(Uo**2+Ao)
            if not(np.isfinite(ee) and np.isfinite(eo)): continue
            L=abs(xe[i]*ky[i]-ye[i]*kx[i])
            ske=1 if xe[i]*kx[i]+ye[i]*ky[i]>=0 else -1; sko=1 if xo*kx[i]+yo*ky[i]>=0 else -1
            z1=GRZ.one_plus_z(rs,(re,ee,0.0,np.sign(Ue[i]) or 1),(ro,eo,0.0,np.sign(Uo) or 1),L,ske,sko)
            out[i]=1/z1
        except Exception: pass
    return out
