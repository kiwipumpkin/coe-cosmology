"""Open placement region for the sequential runs (patched into the ORIGINAL code at import time; originals untouched).
Region (for a reach rmax: 42.9 Gpc for B / outer, 30 Gpc for the A near fix):
   ball of radius rmax around us  UNION  the whole CoE core (|x - CoE| <= 43.9 Gpc, reaching 102.1 Gpc from us)
   * inside our 14.3 Gpc sphere allowed; a sphere may not straddle the CMB ring: wholly inside (|x|+r <= 14.2)
     or wholly outside (|x|-r >= 14.4)  [same no-straddle rule as the all-open fits]
   * nothing within 0.5 Gpc of us (|x|-r >= 0.5)
   * CoE core: no wall, no cap. A body counts as 'in the region' if its centre is within rmax of us,
     or the whole sphere lies inside the core.
Everything else (kernels, masses, packing, rules, seeds) is the original code."""
import os, numpy as np
MODE=os.environ.get("SEQ_OPEN","both")   # tidal = our sphere + shell out to 30 Gpc from us (effective tidal depth), core NOT filled; both | core (core filled, our sphere empty) | inside (our sphere filled, core blocked as before)
RMIN=14.3; C_G=None; A_IN=None
GAP_US=0.5; GAP_RING=0.1; RT=30.0
def setup(C,A):
    global C_G,A_IN; C_G=np.asarray(C,float); A_IN=float(A)
def region_ok(x,r,rmin=None,rmax=42.9):
    x=np.atleast_2d(x); r=np.atleast_1d(r); s=np.linalg.norm(x,axis=1); dc=np.linalg.norm(x-C_G,axis=1)
    noring=(s+r<=RMIN-GAP_RING)|(s-r>=RMIN+GAP_RING)
    if MODE=="core":   return (s-r>=RMIN)&((s<=rmax)|(dc+r<=A_IN))
    if MODE=="inside": return noring&(s-r>=GAP_US)&(s<=rmax)&(dc-r>=A_IN)
    if MODE=="tidalcore": return noring&(s-r>=GAP_US)&(s<=RT)          # ball(30) incl. the part of the core inside it; no wall
    if MODE=="tidal":  return noring&(s-r>=GAP_US)&(s<=RT)&(dc-r>=A_IN)
    reach=(s<=rmax)|(dc+r<=A_IN)
    return noring&(s-r>=GAP_US)&reach
def _dir(rng,n):
    u=rng.normal(size=(n,3)); return u/np.linalg.norm(u,axis=1,keepdims=True)
def rand_in_region(rng,r0,r1,n):
    """uniform over ball(r1) UNION core (r0 ignored: inside our sphere is open)."""
    if MODE in ("inside","tidal","tidalcore"):
        r1=RT if MODE in ("tidal","tidalcore") else r1
        return _dir(rng,n)*(r1*rng.uniform(size=n)**(1/3))[:,None]
    Vb=4/3*np.pi*r1**3; Vc=4/3*np.pi*A_IN**3; out=[]
    while len(out)<n:
        if rng.uniform()<Vb/(Vb+Vc):
            out.append(_dir(rng,1)[0]*r1*rng.uniform()**(1/3))
        else:
            x=C_G+_dir(rng,1)[0]*A_IN*rng.uniform()**(1/3)
            if np.linalg.norm(x)>r1 or rng.uniform()<0: out.append(x)     # core points already inside the ball come from the ball draw
            else:
                # keep uniform density over the union: reject core draws that land inside the ball (already covered)
                continue
    return np.array(out)
def region_volume(rmax,n=2_000_000,seed=0):
    rng=np.random.default_rng(seed); L=np.linalg.norm(C_G)+A_IN
    P=(rng.uniform(-1,1,(n,3)))*L
    s=np.linalg.norm(P,axis=1); dc=np.linalg.norm(P-C_G,axis=1)
    return float((2*L)**3*np.mean((s<=rmax)|(dc<=A_IN)))
