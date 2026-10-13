"""Order / entropy score of a web layout, measured against the starting web (Chris, 10 Oct).
Components (all bodies alive in the modeled region, ball 30 Gpc):
  radial  : counts in 30 equal-volume shells around us (0-30 Gpc; the shell straddling the CMB ring is one of them)
  ring    : counts vs distance from the CMB ring, 24 bins of equal tube volume out to 14.3 Gpc + the rest of the region
  angular : counts in 12 x 24 equal-area sky cells (cos(theta) x phi)
  -> volume-weighted Shannon entropy  H = -sum p log(p / w)  (w = bin's share of region volume, Monte Carlo), so uniform-in-volume = max = 0
  link_cc, link_cv, link_vv : nearest-neighbour link lengths, cluster->cluster, cluster->void, void->void (the link types of
              void_mass_network/census.py, as nearest links), Shannon entropy of the histogram on fixed bins from the start
  ncount  : neighbour count within 3 Gpc (all bodies), Shannon entropy of the integer histogram
A lower entropy than at the start = more ordered. Rule: every component must stay >= its starting value."""
import numpy as np
from scipy.spatial import cKDTree
RMIN=14.3; RT=30.0
_rng=np.random.default_rng(1); _P=_rng.uniform(-RT,RT,(4_000_000,3)); _P=_P[np.linalg.norm(_P,axis=1)<=RT]
def dring(X): X=np.atleast_2d(X); return np.sqrt((np.hypot(X[:,0],X[:,1])-RMIN)**2+X[:,2]**2)
NR=30; NG=24; NC,NPH=12,24
def u_radial(X): return np.clip(NR*(np.linalg.norm(X,axis=1)/RT)**3,0,NR-1e-9)
def u_ring(X):
    d=dring(X); return np.where(d<RMIN,NG*(d/RMIN)**2,NG+0.5)       # equal tube volume (2 pi R * pi d^2) bins; last bin = rest
def u_ang(X):
    s=np.linalg.norm(X,axis=1); ct=np.clip(X[:,2]/np.maximum(s,1e-9),-1,1); ph=np.arctan2(X[:,1],X[:,0])%(2*np.pi)
    return np.clip((ct+1)/2*NC,0,NC-1e-9),np.clip(ph/(2*np.pi)*NPH,0,NPH-1e-9)
def _w():
    w_r=np.bincount(u_radial(_P).astype(int),minlength=NR)/len(_P); w_g=np.bincount(u_ring(_P).astype(int),minlength=NG+1)/len(_P)
    a,b=u_ang(_P); w_a=np.bincount(a.astype(int)*NPH+b.astype(int),minlength=NC*NPH)/len(_P); return w_r,w_g,w_a
W_R,W_G,W_A=_w()
def Hvol(c,w):
    p=c/np.maximum(c.sum(),1e-300); m=p>0; return float(-np.sum(p[m]*np.log(p[m]/w[m])))
def Hsh(c):
    p=c/np.maximum(c.sum(),1e-300); m=p>0; return float(-np.sum(p[m]*np.log(p[m])))
def counts_hard(X):
    a,b=u_ang(X)
    return (np.bincount(u_radial(X).astype(int),minlength=NR).astype(float),np.bincount(u_ring(X).astype(int),minlength=NG+1).astype(float),
            np.bincount(a.astype(int)*NPH+b.astype(int),minlength=NC*NPH).astype(float))
def soft1(u,n):
    """linear interpolation between bin centres (differentiable in u)"""
    c=np.zeros(n); t=u-0.5; i0=np.floor(t).astype(int); f=t-i0; i1=i0+1
    np.add.at(c,np.clip(i0,0,n-1),1-f); np.add.at(c,np.clip(i1,0,n-1),f); return c
def counts_soft(X):
    a,b=u_ang(X); ca=np.zeros(NC*NPH); ta=a-0.5; tb=b-0.5; ia=np.floor(ta).astype(int); ib=np.floor(tb).astype(int); fa=ta-ia; fb=tb-ib
    for da,wa in ((0,1-fa),(1,fa)):
        for db,wb in ((0,1-fb),(1,fb)):
            np.add.at(ca,np.clip(ia+da,0,NC-1)*NPH+(ib+db)%NPH,wa*wb)
    return soft1(u_radial(X),NR),soft1(u_ring(X),NG+1),ca
class Score:
    """fixed bins from the starting layout; components(X, r, isV) -> dict of entropies"""
    def __init__(self,X0,r0,isV):
        L=self.links(X0,r0>0,isV); self.bins={k:np.linspace(0,np.percentile(v,99.5)*1.5,41) for k,v in L.items() if k!="ncount"}
        self.S0=self.components(X0,r0,isV)
    @staticmethod
    def links(X,alive,isV):
        Xa=X[alive]; va=isV[alive]; out={}
        Tc=cKDTree(Xa[~va]); Tv=cKDTree(Xa[va])
        out["link_cc"]=Tc.query(Xa[~va],k=2)[0][:,1]; out["link_cv"]=Tv.query(Xa[~va],k=1)[0]; out["link_vv"]=Tv.query(Xa[va],k=2)[0][:,1]
        Ta=cKDTree(Xa); out["ncount"]=np.array([len(x)-1 for x in Ta.query_ball_point(Xa,3.0)])
        return out
    def components(self,X,r,isV):
        al=(r>0)&(np.linalg.norm(X,axis=1)<=RT); cr,cg,ca=counts_hard(X[al]); S=dict(radial=Hvol(cr,W_R),ring=Hvol(cg,W_G),angular=Hvol(ca,W_A))
        L=self.links(X,r>0,isV)
        for k in ("link_cc","link_cv","link_vv"): S[k]=Hsh(np.histogram(np.clip(L[k],0,self.bins[k][-1]-1e-9),self.bins[k])[0].astype(float))
        S["ncount"]=Hsh(np.bincount(L["ncount"]).astype(float)); return S
    def ok(self,S,tol=0.0): return all(S[k]>=self.S0[k]-tol for k in self.S0)
    def deficit(self,S): return {k:max(0.0,self.S0[k]-S[k]) for k in self.S0}
