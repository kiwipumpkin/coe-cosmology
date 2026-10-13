"""Same fits as ../solution_space, but every object strictly INSIDE our visible sphere (between us and the CMB ring).
Imports ../solution_space/ss_core.py read-only (same R, FLOOR, KG, tidal kernel unit_field == lattice_core.field3).
Geometry (Gpc): object centre distance d from us with  d >= Rcore + GAP_US  (we are not inside any core / void)
                                               and  d + Rcore <= 14.3 - GAP_RING  (core never touches the ring sphere).
Net mass: compensated LOCALLY, i.e. -net spread smoothly over the visible sphere with the exp-lock density shape (FLOC);
the surround compensation FC of the outside fits is used only as a cross-check (re-scored, not refit)."""
import os, sys, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); SS=os.path.join(os.path.dirname(HERE),"solution_space")
sys.path.insert(0,SS); _c=os.getcwd()
import ss_core as S
from ss_core import R,FLOOR,KG,RV,rms,C_G,A_IN_G,GPC,RHO,rho_gpc,TH,unit_field,amps,vol,DELTA_C,GAP
os.chdir(_c)
GAP_US=0.1; GAP_RING=0.1
RLO,RHI=0.5,4.0; MLO=1e49
M_VIS=vol(RV)*GPC**3*RHO          # ~1.04e54 kg visible-sphere matter budget

def _floc(n=4_000_000,seed=1):
    """pattern per 1 kg spread over the visible ball (r<RV) with weight rho_exp-lock. Uniform ball = exactly 0 by the shell
    theorem (field at every ring point is radial & equal, observer term 0) -> only the gradient part is sampled (low variance)."""
    from scipy.stats import qmc
    u=qmc.Sobol(3,seed=seed).random(n); r=RV*u[:,0]**(1/3); ct=2*u[:,1]-1; ph=2*np.pi*u[:,2]; st=np.sqrt(1-ct*ct)
    X=np.c_[r*st*np.cos(ph),r*st*np.sin(ph),r*ct]; w=rho_gpc(X); w=w/w.mean()-1.0
    keep=np.abs(w)>0; F=np.zeros(360)
    for i in range(0,n,50000):
        sl=slice(i,i+50000); F+=(w[sl])@unit_field(X[sl])
    return F/n
# Analytic: a ball (radius a) with linear density gradient rho_bar*(x.n)/L gives, ring minus observer, a PURE dipole
#   -(6pi/5) G g a^2 (n.xhat) = -(9/10) G M /(a L) (n.xhat).  n = toward the CoE. Matches the Sobol MC _floc() to corr 0.9998
#   (MC rms ratio 0.994; m>=2 < 2%). Uniform part contributes exactly zero.
from web_core import L_G
_Ch=C_G/np.linalg.norm(C_G)
FLOC=-0.9/KG*(RV/L_G)*(np.cos(TH)*_Ch[0]+np.sin(TH)*_Ch[1]); FLOC=FLOC-FLOC.mean()

class InModel:
    def __init__(self,kinds,mhi=2e52,comp="local"):
        self.k=np.array(list(kinds)); self.n=len(kinds); self.mhi=mhi; self.comp=comp
        self.lo=np.tile([0,-4*np.pi,-np.pi/2,0],self.n); self.hi=np.tile([1,4*np.pi,np.pi/2,1],self.n)
    def decode(self,p):
        p=p.reshape(self.n,4); u,lon,lat,s=p.T
        M=np.zeros(self.n); Rv=np.zeros(self.n); Rc=np.zeros(self.n); isC=self.k=="C"
        M[isC]=10**(np.log10(MLO)+s[isC]*np.log10(self.mhi/MLO))
        Rv[~isC]=RLO+s[~isC]*(RHI-RLO); Rc[~isC]=Rv[~isC]
        Rc[isC]=(3*M[isC]/(4*np.pi*DELTA_C*RHO*GPC**3))**(1/3)
        dmin=Rc+GAP_US; dmax=RV-GAP_RING-Rc
        d=dmin+u*(dmax-dmin)
        dirs=np.c_[np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]; X=dirs*d[:,None]
        isV=self.k=="V"; M[isV]=-0.70*rho_gpc(X[isV])*vol(Rv[isV])*GPC**3
        return X,M,Rv,Rc
    def pattern(self,p,comp=None):
        comp=comp or self.comp; X,M,Rv,Rc=self.decode(p)
        w=M@unit_field(X)
        if comp=="local": w=w-M.sum()*FLOC
        elif comp=="surround": w=w-(M.sum()/KG)*S.FC
        return w
    def violations(self,p):
        X,M,Rv,Rc=self.decode(p); v=[np.maximum(0,A_IN_G+Rc-np.linalg.norm(X-C_G,axis=1))]
        if self.n>1:
            i,j=np.triu_indices(self.n,1); d=np.linalg.norm(X[i]-X[j],axis=1); v.append(np.maximum(0,Rc[i]+Rc[j]+GAP-d))
        v.append([max(0.0,M[M>0].sum()/M_VIS-1.0)])        # concentrated mass cannot exceed the sphere's whole budget
        return np.concatenate(v)
    def leftover(self,p,comp=None): return rms(R-self.pattern(p,comp))/FLOOR

def check():
    """kernel sanity for interior points: unit_field (field3) vs direct Newtonian a(ring)-a(us) projected on the sightline."""
    rng=np.random.default_rng(0); out=[]
    for _ in range(200):
        x=rng.normal(size=3); x/=np.linalg.norm(x); X=x*rng.uniform(0.5,13.5)
        ring=np.c_[np.cos(TH),np.sin(TH),np.zeros(360)]*RV
        G=6.674e-11; dr=X-ring; a_r=G*(dr*GPC)/np.linalg.norm(dr*GPC,axis=1)[:,None]**3
        a_o=G*X*GPC/np.linalg.norm(X*GPC)**3
        proj=((a_r-a_o)*ring/RV).sum(1); proj-=proj.mean()
        f=unit_field(X[None])[0]; c=np.dot(f,proj)/np.linalg.norm(f)/np.linalg.norm(proj); out.append(c)
    return min(out)
