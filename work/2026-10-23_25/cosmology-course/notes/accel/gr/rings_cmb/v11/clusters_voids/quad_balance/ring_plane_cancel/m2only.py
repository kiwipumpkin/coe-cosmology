# Handful of in-plane point bodies (paper Eq.2 tidal kernel via seq_core.pattern) fitted to supply ONLY the m=2 part of the smooth ring R.
import sys,os,json,numpy as np
from scipy.optimize import least_squares
sys.path.insert(0,'/workspace/clump_scan/sequential');H=os.getcwd();os.chdir('/workspace/clump_scan/sequential');import seq_core as S;os.chdir(H)
th=np.radians(np.arange(360)) if not hasattr(S,'TH') else S.TH;R=S.R;F=S.FLOOR
def proj(v,ms):
    o=np.zeros_like(v)
    for m in ms:o+=2*np.mean(v*np.cos(m*th))*np.cos(m*th)+2*np.mean(v*np.sin(m*th))*np.sin(m*th)
    return o
def mm(v):m=np.arange(1,61)[:,None];a=np.hypot(2*np.mean(v*np.cos(m*th),-1),2*np.mean(v*np.sin(m*th),-1));return a/np.sqrt(2)/F
R2=proj(R,[2]);KG=S.KG if hasattr(S,'KG') else 1
def lay(p,N):
    q=p.reshape(N,3);rad=14.3+15.7/(1+np.exp(-q[:,0]));X=np.c_[rad*np.cos(q[:,1]),rad*np.sin(q[:,1]),np.zeros(N)];M=q[:,2]*1e52;return X,M
def ok(X):return np.all(np.linalg.norm(X-S.C_G,axis=1)>=S.A_IN_G) and np.all(np.linalg.norm(X,axis=1)<=30)
rng=np.random.default_rng(1);out={}
for N in (1,2,4,8):
    best=None
    for t in range(15):
        p0=np.c_[np.log(rng.uniform(0.3,15,N)),rng.uniform(0,2*np.pi,N),rng.normal(0,1,N)].ravel()
        f=lambda p:(lambda XM:(S.pattern(*XM)-R2)/F/np.sqrt(360)*(1 if ok(XM[0]) else 30))(lay(p,N))
        try:r=least_squares(f,p0,max_nfev=300)
        except Exception:continue
        if best is None or r.cost<best.cost:best=r
    X,M=lay(best.x,N);w=S.pattern(X,M);res=R-w;a=mm(res);aw=mm(w)
    out[N]=dict(m2_after=float(a[1]),others_max_1_12=float(np.max(np.delete(a[:12],1))),per_m_after=a[:12].tolist(),body_leak_per_m=np.delete(aw[:12],1).tolist(),total=float(S.rms(res)/F),X=X.tolist(),M=M.tolist(),valid=bool(ok(X)))
    print(N,'m2',round(a[1],3),'per-m',np.round(a[:12],2).tolist(),'leak(m!=2) max',round(np.max(np.delete(aw[:12],1)),3),'valid',ok(X),flush=True)
json.dump(out,open('m2only_in30.json','w'),indent=1)
