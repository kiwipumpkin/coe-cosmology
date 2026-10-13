"""Diagnostic: every body confined to the CoE core (within 43.9 Gpc of the CoE, wholly outside our 14.3 Gpc sphere).
Same masses/voids, same FC compensation as ss_core. How far can core bodies alone bring the leftover down?
usage: python3 run_coreonly.py NPROC NSTARTS -> raw_coreonly/<kinds>.json"""
import os,sys,json,time,numpy as np
os.environ.setdefault("OMP_NUM_THREADS","1")
from multiprocessing import Pool
from co_core import *
import fit_run
class CoreOnly(CO):
    def decode(self,p):
        p=p.reshape(self.n,4); u,lon,lat,s=p.T; isC=self.k=="C"
        M=np.zeros(self.n); Rv=np.zeros(self.n); Rc=np.zeros(self.n)
        M[isC]=10**(np.log10(MLO)+s[isC]*np.log10(MHI/MLO)); Rv[~isC]=RLO+s[~isC]*(RHI-RLO); Rc[~isC]=Rv[~isC]
        Rc[isC]=(3*M[isC]/(4*np.pi*DELTA_C*RHO*GPC**3))**(1/3)
        U=np.c_[np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]      # direction seen from the CoE
        X=C_G+U*((A_IN_G-Rc)*u**(1/3))[:,None]
        isV=self.k=="V"; M[isV]=-0.70*rho_gpc(X[isV])*vol(Rv[isV])*GPC**3
        return X,M,Rv,Rc
    def violations(self,p):
        X,M,Rv,Rc=self.decode(p); v=[np.maximum(0,RMIN+Rc-np.linalg.norm(X,axis=1))]
        if self.n>1:
            i,j=np.triu_indices(self.n,1); v.append(np.maximum(0,Rc[i]+Rc[j]+GAP-np.linalg.norm(X[i]-X[j],axis=1)))
        return np.concatenate(v)
CFG=["C","CC","CCCC","CCVV","CCCCCC","CCCCVV","CCCCVVVV","CCCCCCVVVV","CCCCCCVVVVVV"]
def job(a):
    k,p0,tag=a; m=CoreOnly(k); t0=time.time(); p=fit_run.fit_one(m,p0); p=fit_run.fit_one(m,p)
    X,M,Rv,Rc=m.decode(p); w=m.pattern(p)
    return dict(scen="coreonly",name=k,start=tag,left=float(rms(R-w)/FLOOR),ok=bool(m.violations(p).max(initial=0)<=1e-3),p=p.tolist(),X=X.tolist(),M=M.tolist(),Rv=Rv.tolist(),
                resid=(R-w).tolist(),amps=(amps(R-w)/FLOOR).tolist(),core_depth=core_depth(X).tolist(),dist=np.linalg.norm(X,axis=1).tolist(),secs=time.time()-t0)
if __name__=="__main__":
    npr,ns=int(sys.argv[1]),int(sys.argv[2]); os.makedirs("raw_coreonly",exist_ok=True); tasks=[]
    for ci,k in enumerate(CFG):
        rng=np.random.default_rng(77+ci); n=len(k)
        for i in range(ns):
            p=np.zeros((n,4)); p[:,0]=rng.uniform(0,1,n); p[:,1]=rng.uniform(-np.pi,np.pi,n); p[:,2]=np.arcsin(rng.uniform(-1,1,n)); p[:,3]=rng.uniform(0,1,n)
            tasks.append((k,p.ravel(),f"cold{i}"))
    tasks.sort(key=lambda t:-len(t[0])); res={}; t0=time.time()
    with Pool(npr) as pool:
        for r in pool.imap_unordered(job,tasks):
            res.setdefault(r["name"],[]).append(r)
            print(f"{time.time()-t0:6.0f}s {r['name']:13s} {r['start']:7s} left={r['left']:.2f} ok={r['ok']} dist={np.round(r['dist'],1).tolist()}",flush=True)
            json.dump(res[r["name"]],open(f"raw_coreonly/{r['name']}.json","w"))
    print("done",time.time()-t0)
