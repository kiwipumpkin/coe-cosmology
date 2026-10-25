# McVittie (H R/c=0.3, matter-only) + observer and shell moving through space at peculiar v0 toward 0 deg.
# Shell round in its own frame: start events built with the flat boost formula about the observer's start point
# (coordinate time and coordinates, a(t0)=1):  emitter  t = t0 + g v0 (1+e_x),  X = (2,0) + (g(e_x+v0^2), e_y)
#                                               observer t = t0 + g v0,        X = (2,0) + (g v0^2, 0)
import sys,pickle; sys.path.insert(0,'.'); import numpy as np, iso, multiprocessing as mp
from grcore import relsum
th=np.radians(np.arange(0,361,1.0)); H0=0.3; t0=2/(3*H0); C0=np.array([2.0,0])
def setup(m,u,v0,T):
    Mt=iso.Metric(m,t0); g=1/np.sqrt(1-v0*v0); ems=[]
    for t in th:
        e=np.array([np.cos(t),np.sin(t)]); ts=t0+g*v0*(1+e[0])
        at,h=iso.worldline(Mt,C0+np.array([g*(e[0]+v0*v0),e[1]]),relsum(v0,u,e),ts,T); ems.append((at,ts))
    oat,_=iso.worldline(Mt,C0+np.array([g*v0*v0,0]),np.array([v0,0]),t0+g*v0,T)
    return Mt,ems,oat
def ok180(m,v0,T):
    Mt=iso.Metric(m,t0); g=1/np.sqrt(1-v0*v0)
    oat,_=iso.worldline(Mt,C0+np.array([g*v0*v0,0]),np.array([v0,0]),t0+g*v0,T); Xo,vo=oat(T)
    so=iso.photon(Mt,Xo,0.0,T,t0-0.1)
    if so.t[-1]>t0: return False
    return so.sol(t0)[0] < C0[0]-1/g
def run(m,u,v0,T):
    Mt,ems,oat=setup(m,u,v0,T); Xo,vo=oat(T); return iso.observe(Mt,Xo,vo,ems,T)[0]
def job(v0):
    g=1/np.sqrt(1-v0*v0); lo=t0+2*g*v0+0.05; hi=lo+2
    while not ok180(0.1,v0,hi):
        hi=t0+(hi-t0)*1.5
        if hi-t0>2000: return v0,None,None,"light from behind never caught up by dt=2000"
    for _ in range(25):
        mid=(lo+hi)/2
        if ok180(0.1,v0,mid): hi=mid
        else: lo=mid
    dt=hi-t0; step=0.1 if dt<5 else (0.5 if dt<30 else 1.0); dt=max(2.0,np.ceil(dt/step-1e-9)*step) if v0==0 else np.ceil(dt/step-1e-9)*step
    for k in range(20):
        res={(m,u):run(m,u,v0,t0+dt) for m in (0.0,0.1) for u in (0.9999,)}
        if not any(np.isnan(r).any() for r in res.values()): return v0,dt,res,None
        dt+=step
    return v0,dt,res,"still missing directions"
if __name__=='__main__':
    with mp.Pool(6) as P:
        out={}
        for v0,dt,res,flag in P.imap_unordered(job,(0.0,0.1,0.3,0.5,0.9,0.99)):
            out[v0]=(dt,res,flag); pickle.dump(out,open('mcv_moving.pkl','wb'))
            def s(r):
                if r is None: return "-"
                if np.isnan(r).any(): return "missing %d"%np.isnan(r).sum()
                d=r/r.mean()-1; f=np.fft.rfft(d[:-1])/360; return "%.4f %.4f %.4f"%(2*abs(f[1]),2*abs(f[2]),d.min())
            print(v0,"dt",dt,flag,"" if res is None else " | ".join(f"m{m}: "+s(r) for (m,u),r in res.items()),flush=True)
