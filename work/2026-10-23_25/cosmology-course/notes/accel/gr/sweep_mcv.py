import sys,pickle; sys.path.insert(0,'.'); import numpy as np, iso, multiprocessing as mp
th=np.radians(np.arange(0,361,1.0))
def one(m,u,H0,dt):
    t0=2/(3*H0); Mt=iso.Metric(m,t0); Xo0=np.array([2.0,0]); ems=[]
    for t in th:
        e=np.array([np.cos(t),np.sin(t)]); at,h=iso.worldline(Mt,Xo0+e,u*e,t0,t0+dt); ems.append((at,t0))
    oat,_=iso.worldline(Mt,Xo0,np.zeros(2),t0,t0+dt); Xo,vo=oat(t0+dt)
    return iso.observe(Mt,Xo,vo,ems,t0+dt)[0]
def job(v):
    H0=v/10; dt=2.0
    while True:
        res={(m,u):one(m,u,H0,dt) for m in (0.0,0.1) for u in (0.0,0.9999)}
        if not any(np.isnan(r).any() for r in res.values()) or dt>200: return v,dt,res
        dt*=1.25
if __name__=='__main__':
    with mp.Pool(6) as P:
        out={}
        for v,dt,res in P.imap_unordered(job,(1,1.5,2,3,5,10)):
            out[v]=(dt,res); pickle.dump(out,open('sweep_mcv.pkl','wb'))
            def s(r):
                if np.isnan(r).any(): return "missing %d"%np.isnan(r).sum()
                d=r/r.mean()-1; f=np.fft.rfft(d[:-1])/360; return "%.4f %.4f %.4f"%(2*abs(f[1]),2*abs(f[2]),d.min())
            print(v,"dt",dt," | ".join(f"m{m} u{u}: "+s(r) for (m,u),r in res.items()),flush=True)
