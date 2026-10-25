import sys,pickle; sys.path.insert(0,'.'); import numpy as np, iso, multiprocessing as mp
th=np.radians(np.arange(0,361,1.0))
H0=0.3; t0=2/(3*H0)   # assumption: observer at D=10R from the centre of expansion, so 3c = H0*D -> H0 R/c = 0.3
def run(c):
    m,u,dt=c; Mt=iso.Metric(m,t0); Xo0=np.array([2.0,0]); ems=[];hz=0
    for t in th:
        e=np.array([np.cos(t),np.sin(t)]); at,h=iso.worldline(Mt,Xo0+e,u*e,t0,t0+dt); ems.append((at,t0)); hz+=h
    oat,_=iso.worldline(Mt,Xo0,np.zeros(2),t0,t0+dt); Xo,vo=oat(t0+dt)
    r,info=iso.observe(Mt,Xo,vo,ems,t0+dt); return c,r,hz
if __name__=='__main__':
    dt=float(sys.argv[1])
    C=[(m,u,dt) for m in (0.0,0.1) for u in (0.0,0.9999)]
    with mp.Pool(4) as P: R={c:(r,h) for c,r,h in P.map(run,C)}
    pickle.dump(R,open(f'mcv_{dt}.pkl','wb'))
    st=lambda r:(lambda m:(2*abs(m[1]),2*abs(m[2]),(r/r.mean()-1).min()))(np.fft.rfft((r/r.mean()-1)[:-1])/360)
    for c,(r,h) in sorted(R.items()):
        print(c, "missing %d"%np.isnan(r).sum() if np.isnan(r).any() else "dip %.4f quad %.4f min %.4f"%st(r), "mean T_o/T_e %.4f"%np.nanmean(r),"hz",h,flush=True)
