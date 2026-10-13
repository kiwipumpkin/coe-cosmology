import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import least_squares
exec(open("scan2.py").read().split("out={")[0])
KG=7.7e51/0.0039; floor=1e-5
amps=lambda x:np.array([np.hypot(2*np.mean(x*np.cos(k*th)),2*np.mean(x*np.sin(k*th))) for k in range(1,16)])
def model(q,N):
    k=q[:N]; a=q[N:2*N]; d=q[2*N:]
    return sum(k[i]*field(np.radians(a[i]),d[i]) for i in range(N))
out={}
for Mcap in (1e52,):
  for N in (12,):
    kc=Mcap/KG; best=None
    for s in range(4):
        rng=np.random.default_rng(s)
        q0=np.r_[rng.uniform(-kc,kc,N)*0.3,np.linspace(0,360,N,endpoint=False)+rng.uniform(-5,5,N),rng.uniform(1.15,2.5,N)]
        lb=np.r_[[-kc]*N,[-720]*N,[1.05]*N]; ub=np.r_[[kc]*N,[1080]*N,[5.0]*N]
        wnet=1.0/(kc*N)*np.sqrt(len(th))*1e-3  # soft zero-net-mass penalty
        fun=lambda q:np.r_[model(q,N)-R, wnet*np.sum(q[:N])*1e-0]
        f=least_squares(fun,q0,bounds=(lb,ub),max_nfev=600,x_scale=np.r_[[kc]*N,[30]*N,[0.5]*N])
        if best is None or f.cost<best.cost: best=f
    q=best.x; left=R-model(q,N); r=rms(left)
    out[f"{N}_{Mcap:.0e}"]=dict(N=N,Mcap=Mcap,rms=r,x_floor=r/floor,net_kg=float(q[:N].sum()*KG),masses=(q[:N]*KG).tolist(),ang=(q[N:2*N]%360).tolist(),dist_R=q[2*N:].tolist(),amps=(amps(left)/floor).tolist())
    print(N,f"{Mcap:.0e}",f"rms={r:.2e} ({r/floor:.1f}x)",f"net={q[:N].sum()*KG:.1e}","at cap:",int(np.sum(np.abs(q[:N])>0.98*kc)),flush=True)
json.dump(out,open("web12.json","w"),indent=1)
