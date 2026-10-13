import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import minimize
G=6.674e-11; c=2.998e5; Mpc=3.0857e22; k=70/Mpc*1e3  # 1/s
D0=3200*Mpc; mu1=G*6.08e51/(k**2*D0**3); vu=70*3200/c   # dimless vel unit -> c at D=3.2Gpc
ang=np.linspace(0,360,1441); t=np.radians(ang); n=np.stack([np.cos(t),np.sin(t)],1); w=np.abs(np.sin(t))
def fieldA(mu,tau):  # start v=-K r0, run T; units D=1, no-mass final slope -1
    K=1-tau; T=tau/(1-tau); r=np.linspace(0.02,3.0/K,3000); v=-K*r; N=1500; h=T/N; alive=np.ones_like(r,bool)
    f=lambda r,v:(v,-mu/np.maximum(r,1e-4)**2)
    for _ in range(N):
        a1=f(r,v);a2=f(r+h/2*a1[0],v+h/2*a1[1]);a3=f(r+h/2*a2[0],v+h/2*a2[1]);a4=f(r+h*a3[0],v+h*a3[1])
        r=r+h/6*(a1[0]+2*a2[0]+2*a3[0]+a4[0]); v=v+h/6*(a1[1]+2*a2[1]+2*a3[1]+a4[1]); alive&=r>1e-3
    return r,v,alive
def fieldB(rs):  # stall: v=-sqrt(2mu(1/r-1/rs)); mu cancels in sign; scale so v(D)=-1
    r=np.linspace(0.01,rs*0.9999,4000); v=-np.sqrt(1/r-1/rs); return r,v/np.sqrt(1-1/rs),np.ones_like(r,bool)
def evalf(r,v,alive,RD):
    lo,hi=1-RD-.01,1+RD+.01; sel=(r>lo)&(r<hi)
    if r[alive].min()>lo or r.max()<hi: return None
    i=np.where(sel)[0]; ok=np.all(v[sel]<0) and np.all(np.diff(r[i[0]-1:i[-1]+2])>0) and alive[sel].all()
    F=lambda x:np.interp(x,r,v); xo=np.array([1.,0])
    def vl(R): xp=xo+R*n; rp=np.linalg.norm(xp,axis=1); return ((F(rp)[:,None]*xp/rp[:,None]-F(1)*xo)*n).sum(1)
    out={}
    for nm,R in (("max",RD),("small",0.1*RD)):
        e=1e-3*RD; y=(vl(R+e)-vl(R-e))/(2*e); mn=(w*y).sum()/w.sum()
        out[nm]=dict(y=y,win=100*(w*(y>0)).sum()/w.sum(),plain=100*(y>0).mean(),mean=mn,sd=np.sqrt((w*(y-mn)**2).sum()/w.sum()))
    out["ok"]=ok; out["vmax"]=np.abs(v[sel]).max()*vu; return out
# --- family A grid
mus=np.linspace(0,2.0,26); taus=np.linspace(0,0.9,19); RDs=[0.1,0.246,0.4,0.6]
res={}
for RD in RDs:
  W=np.full((len(taus),len(mus)),np.nan); V=W.copy()
  for i,ta in enumerate(taus):
    for j,mu in enumerate(mus):
      fa=fieldA(mu,ta); o=evalf(*fa,RD)
      if o and o["ok"]: W[i,j]=o["max"]["win"]; V[i,j]=o["vmax"]
  res[RD]=(W,V); print("A RD",RD,"best",np.nanmax(W),flush=True)
np.save("maxwin_A.npy",dict(mus=mus,taus=taus,res=res),allow_pickle=True)
# --- family B grid
rss=np.linspace(1.05,6,50); RDb=np.linspace(0.02,0.95,48); WB=np.full((len(rss),len(RDb)),np.nan)
for i,rs in enumerate(rss):
  fb=fieldB(rs)
  for j,RD in enumerate(RDb):
    if rs<=1+RD+.02: continue
    o=evalf(*fb,RD); WB[i,j]=o["max"]["win"] if o and o["ok"] else np.nan
np.save("maxwin_B.npy",dict(rss=rss,RDb=RDb,WB=WB),allow_pickle=True); print("B best",np.nanmax(WB))
