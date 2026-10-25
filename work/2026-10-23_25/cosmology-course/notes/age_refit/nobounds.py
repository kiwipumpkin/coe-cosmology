import sys,json,numpy as np
from scipy.optimize import differential_evolution, minimize
import joint_vgtc as J
C=J.C; RPH=np.array([350.0])
DFIX=float(sys.argv[1]) if len(sys.argv)>1 else None
# params: log10 U1, log10 q, log10 mu, log10(s-1+1e-9)?, phi, p, log10 H
B=[(-3,3),(-4,2),(-12,2),(-2,4),(0,360),(-20,20),(0,4)]
if DFIX: B=B[:6]
def unpack(z):
    U1,q,mu,s=10**z[0],10**z[1],10**z[2],10**z[3]; phi,p=z[4],z[5]
    vo2=J.v2(1,s,U1,mu,1)
    if vo2<=0: return None
    H=(np.sqrt(vo2)*C/(DFIX*1e3)) if DFIX else 10**z[6]
    return [U1,q,mu,s,phi%360,p,H]
def ev(z):
    x=unpack(z)
    if x is None: return None
    b1,q,mu,s,phi,p,Hp=x
    try:
        S=J.Smodel(b1,q,mu,s,phi,p)
        if S is None: return None
        rms=J.cc.rms(J.cc.resid(J.cc.G_FULL,S)); age=J.ageF(b1,mu,s)*977.8/Hp
        vo=np.sqrt(J.v2(1,s,b1,mu,1)); D=vo*C/Hp; d1=D/s; k=Hp/(vo/s)
        Hb2=0;Hc2=0;Hs=0
        for r in RPH:
            hb=J.Hdir(J.LB-np.radians(phi),b1,mu,s,p,r/d1); hc=J.Hdir(np.radians(J.cc.ELL-phi),b1,mu,s,p,r/d1); hs=J.Hdir(J.GA,b1,mu,s,p,r/d1)
            if hb is None or hc is None or hs is None or not np.all(np.isfinite(hc)): return None
            Hb2=Hb2+(hb*k)**2/len(RPH); Hc2=Hc2+(hc*k)**2/len(RPH); Hs=Hs+hs*k/len(RPH)
        A=Hb2/np.mean(Hc2); chi=float(np.sum(J.W*(A-J.Y)**2))
        Hm=float(np.sum(Hs*np.sin(J.GA))/np.sum(np.sin(J.GA)))
    except Exception: return None
    if not np.isfinite(chi+rms+age): return None
    return dict(x=x,rms=rms,age=float(age),chi2=chi,Hmean=Hm,D=float(D),vo=float(vo),M=float(mu*(C*1e3)**2*d1*3.0857e22/6.674e-11),
                cancel=100*(1-rms/J.cc.rms(J.cc.resid(J.cc.G_FULL,np.ones(360)))),A=A.tolist())
def cost(z):
    e=ev(z)
    if e is None: return 1e6
    return (max(0,e['chi2']-11.6)/45.14*10 + e['chi2']/45.14 + 10*max(0,e['rms']/2.37e-3-1)**2
            +100*max(0,(13-e['age'])/0.5)**2 + 100*(max(0,67-e['Hmean'],e['Hmean']-73))**2)
best=None
for sd in (1,2):
    r=differential_evolution(cost,B,seed=sd,popsize=30,maxiter=700,tol=1e-10,polish=False)
    r2=minimize(cost,r.x,method='Nelder-Mead',options=dict(maxiter=20000,xatol=1e-10,fatol=1e-12,adaptive=True))
    z,f=(r2.x,r2.fun) if r2.fun<r.fun else (r.x,r.fun); print(sd,f,list(z),flush=True)
    if best is None or f<best[1]: best=(z,f)
z=best[0]; e=ev(z)
edges=[i for i,(v,(lo,hi)) in enumerate(zip(z,B)) if i!=4 and (abs(v-lo)<1e-3*(hi-lo) or abs(v-hi)<1e-3*(hi-lo))]
e.update(cost=best[1],z=list(z),edges=edges,DFIX=DFIX)
json.dump(e,open(f'nobounds_{DFIX or "free"}.json','w'),indent=1,default=float)
print({k:v for k,v in e.items() if k!='A'})
