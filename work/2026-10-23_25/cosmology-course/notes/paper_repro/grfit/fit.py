import sys
from core import *
from scipy.optimize import differential_evolution, minimize
model,sn=sys.argv[1],sys.argv[2]=="sine"
fn=ltb if model=="ltb" else nsky
B=[(np.log10(20),np.log10(150)),(np.log10(1.7e-4),np.log10(1.7e-2)),(-3,np.log10(9)),(-3,1),(0,360)]
def par(u): return 10**u[0],10**u[1],1+10**u[2],1+10**u[3],u[4]%360
def ev(u):
    D,mcc,kin,kout,l0=par(u); af=fn(D,mcc,kin,kout)
    if af is None: return 800.0
    c=shifted(af,l0)
    if sn: c=sine(c,l0)[0]
    return chi(c)
t0=time.time(); rng=np.random.default_rng(1)
de=differential_evolution(ev,B,seed=1,popsize=12 if model=="ltb" else 8,maxiter=200 if model=="ltb" else 30,tol=1e-7,polish=False)
mins=[]
starts=[de.x]+[np.array([rng.uniform(*b) for b in B]) for _ in range(23 if model=="ltb" else 19)]
for s in starts:
    r=minimize(ev,s,method="Nelder-Mead",options=dict(maxiter=1500 if model=="ltb" else 250,xatol=1e-5,fatol=1e-6,adaptive=True))
    mins.append((float(r.fun),r.x.tolist()))
best=min(m[0] for m in mins+[(float(de.fun),de.x.tolist())])
rows=[]
for f,u in mins+[(float(de.fun),de.x.tolist())]:
    if f<best+1:
        D,mcc,kin,kout,l0=par(np.array(u)); af=fn(D,mcc,kin,kout); c=shifted(af,l0); extra={}
        if sn: c,a,ph=sine(c,l0); extra=dict(amp=float(a),phi_deg=float(np.degrees(ph)))
        rows.append(dict(D=D,mcc=mcc,kin=kin,kout=kout,l0=l0,Rin=D-kin*RV,Rout=D+kout*RV,chi2=chi(c),R2=r2(c),**extra))
json.dump(dict(model=model,sine=sn,de=float(de.fun),nm=[m[0] for m in mins],rows=rows,seconds=time.time()-t0),open(f"fit_{model}_{sys.argv[2]}.json","w"),indent=1)
print(model,sys.argv[2],"done",time.time()-t0,best)
