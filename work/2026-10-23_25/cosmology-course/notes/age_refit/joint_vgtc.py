import json,numpy as np,importlib.util
from scipy.integrate import quad
from scipy.optimize import differential_evolution, minimize
src=open('../evidence.py').read().split('def chi2')[0]; ns={}; import os; os.chdir('..'); exec(src,ns); os.chdir('/workspace/cosmology-course/notes/age_refit')
Y,W,LB,CH0=ns['Y'],ns['W'],ns['LBIN'],ns['CHI2_NULL']
sp=importlib.util.spec_from_file_location('cc','../cmb-cancel.py'); cc=importlib.util.module_from_spec(sp); sp.loader.exec_module(cc)
C=299792.458; RCL=350.0  # Mpc, representative Migkas cluster distance (z~0.08)
def v2(r0,r,b1,mu,p): return b1**2*r0**(2*p)+2*mu*(1/r-1/r0)
def Smodel(b1,q,mu,s,phi,p):
    th=np.radians(cc.ELL-phi); xe0=1+q*np.cos(th); ye0=q*np.sin(th); re0=np.hypot(xe0,ye0)
    ro,re=s,s*re0; vo2=v2(1,ro,b1,mu,p); ve2=v2(re0,re,b1,mu,p)
    if vo2<=0 or np.any(ve2<=0) or 2*mu/ro>=1 or np.any(2*mu/re>=1): return None
    vo,ve=np.sqrt(vo2/(1+vo2)),np.sqrt(ve2/(1+ve2)); xe,ye=s*xe0,s*ye0; dx,dy=ro-xe,-ye
    a=vo**2-1;b=2*dx*vo;c=dx**2+dy**2; tf=(-b-np.sqrt(b*b-4*a*c))/(2*a); nx,ny=(dx+vo*tf)/tf,dy/tf
    Dm=np.sqrt(1+vo2)/np.sqrt(1+ve2)*(1-vo*nx)/(1-ve*(xe*nx+ye*ny)/re)
    return Dm*np.sqrt((1-2*mu/re)/(1-2*mu/ro))
def ageF(b1,mu,s):
    Uo2=v2(1,s,b1,mu,1); vo=np.sqrt(Uo2)
    return vo/s*quad(lambda r:np.sqrt(r)/np.sqrt(Uo2*r+2*mu*(1-r/s)),0,s,limit=100)[0]
def Hdir(ang,b1,mu,s,p,rl):  # finite-diff LOS H(ang), ang from CoE->us axis; units c/d1
    nx,ny=np.cos(ang),np.sin(ang); x,y=s+rl*nx,rl*ny; R=np.hypot(x,y)
    from scipy.optimize import brentq
    out=[]
    for Ri in np.atleast_1d(R):
        f=lambda r0:r0-Ri/s  # common expansion r=s r0
        out.append(Ri/s)
    r0=np.array(out); ve2=v2(r0,R,b1,mu,p)
    if np.any(ve2<=0): return None
    vo=np.sqrt(v2(1,s,b1,mu,p)); ve=np.sqrt(ve2)
    return (ve*(x*nx+y*ny)/R - vo*nx)/rl
GA=np.linspace(0,np.pi,181)
def evalp(x):
    b1,q,mu,s,phi,p,Hp=x
    S=Smodel(b1,q,mu,s,phi,p)
    if S is None: return None
    rms=cc.rms(cc.resid(cc.G_FULL,S)); F=ageF(b1,mu,s); age=F*977.8/Hp
    vo=np.sqrt(v2(1,s,b1,mu,1)); D=vo*C/Hp; d1=D/s; rl=RCL/d1; k=Hp/(vo/s)  # (km/s/Mpc) per (c/d1)
    Hb=Hdir(LB-np.radians(phi),b1,mu,s,p,rl); Hs=Hdir(GA,b1,mu,s,p,rl); Hc=Hdir(np.radians(cc.ELL-phi),b1,mu,s,p,rl)
    if Hb is None or Hs is None or Hc is None or np.any(Hc<=0): return None
    Hb,Hs,Hc=Hb*k,Hs*k,Hc*k
    A=Hb**2/np.mean(Hc**2); chi=float(np.sum(W*(A-Y)**2))
    Hm=float(np.sum(Hs*np.sin(GA))/np.sum(np.sin(GA)))
    return dict(vo=float(vo),rms=rms,age=age,chi2=chi,Hmean=Hm,D=D,M=mu*(C*1e3)**2*d1*3.0857e22/6.674e-11,A=A,Hc=Hc,
                Hamp=float((Hc.max()-Hc.min())/2/Hc.mean()),cancel=100*(1-rms/cc.rms(cc.resid(cc.G_FULL,np.ones(360)))))
B=[(1e-3,30),(1e-3,3),(0,2),(1,50),(0,360),(-1,3),(20,200)]
WT=dict(mig=1.0,cmb=10.0,age=100.0,H=100.0)
def cost(x):
    if any(v<lo or v>hi for v,(lo,hi) in zip(x,B)): return 1e6
    try: e=evalp(x)
    except Exception: return 1e6
    if e is None: return 1e6
    return (WT['mig']*e['chi2']/CH0 + WT['cmb']*max(0,e['rms']/2.37e-3-1)**2
            + WT['age']*max(0,(13-e['age'])/0.5)**2 + WT['H']*(max(0,67-e['Hmean'],e['Hmean']-73)/1.0)**2)
if __name__=='__main__':
    x0=[0.6106,0.2599,0.09395,3,15.17,0.903,55.54]; print('x0 cost',cost(x0))
    best=None
    for sd in (1,2):
        r=differential_evolution(cost,B,seed=sd,popsize=25,maxiter=600,tol=1e-10,polish=False,x0=x0)
        r2=minimize(cost,r.x,method='Nelder-Mead',options=dict(maxiter=20000,xatol=1e-10,fatol=1e-12,adaptive=True))
        x,f=(r2.x,r2.fun) if r2.fun<r.fun else (r.x,r.fun); print(sd,f,list(x),flush=True)
        if best is None or f<best[1]: best=(x,f)
    e=evalp(best[0]); json.dump(dict(x=list(best[0]),cost=best[1],**{k:(v.tolist() if hasattr(v,'tolist') else v) for k,v in e.items()}),open('joint_vgtc.json','w'),indent=1)
    print({k:v for k,v in e.items() if k not in('A','Hc')})
