import json,numpy as np
import importlib.util
sp=importlib.util.spec_from_file_location('cc','../cmb-cancel.py'); cc=importlib.util.module_from_spec(sp); sp.loader.exec_module(cc)
from scipy.integrate import quad
from scipy.optimize import differential_evolution, minimize
# v = proper (recession-style) speed U = dr/dtau, may exceed c. Doppler uses beta=U/sqrt(1+U^2), gamma=sqrt(1+U^2).
def S(U1,q,mu,s,phi):
    th=np.radians(cc.ELL-phi); xe0=1+q*np.cos(th); ye0=q*np.sin(th); re0=np.hypot(xe0,ye0)
    if np.any(re0<1e-6): return None
    ro,re=s,s*re0
    Uo2=U1**2+2*mu*(1/ro-1); Ue2=U1**2*re0**2+2*mu*(1/re-1/re0)
    if Uo2<=0 or np.any(Ue2<=0) or 2*mu/ro>=1 or np.any(2*mu/re>=1): return None
    bo=np.sqrt(Uo2/(1+Uo2)); be=np.sqrt(Ue2/(1+Ue2))
    xe,ye=s*xe0,s*ye0; dx,dy=ro-xe,-ye; a=bo**2-1; b=2*dx*bo; c=dx**2+dy**2
    tf=(-b-np.sqrt(b*b-4*a*c))/(2*a); nx,ny=(dx+bo*tf)/tf,dy/tf
    go,ge=np.sqrt(1+Uo2),np.sqrt(1+Ue2)
    Dm=go*(1-bo*nx)/(ge*(1-be*(xe*nx+ye*ny)/re)); Dg=np.sqrt((1-2*mu/re)/(1-2*mu/ro))
    return Dm*Dg
def obj(p):
    if any(v<lo or v>hi for v,(lo,hi) in zip(p,B)): return 1.0
    o=S(*p); return 1.0 if o is None else cc.rms(cc.resid(cc.G_FULL,o))
B=[(1e-3,20),(1e-4,3),(0,2),(1,3),(0,360)]
best=None
for sd in (1,2,3):
    r=differential_evolution(obj,B,seed=sd,popsize=30,maxiter=1500,tol=1e-12,polish=False,x0=[0.7255,0.2304,0.09327,3,15.16])
    r2=minimize(obj,r.x,method='Nelder-Mead',options=dict(xatol=1e-11,fatol=1e-15,maxiter=30000,adaptive=True))
    x,f=(r2.x,r2.fun) if r2.fun<r.fun else (r.x,r.fun); print(sd,f,x)
    if best is None or f<best[1]: best=(x,f)
U1,q,mu,s,phi=best[0]; Uo=np.sqrt(U1**2+2*mu*(1/s-1))
I=quad(lambda r:np.sqrt(r)/np.sqrt(Uo**2*r+2*mu*(1-r/s)),0,s,limit=200)[0]
ageHperp=Uo/s*I; Hperp=ageHperp*977.8/13; ratio=(U1**2+mu*(s-1)/s)/Uo**2; Hpar=ratio*Hperp
c=299792.458; D=Uo*c/Hperp; M=mu*(c*1e3)**2*(D/s)*3.0857e22/6.674e-11
S0=cc.rms(cc.resid(cc.G_FULL,np.ones(360)))
out=dict(U1=U1,beta1_kin=U1/np.sqrt(1+U1**2),q=q,mu=mu,s=s,phi=phi%360,vo_over_c=Uo,Hperp=Hperp,Hpar=Hpar,Hmean=Hpar/3+2*Hperp/3,age=13.0,rms=best[1],floor=best[1]/1e-5,cancel=100*(1-best[1]/S0),D_Mpc=D,M_kg=M)
print(json.dumps(out,indent=1)); json.dump(out,open('refit_vgtc.json','w'),indent=1)
