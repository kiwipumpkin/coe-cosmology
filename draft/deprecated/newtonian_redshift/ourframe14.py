"""Fixed physical central mass, shells launched together from the CoE (hyperbolic Kepler, exact), t0=1/70.
Cosmological z traced in OUR frame (we at rest at origin, matter velocity relative to us), also CoE frame for comparison.
Classical: present relative v_los/R. Relativistic: rapidity Doppler (CoE-frame speeds)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, pickle, sys
c=2.998e8; GPC=3.0857e25; G=6.674e-11; Msun=1.989e30
H0=70e3/3.0857e22; t0=1/H0
A=np.arange(0,181,1.0); th=np.radians(A); w=np.sin(th); n=np.stack([np.cos(th),np.sin(th)],1)
def eta_of(s):
    e=np.where(s<1,np.cbrt(6*s),np.arcsinh(s)+0.5)
    for _ in range(60): e=e-(np.sinh(e)-e-s)/np.maximum(np.cosh(e)-1,1e-30); e=np.maximum(e,1e-12)
    return e
def shells(GM,u,t):  # u: asymptotic speed label; returns r,v
    t=np.asarray(t)[...,None]
    if GM==0: return u*t+0*u, u+0*t
    s=t*u**3/GM; e=eta_of(s); r=GM/u**2*(np.cosh(e)-1); v=u*np.sinh(e)/(np.cosh(e)-1)
    small=s<1e-6; r=np.where(small,(4.5*GM*t**2)**(1/3)+0*u,r); v=np.where(small,(2/3)*r/t,v); return r,v
def run(D,rs,NS=2000,sub=6,RP=14.3):
    GM=rs*GPC*c**2/2
    umax=3*(D+18)*GPC/t0+3*np.sqrt(GM/(0.1*GPC)); u=np.geomspace(1e-5,1,4000)*umax
    ts=np.linspace(t0*1e-4,t0,NS); R,V=shells(GM,u,ts); dVdR=np.nan_to_num(np.gradient(V,axis=1)/np.gradient(R,axis=1),nan=0.0,posinf=0.0,neginf=0.0)
    uo=np.interp(D*GPC,R[-1],u); Ro,Vo=shells(GM,np.array([uo]),ts); Ro=Ro[:,0]; Vo=Vo[:,0]
    def fld(t,rr):
        j=np.clip((t-ts[0])/(ts[1]-ts[0]),0,NS-1.0001); a=int(j); f=j-a
        Rr=R[a]*(1-f)+R[a+1]*f; Vv=V[a]*(1-f)+V[a+1]*f; Gg=dVdR[a]*(1-f)+dVdR[a+1]*f
        X=Ro[a]*(1-f)+Ro[a+1]*f; U=Vo[a]*(1-f)+Vo[a+1]*f
        return np.interp(rr,Rr,Vv),np.interp(rr,Rr,Gg),Rr,X,U
    res={}
    for frame in ("ours","coe"):
        MINR=np.full(len(A),1e40); xp=np.zeros((len(A),2)); lnz=np.zeros(len(A)); out=np.full(len(A),np.nan)
        nt=(NS-1)*sub; dt=(ts[-1]-ts[0])/nt; prevd=np.zeros(len(A)); prevz=np.zeros(len(A))
        for i in range(nt):
            t=ts[-1]-i*dt
            def step(t,xp):
                # absolute position from our-frame position
                _,_,_,X,U=fld(t,np.ones(1))
                xa=xp+np.array([X,0.]) if frame=="ours" else xp
                rr=np.maximum(np.linalg.norm(xa,axis=1),1e-4*GPC); vr,dv,_,_,_=fld(t,rr); ua=vr[:,None]*xa/rr[:,None]
                urel=ua-np.array([U,0.]) if frame=="ours" else ua      # matter velocity in chosen frame
                cs=(xa*n).sum(1)/rr; grad=dv*cs**2+vr/rr*(1-cs**2)     # n.grad(u).n (same in both frames)
                return (n*c-urel), grad
            if i==0 and frame=="coe": xp=np.tile([Ro[-1],0.],(len(A),1))
            d1,g1=step(t,xp); xm=xp+.5*dt*d1; d2,g2=step(t-.5*dt,xm)
            xp=xp+dt*d2; lnz=lnz+dt*g2
            # emitter present position: label at its current radius, moved to now along same radial direction
            tt=t-dt; _,_,Rr,X,_=fld(tt,np.ones(1))
            xa=xp+np.array([X,0.]) if frame=="ours" else xp; rr=np.linalg.norm(xa,axis=1); MINR=np.where(np.isnan(out),np.minimum(MINR,rr),MINR)
            lab=np.interp(rr,Rr,u); xnow=xa/rr[:,None]*np.interp(lab,u,R[-1])[:,None]
            dnow=np.linalg.norm(xnow-np.array([Ro[-1],0.]),axis=1)/GPC
            hit=(prevd<RP)&(dnow>=RP)&np.isnan(out); fr=(RP-prevd[hit])/(dnow[hit]-prevd[hit]); out[hit]=prevz[hit]+fr*(lnz[hit]-prevz[hit])
            prevd=dnow; prevz=lnz.copy()
            if not np.isnan(out).any(): break
        res[frame]=c/1e3*out/(RP*1e3)
    # classical & relativistic (present day)
    xe=np.array([D*GPC,0.])+RP*GPC*n; re=np.linalg.norm(xe,axis=1); ve=np.interp(re,R[-1],V[-1]); rh=xe/re[:,None]
    vlos=((ve[:,None]*rh-np.array([Vo[-1],0.]))*n).sum(1); res["classical"]=vlos/(RP*1e3)/1e3
    ee=ve/c; eo=Vo[-1]/c; opz=(np.cosh(ee)+np.sinh(ee)*(rh*n).sum(1))/(np.cosh(eo)+np.sinh(eo)*n[:,0]); res["relativistic"]=c/1e3*np.log(opz)/(RP*1e3)
    res["minr"]=MINR/GPC; res["Hlocal"]=Vo[-1]/Ro[-1]*3.0857e22/1e3; res["vus"]=Vo[-1]/c; res["M"]=GM/G
    return res
if __name__=="__main__":
    D=int(sys.argv[1]); RES={}
    for rs in [float(x) for x in sys.argv[2].split(',')]:
        r_=run(D,rs); RES[rs]=r_
        for key in ("ours","coe","classical","relativistic"):
            y=r_[key]; mn=(w*y).sum()/w.sum(); q=y/mn
            print("D=%d rs=%g Gpc M=%.2e kg (%.1e Msun) us v=%.2fc v/D=%.1f | %-12s H 0:%.1f 45:%.1f 90:%.1f 135:%.1f 180:%.1f mean %.1f | A/Aall %.3f %.3f %.3f %.3f %.3f"%(D,rs,r_["M"],r_["M"]/Msun,r_["vus"],r_["Hlocal"],key,y[0],y[45],y[90],y[135],y[180],mn,*q[[0,45,90,135,180]]),flush=True)
    pickle.dump(RES,open("ourframe14_D%d_%s.pkl"%(D,sys.argv[2]),"wb"))
