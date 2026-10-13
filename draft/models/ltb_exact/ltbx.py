"""Exact LTB (GR dust) model, simultaneous bang t_B=0, softened central mass, off-centre observer.
Units: Gpc, Gyr. c=0.306601 Gpc/Gyr. Label r: E(r)=u^2/2, u=r/t0  (no mass -> Milne, H0=70 exact).
M(r)=Mc r^3/(r^2+rc^2)^1.5  (dust core, no singular point). Exact null geodesics with conserved J."""
import numpy as np
C=0.306601; KM=977.79; t0=KM/70.0
def eta_solve(s):
    e=np.where(s<1,np.cbrt(6*s),np.arcsinh(s)+0.5)
    for _ in range(80):
        f=np.where(e<1e-2,e**3/6*(1+e**2/20+e**4/840),np.sinh(e)-e)-s
        e=np.maximum(e-f/np.maximum(2*np.sinh(e/2)**2,1e-300),1e-12)
    return e
def RRd(r,t,GM):
    u=r/t0
    if np.all(GM==0): return u*t+0*r, u+0*r
    s=t*u**3/GM; e=eta_solve(s); ch=2*np.sinh(e/2)**2
    return GM/u**2*ch, u*np.sinh(e)/ch
class Model:
    def __init__(s,D,rs,rc=0.1,NQ=3000,NT=1500,rmax=None):
        s.D=D; GMc=rs*C**2/2; s.GMc=GMc
        rmax=rmax or 40*(D+16)
        s.r=np.geomspace(1e-4,rmax,NQ); s.GM=GMc*s.r**3/(s.r**2+rc**2)**1.5
        s.E2=(s.r/t0)**2/C**2   # 2E/c^2
        s.t=np.concatenate([np.geomspace(1e-5*t0,0.02*t0,NT//3,endpoint=False),np.linspace(0.02*t0,t0,NT-NT//3)])
        R=np.empty((len(s.t),NQ)); Rd=np.empty_like(R)
        for i,tt in enumerate(s.t): R[i],Rd[i]=RRd(s.r,tt,s.GM)
        s.R=R; s.Rd=Rd; s.Rp=np.gradient(R,s.r,axis=1); s.Rdp=np.gradient(Rd,s.r,axis=1)
        assert (s.Rp>0).all(), "shell crossing"
        s.ro=np.interp(D,R[-1],s.r)
    def f(s,t,rq):
        i=np.clip(np.searchsorted(s.t,t)-1,0,len(s.t)-2); w=(t-s.t[i])/(s.t[i+1]-s.t[i])
        out=[]
        for A in (s.R,s.Rd,s.Rp,s.Rdp):
            out.append(np.interp(rq,s.r,A[i])*(1-w)+np.interp(rq,s.r,A[i+1])*w)
        return out
    def trace(s,theta_deg,Rlist,sub=3,rmin=None):
        th=np.radians(theta_deg); n=len(th); rmin=rmin or 2*s.r[0]
        r=np.full(n,s.ro); ph=np.zeros(n); lnE=np.zeros(n); sg=np.where(np.cos(th)>=0,1.,-1.)
        J=s.D*np.sin(th); dist=np.zeros(n); dead=np.zeros(n,bool)
        out={R:np.full(n,np.nan) for R in Rlist}; tem={R:np.full(n,np.nan) for R in Rlist}; rem={R:np.full(n,np.nan) for R in Rlist}
        sq=lambda rq:1/np.sqrt(1+np.interp(rq,s.r,s.E2))
        R0=lambda rq:np.interp(rq,s.r,s.R[-1]); X0=lambda rq:np.interp(rq,s.r,s.Rp[-1])*sq(rq)
        ts=s.t[::-1]; minr=np.full(n,np.inf); s.maxd=np.zeros(n)
        def rhs(t,r,lnE,sg):
            R,Rd,Rp,Rdp=s.f(t,r); E=np.exp(lnE); sn=np.clip(J/(R*E),-1,1); cs2=1-sn**2
            X=Rp*sq(r)
            return sg*C*np.sqrt(cs2)/X, C*sn/R, Rdp/Rp*cs2+Rd/R*sn**2, cs2
        for k in range(len(ts)-1):
            for j in range(sub):
                ta=ts[k]+(ts[k+1]-ts[k])*j/sub; dt=(ts[k]-ts[k+1])/sub
                dr1,dp1,dl1,_=rhs(ta,r,lnE,sg)
                rm=np.maximum(r+.5*dt*dr1,rmin); dr2,dp2,dl2,cs2=rhs(ta-.5*dt,rm,lnE+.5*dt*dl1,sg)
                rn=r+dt*dr2; phn=ph+dt*dp2; lnEn=lnE+dt*dl2
                turn=(cs2<1e-10)&(sg<0); sg[turn]=1
                dead|=rn<=rmin; rn=np.maximum(rn,rmin)
                rmid=.5*(r+rn); seg=np.sqrt((X0(rmid)*(rn-r))**2+(R0(rmid)*(phn-ph))**2)
                dn=dist+seg; tn=ta-dt
                for Rr in Rlist:
                    hit=(dist<Rr)&(dn>=Rr)&np.isnan(out[Rr])&~dead; fr=(Rr-dist[hit])/np.maximum(seg[hit],1e-30)
                    out[Rr][hit]=lnE[hit]+fr*(lnEn[hit]-lnE[hit]); tem[Rr][hit]=tn; rem[Rr][hit]=np.interp(rn[hit],s.r,s.f(tn,s.r)[0]) if hit.any() else 0
                minr=np.where(np.isnan(out[max(Rlist)]),np.minimum(minr,rn),minr)
                s.maxd=np.where(dead,s.maxd,dn)
                r,ph,lnE,dist=rn,phn,lnEn,dn
            if all((~np.isnan(out[R])|dead).all() for R in Rlist): break
        H={R:C*out[R]/R*KM for R in Rlist}
        return H,tem,rem,minr
    def classical(s,theta_deg,R):
        th=np.radians(theta_deg); n=np.stack([np.cos(th),np.sin(th)],1)
        xe=np.array([s.D,0.])+R*n; re=np.linalg.norm(xe,axis=1)
        Rn=s.R[-1]; Vn=s.Rd[-1]; vo=np.interp(s.D,Rn,Vn); ve=np.interp(re,Rn,Vn)
        vl=((ve[:,None]*xe/re[:,None]-np.array([vo,0.]))*n).sum(1); return vl/R*KM
    def local(s,theta_deg):  # exact local limit: H_r cos^2 + H_perp sin^2 at observer
        R,Rd,Rp,Rdp=[np.interp(s.ro,s.r,A[-1]) for A in (s.R,s.Rd,s.Rp,s.Rdp)]
        th=np.radians(theta_deg); return (Rdp/Rp*np.cos(th)**2+Rd/R*np.sin(th)**2)*KM

def trace2(s,theta_deg,Rlist,sub=6):
    """Same as trace but integrates the radial momentum ODE (no sqrt/turning-point logic):
    dd_r/dtau = (Rdot'/R') d_r + c (R'/(X R)) (J/R)^2 / E ,  E = sqrt(d_r^2 + (J/R)^2)."""
    th=np.radians(theta_deg); n=len(th)
    r=np.full(n,s.ro); ph=np.zeros(n); dr_=np.cos(th).copy(); J=s.D*np.sin(th); dist=np.zeros(n); dead=np.zeros(n,bool)
    out={R:np.full(n,np.nan) for R in Rlist}
    sq=lambda rq:1/np.sqrt(1+np.interp(rq,s.r,s.E2))
    R0=lambda rq:np.interp(rq,s.r,s.R[-1]); X0=lambda rq:np.interp(rq,s.r,s.Rp[-1])*sq(rq)
    def rhs(t,r,d):
        R,Rd,Rp,Rdp=s.f(t,r); X=Rp*sq(r); pt=J/R; E=np.sqrt(d**2+pt**2)
        return C*d/(E*X), C*pt/(E*R), Rdp/Rp*d + C*(Rp/(X*R))*pt**2/E, E
    ts=s.t[::-1]; rmin=2*s.r[0]; lnE=np.zeros(n); s.maxd=np.zeros(n); s.edge=np.zeros(n,bool)
    for k in range(len(ts)-1):
        for j in range(sub):
            ta=ts[k]+(ts[k+1]-ts[k])*j/sub; dt=(ts[k]-ts[k+1])/sub
            a1,b1,c1,_=rhs(ta,r,dr_); rm=np.maximum(r+.5*dt*a1,rmin); a2,b2,c2,_=rhs(ta-.5*dt,rm,dr_+.5*dt*c1)
            rn=r+dt*a2; phn=ph+dt*b2; dn_=dr_+dt*c2
            dead|=rn<=rmin; rn=np.maximum(rn,rmin)
            Rn=s.f(ta-dt,rn)[0]; lnEn=np.log(np.sqrt(dn_**2+(J/Rn)**2))
            rmid=.5*(r+rn); seg=np.sqrt((X0(rmid)*(rn-r))**2+(R0(rmid)*(phn-ph))**2); dn=dist+seg
            for Rr in Rlist:
                hit=(dist<Rr)&(dn>=Rr)&np.isnan(out[Rr])&~dead; fr=(Rr-dist[hit])/np.maximum(seg[hit],1e-30)
                out[Rr][hit]=lnE[hit]+fr*(lnEn[hit]-lnE[hit])
            s.maxd=np.where(dead,s.maxd,dn); s.edge|=rn>=0.9*s.r[-1]
            r,ph,dr_,lnE,dist=rn,phn,dn_,lnEn,dn
        if all((~np.isnan(out[R])|dead).all() for R in Rlist): break
    return {R:C*out[R]/R*KM for R in Rlist}
Model.trace2=trace2
