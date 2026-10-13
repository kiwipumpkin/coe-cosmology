"""Independent cosmological redshift: exact LTB (relativistic) version of the same shell model.
Metric ds^2=-c^2dt^2 + R'^2/(1+2E/c^2) dq^2 + R^2 dOmega^2, R(q,t) from the same Newtonian shell eq (exact in LTB).
Null geodesics in label space; 1+z from photon energy measured by comoving observers."""
import numpy as np, sys, pickle
exec(open('integ_h70D.py').read().split('ang=')[0])
kT40=199.0; k=H*(1+kT40); T=kT40/k
A=np.arange(0,181,1.0); th=np.radians(A); w=np.sin(th)
def model(D,m,N=20000,S=20,NQ=4001,flat=False):
    Msc=(k/kp)**2*(D/58.2)**3*((1+kTp)/(1+kT40))**3; GM=m*Msc*GM1
    q=np.linspace(1e-4,4.0,NQ)*D*GPC/(1+kT40); h=T/N
    r=q.copy(); v=k*q; f=lambda r:-GM/np.maximum(np.abs(r),1e18)**2; sR=[r.copy()]; sV=[v.copy()]
    for i in range(N):
        a1,q1=v,f(r); a2,q2=v+.5*h*q1,f(r+.5*h*a1); a3,q3=v+.5*h*q2,f(r+.5*h*a2); a4,q4=v+h*q3,f(r+h*a3)
        r=r+h/6*(a1+2*a2+2*a3+a4); v=v+h/6*(q1+2*q2+2*q3+q4)
        if (i+1)%S==0: sR.append(r.copy()); sV.append(v.copy())
    sR=np.array(sR); sV=np.array(sV); ts=np.arange(len(sR))*h*S
    E2=(k*q)**2/c**2-2*GM/(q*c**2)          # 2E/c^2
    g=np.ones_like(q) if flat else 1/np.sqrt(1+E2)   # A = R' * g
    sRp=np.gradient(sR,q,axis=1); sVp=np.gradient(sV,q,axis=1)
    return dict(q=q,sR=sR,sV=sV,sRp=sRp,sVp=sVp,g=g,ts=ts,D=D)
def fields(M,t,qq):
    j=np.clip(t/M['ts'][1],0,len(M['ts'])-1.000001); a=int(j); fr=j-a; o={}
    for key in ('sR','sV','sRp','sVp'):
        X=M[key][a]*(1-fr)+M[key][a+1]*fr; o[key]=np.interp(qq,M['q'],X)
    o['g']=np.interp(qq,M['q'],M['g']); return o
def trace(M,Rlist,sub=4):
    q0=np.interp(M['D']*GPC,M['sR'][-1],M['q']); nt=(len(M['ts'])-1)*sub; dt=M['ts'][-1]/nt
    qq=np.full(len(A),q0); ph=np.zeros(len(A)); lnp=np.zeros(len(A)); s=np.sign(np.cos(th)); s[s==0]=1
    L=M['D']*GPC*np.sin(th)    # p=1 at obs, p_t*R conserved
    dist=np.zeros(len(A)); out={R:np.full(len(A),np.nan) for R in Rlist}; zmax=np.zeros(len(A))
    fin=fields(M,M['ts'][-1],qq)
    def rhs(t,qq,lnp,s):
        F=fields(M,t,qq); p=np.exp(lnp); pt=L/F['sR']; pr2=np.maximum(p**2-pt**2,0); pr=s*np.sqrt(pr2)
        Aq=F['sRp']*F['g']; cos2=pr2/p**2; sin2=1-cos2
        dq=c*pr/(p*Aq); dph=c*pt/(p*F['sR']); dl=(F['sVp']/F['sRp'])*cos2+(F['sV']/F['sR'])*sin2
        return dq,dph,dl,pr2,p,pt
    for i in range(nt):
        t=M['ts'][-1]-i*dt
        dq,dph,dl,pr2,p,pt=rhs(t,qq,lnp,s)
        qm=qq+.5*dt*dq; lm=lnp+.5*dt*dl
        dq2,dph2,dl2,pr2m,_,_=rhs(t-.5*dt,qm,lm,s)
        qn=qq+dt*dq2; phn=ph+dt*dph2; lnp_n=lnp+dt*dl2
        # present-slice length of the label path
        F0=fields(M,M['ts'][-1],.5*(qq+qn)); seg=np.sqrt((F0['sRp']*F0['g']*(qn-qq))**2+(F0['sR']*(phn-ph))**2)
        dn=dist+seg
        for R in Rlist:
            hit=(dist<R*GPC)&(dn>=R*GPC)&np.isnan(out[R]); fr=(R*GPC-dist[hit])/seg[hit]
            out[R][hit]=lnp[hit]+fr*(lnp_n[hit]-lnp[hit])
        qq,ph,lnp,dist=qn,phn,lnp_n,dn
        turn=(pr2m<=1e-12)&(s<0); s[turn]=1   # radial turning point
        if all(~np.isnan(out[R]).any() for R in Rlist): break
    return {R:c/1e3*out[R]/(R*1e3) for R in Rlist}   # km/s/Mpc
if __name__=="__main__":
    RL=[3.0,6.0,10.0,14.3]; RES={}
    flat='flat' in sys.argv
    for D in [int(x) for x in sys.argv[1].split(',')]:
        for m in [float(x) for x in sys.argv[2].split(',')]:
            M=model(D,m,flat=flat); Hs=trace(M,RL); RES[(D,m)]=Hs
            for R in RL:
                y=Hs[R]; mn=np.nansum(w*y)/w[~np.isnan(y)].sum()
                print("%s D=%d M=%gx R=%.1f  0° %.1f 45° %.1f 90° %.1f 135° %.1f 180° %.1f mean %.1f | A/Aall 0 %.3f 90 %.3f 180 %.3f"%("FLAT" if flat else "LTB",D,m,R,y[0],y[45],y[90],y[135],y[180],mn,y[0]/mn,y[90]/mn,y[180]/mn),flush=True)
    pickle.dump(RES,open("ltb14_res%s_%s_%s.pkl"%("_flat" if flat else "",sys.argv[1],sys.argv[2]),"wb"))
