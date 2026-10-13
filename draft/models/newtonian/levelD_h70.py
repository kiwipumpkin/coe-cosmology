"""Level sets over observer distance D for H=70 runs (expanding: integ_h70.py; collapsing: integ_h70_collapse.py). Same physical masses at every D."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
c=2.998e8; GPC=3.0857e25; KMSMPC=299792.458/1000
b1,mu=0.5866,0.09332; GM1=mu*c**2*(58.2/3)*GPC; k0=c*np.sqrt(b1**2+2*mu*(1/3-1))/58.2/GPC
kp=3*k0; Tp=0.247*(58.2*GPC)**2*k0*14.3*GPC/GM1; kTp=kp*Tp; H=70/KMSMPC*c/GPC
Mscale=(H*(1+kTp)/kp)**2*(3.2/58.2)**3
MODES={"expand":(H*(1+kTp),kTp/(H*(1+kTp)),+1,np.linspace(0.02,5,4000)),"collapse":(H*0.9,0.1/(H*0.9),-1,np.linspace(0.05,14,6000))}
DS=[0.8,1.6,2.4,3.2]; Ms=[0,1,4,10]
ang=np.linspace(0,360,1441); t=np.radians(ang); n=np.stack([np.cos(t),np.sin(t)],1); w=np.abs(np.sin(t))
def field(mode,m):
    k,T,s,g=MODES[mode]; ri=g*GPC; GM=m*Mscale*GM1; r=ri.copy(); v=s*k*ri; N=4000; h=T/N; f=lambda r:-GM/np.maximum(r,1e18)**2
    for _ in range(N):
        a1,q1=v,f(r); a2,q2=v+.5*h*q1,f(r+.5*h*a1); a3,q3=v+.5*h*q2,f(r+.5*h*a2); a4,q4=v+h*q3,f(r+h*a3)
        r=r+h/6*(a1+2*a2+2*a3+a4); v=v+h/6*(q1+2*q2+2*q3+q4)
    return r/GPC,v
plt.rcParams.update({"font.size":17,"axes.titlesize":18,"axes.labelsize":18,"legend.fontsize":14,"xtick.labelsize":16,"ytick.labelsize":16})
cm=plt.cm.plasma(np.linspace(0,.85,len(DS)))
for mode in MODES:
    F={}; s=MODES[mode][2]; OK={}
    for m in Ms:
        rf,v=field(mode,m)
        for D in DS:
            sel=(rf>D*0.754-0.01)&(rf<D*1.246+0.01); i1=np.where(sel)[0].max()+1
            ok=np.all(np.diff(rf[:i1+1])>0) and np.all(s*v[sel]>0) and np.all(abs(v[sel])<c)
            OK[(m,D)]=ok; print(mode,"M=%g D=%.1f ok=%s vrange %.3f..%.3fc"%(m,D,ok,v[sel].min()/c,v[sel].max()/c))
        F[m]=(lambda rf,v:(lambda x:np.interp(x,rf,v)))(rf,v)
    def vlos(m,D,Rv,nn):
        xo=np.array([D,0.]); xp=xo+Rv*nn; rp=np.linalg.norm(xp,axis=-1); vr=F[m](rp)[...,None]*xp/rp[...,None]-F[m](np.array(D))*xo/D; return (vr*nn).sum(-1)/1e3
    T="Expanding" if s>0 else "Collapsing"
    for m in Ms[1:]:
        for a in (180,90):
            nn=np.array([np.cos(np.radians(a)),np.sin(np.radians(a))]); f,ax=plt.subplots(figsize=(8,6))
            for D,cl in zip(DS,cm):
                if not OK[(m,D)]: continue
                Rs=np.linspace(0.005,0.246*D,100); ax.plot(Rs*1000,[vlos(m,D,r,nn) for r in Rs],color=cl,lw=3.5,label="D=%.1f Gpc"%D)
            Rs=np.linspace(0,0.246*DS[-1],50); ax.plot(Rs*1000,s*70*Rs*1000,"k--",lw=2,label="no mass")
            ax.set_xlabel("galaxy distance R (Mpc)"); ax.set_ylabel("v_LOS (km/s)"); ax.set_xlim(0,None)
            ax.set_title("%s, M = %g×, angle %d°"%(T,m,a)); ax.legend(frameon=False,loc="upper left" if s>0 else "lower left")
            f.tight_layout(); f.savefig("h70_%s_levelD_M%g_%d.png"%(mode,m,a),dpi=110); plt.close(f)
    res={m:[] for m in Ms}
    for m in Ms:
        for D in DS:
            Rv=0.246*D; hh=0.002; y=(vlos(m,D,Rv+hh,n)-vlos(m,D,Rv-hh,n))/(2*hh)/1000
            mn=(w*y).sum()/w.sum(); sd=np.sqrt((w*(y-mn)**2).sum()/w.sum()); pos=100*(w*(y>0)).sum()/w.sum()
            res[m].append((mn,sd) if OK[(m,D)] else (np.nan,np.nan)); print("%s M=%g D=%.1f Rmax=%.0f: 0° %.1f 90° %.1f 180° %.1f mean %.1f sd %.1f pos %.0f%%"%(mode,m,D,Rv*1000,y[0],y[360],y[720],mn,sd,pos))
    for j,(lab,tag) in enumerate((("sky-mean dv/dR (km/s/Mpc)","mean"),("sky std dev of dv/dR (km/s/Mpc)","std"))):
        f,ax=plt.subplots(figsize=(8,6)); cc=plt.cm.viridis(np.linspace(.1,.9,3)); i=0
        for m in Ms:
            yy=[r[j] for r in res[m]]
            if m==0: ax.plot(DS,yy,"k--",lw=2.5,label="no mass")
            else: ax.plot(DS,yy,"o-",color=cc[i],lw=3.5,ms=9,label="%g×"%m); i+=1
        ax.set_xlabel("observer distance from CoE D (Gpc)"); ax.set_ylabel(lab); ax.set_title("%s, at R = 0.246 D"%T)
        ax.legend(frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(.5,-.18)); f.tight_layout(); f.savefig("h70_%s_levelD_sky%s.png"%(mode,tag),dpi=110); plt.close(f)
