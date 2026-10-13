# REAL geometry D=58.2 Gpc, R<=14.3 Gpc; v=H r is the only isotropic no-mass field
"""H0=70 tuned expanding run = exact dimensionless rescaling of integ v2 (k=3k0, kT=1.286, R/D=0.246, GM/(k^2 D^3) matched per mass level)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
c=2.998e8; GPC=3.0857e25; GYR=3.156e16; KMSMPC=299792.458/1000
b1,mu=0.5866,0.09332; GM1=mu*c**2*(58.2/3)*GPC; k0=c*np.sqrt(b1**2+2*mu*(1/3-1))/58.2/GPC
kp=3*k0; Tp=0.247*(58.2*GPC)**2*k0*14.3*GPC/GM1; Dp=58.2; kTp=kp*Tp
H=70/KMSMPC*c/GPC                       # 1/s
k=H*(1+kTp); T=kTp/k; D=58.2; RMAX=14.3
Mscale=(k/kp)**2*(D/Dp)**3              # GM_new/GM_prev for same dimensionless strength
print("H %.4e 1/s, k %.4f c/Gpc (%.1f km/s/Mpc), T %.2f Gyr, D %.2f Gpc, Rmax %.3f Gpc, mass scale %.3e"%(H,k*GPC/c,k*GPC/c*KMSMPC,T/GYR,D,RMAX,Mscale))
ang=np.linspace(0,360,1441); t=np.radians(ang); n=np.stack([np.cos(t),np.sin(t)],1); w=np.abs(np.sin(t))
def field(m):
    ri=np.linspace(0.3,90,3000)*GPC; GM=m*Mscale*GM1; r=ri.copy(); v=k*ri; N=4000; h=T/N; f=lambda r:-GM/np.maximum(r,1e18)**2
    for _ in range(N):
        a1,b1_=v,f(r); a2,b2=v+.5*h*b1_,f(r+.5*h*a1); a3,b3=v+.5*h*b2,f(r+.5*h*a2); a4,b4=v+h*b3,f(r+h*a3)
        r=r+h/6*(a1+2*a2+2*a3+a4); v=v+h/6*(b1_+2*b2+2*b3+b4)
    rf=r/GPC; assert np.all(np.diff(rf)>0); return lambda x:np.interp(x,rf,v)
Ms=[0,1,4,10]; F={m:field(m) for m in Ms}; xo=np.array([D,0.])
for m in Ms:
    sv=F[m](np.linspace(D-RMAX,D+RMAX,300))/c; assert sv.min()>0; print("M=%g (=%.3e kg): our v %.3fc, shells in view %.3f..%.3fc"%(m,m*Mscale*7.5e52,F[m](D)/c,sv.min(),sv.max()))
def vlos(m,Rv,nn):
    xp=xo+Rv*nn; rp=np.linalg.norm(xp,axis=-1); vr=F[m](rp)[...,None]*xp/rp[...,None]-F[m](np.array(D))*xo/D; return (vr*nn).sum(-1)/c
plt.rcParams.update({"font.size":17,"axes.titlesize":18,"axes.labelsize":18,"legend.fontsize":15,"xtick.labelsize":16,"ytick.labelsize":16})
Rs=np.linspace(0.01,RMAX,120); dirs=[0,45,90,135,180]; cols=["#b2182b","#ef8a62","#4d4d4d","#67a9cf","#2166ac"]
vmax=max(vlos(m,RMAX,np.array([np.cos(np.radians(a)),np.sin(np.radians(a))])) for m in Ms for a in dirs)*c/1e3
for m in Ms:
    f,ax=plt.subplots(figsize=(8,6))
    for a_,cl in zip(dirs,cols):
        nn=np.array([np.cos(np.radians(a_)),np.sin(np.radians(a_))]); ax.plot(Rs*1000,[vlos(m,r,nn)*c/1e3 for r in Rs],color=cl,lw=3.5,label="%d°"%a_)
    ax.plot(Rs*1000,[vlos(0,r,np.array([1.,0]))*c/1e3 for r in Rs],"k--",lw=2,label="no mass")
    ax.set_xlabel("galaxy distance R (Mpc)"); ax.ticklabel_format(style="sci",scilimits=(0,0)); ax.set_ylabel("v_LOS (km/s)"); ax.set_ylim(0,vmax*1.05); ax.set_xlim(0,RMAX*1000)
    ax.set_title("Expanding, H₀ = 70, "+("no mass" if m==0 else "M = %g×"%m)); ax.legend(frameon=False,loc="upper left")
    f.tight_layout(); f.savefig("integ_h70real_diagram_M%g.png"%m,dpi=110); plt.close(f)
hh=0.002; mean=lambda y:(w*y).sum()/w.sum(); cm=plt.cm.viridis(np.linspace(.1,.9,3))
Rsmall=1.0
for Rv,tag in ((Rsmall,"Rsmall"),(RMAX,"Rmax")):
    f,ax=plt.subplots(figsize=(8,6)); j=0; allv=[]
    for m in Ms:
        y=(vlos(m,Rv+hh,n)-vlos(m,Rv-hh,n))/(2*hh)*KMSMPC; mn=mean(y); an=100*(y.max()-y.min())/mn
        print("R=%.0f Mpc M=%g: 0° %.1f 90° %.1f 180° %.1f | sky mean %.1f | max %.1f@%.0f min %.1f@%.0f | anisotropy %.0f%%"%(Rv*1000,m,y[0],y[360],y[720],mn,y.max(),ang[y.argmax()],y.min(),ang[y[:721].argmin()],an))
        if m==0: ax.plot(ang,y,"k--",lw=2.5,label="no mass")
        else: ax.plot(ang,y,color=cm[j],lw=3.5,label="%g×"%m); j+=1
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,90)); ax.set_xlabel("angle (deg), 0° = away from CoE"); ax.set_ylabel("dv/dR (km/s/Mpc)")
    ax.set_title("Hubble rate vs angle, R = %.0f Mpc"%(Rv*1000)); ax.legend(frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(.5,-.18),columnspacing=1)
    f.tight_layout(); f.savefig("integ_h70real_slope_%s.png"%tag,dpi=110); plt.close(f)
