"""H0=70 expanding, v=k r launch (only isotropic straight no-mass field), D=20..60 Gpc, R up to 14.3 Gpc. Self-similar mass scaling."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
c=2.998e8; GPC=3.0857e25; KMSMPC=299792.458/1000
b1,mu=0.5866,0.09332; GM1=mu*c**2*(58.2/3)*GPC; k0=c*np.sqrt(b1**2+2*mu*(1/3-1))/58.2/GPC
kp=3*k0; Tp=0.247*(58.2*GPC)**2*k0*14.3*GPC/GM1; kTp=kp*Tp; H=70/KMSMPC*c/GPC; k=H*(1+kTp); T=kTp/k
ang=np.linspace(0,360,721); t=np.radians(ang); n=np.stack([np.cos(t),np.sin(t)],1); w=np.abs(np.sin(t))
plt.rcParams.update({"font.size":17,"axes.titlesize":18,"axes.labelsize":18,"legend.fontsize":14,"xtick.labelsize":16,"ytick.labelsize":16})
DS=[20,30,40,50,60]; Ms=[0,1,4,10]; RES={}
for D in DS:
    Msc=(k/kp)**2*(D/58.2)**3; RM=14.3
    ri=np.linspace(0.005,1.6,3000)*D*GPC; F={}
    for m in Ms:
        GM=m*Msc*GM1; r=ri.copy(); v=k*ri; N=3000; h=T/N; f=lambda r:-GM/np.maximum(r,1e18)**2
        for _ in range(N):
            a1,q1=v,f(r); a2,q2=v+.5*h*q1,f(r+.5*h*a1); a3,q3=v+.5*h*q2,f(r+.5*h*a2); a4,q4=v+h*q3,f(r+h*a3)
            r=r+h/6*(a1+2*a2+2*a3+a4); v=v+h/6*(q1+2*q2+2*q3+q4)
        rf=r/GPC; F[m]=(lambda rf,v:lambda x:np.interp(x,rf,v))(rf,v)
        sv=F[m](np.linspace(D-RM,D+RM,200))/c
        print("D=%d M=%g (%.2e kg): our %.2fc, shells in view %.2f..%.2fc, outward %s"%(D,m,m*Msc*7.5e52,F[m](D)/c,sv.min(),sv.max(),sv.min()>0))
    xo=np.array([D,0.])
    def vl(m,R): xp=xo+R*n; rp=np.linalg.norm(xp,axis=1); return ((F[m](rp)[:,None]*xp/rp[:,None]-F[m](D)*xo/D)*n).sum(1)
    fig,ax=plt.subplots(figsize=(8,6)); cc=plt.cm.viridis(np.linspace(.1,.9,3)); j=0
    for m in Ms:
        for Rv,lab in ((1.0,"1"),(RM,"14.3"),(0.246*D,"RD")):
            y=(vl(m,Rv+.002)-vl(m,Rv-.002))/4*KMSMPC/1000*1000/1000
            y=(vl(m,Rv+.002)-vl(m,Rv-.002))/0.004/1e3/1e3   # km/s per Mpc
            mn=(w*y).sum()/w.sum(); sd=np.sqrt((w*(y-mn)**2).sum()/w.sum()); an=100*(y.max()-y.min())/mn
            rel=np.abs(vl(m,Rv)).max()/c
            RES[(D,m,lab)]=(y[0],y[180],y[360],mn,sd,an,rel)
            print("   R=%s: 0° %.1f 90° %.1f 180° %.1f mean %.1f sd %.1f anis %.0f%% max|v_rel| %.2fc"%(lab,y[0],y[180],y[360],mn,sd,an,rel))
            if lab=="14.3":
                if m==0: ax.plot(ang,y,"k--",lw=2.5,label="no mass")
                else: ax.plot(ang,y,color=cc[j],lw=3.5,label="%g×"%m); j+=1
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,90)); ax.set_xlabel("angle (deg), 0° = away from CoE"); ax.set_ylabel("dv/dR (km/s/Mpc)")
    ax.set_title("Expanding, D = %d Gpc, R = 14.3 Gpc"%D); ax.legend(frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(.5,-.18))
    fig.tight_layout(); fig.savefig("integ_h70D_slope_D%d.png"%D,dpi=110); plt.close(fig)
for idx,(yl,tag) in ((3,("sky-mean dv/dR (km/s/Mpc)","mean")),(5,("anisotropy (max−min)/mean, %","anis"))):
    fig,ax=plt.subplots(figsize=(8,6)); j=0
    for m in Ms:
        y=[RES[(D,m,"14.3")][idx] for D in DS]
        if m==0: ax.plot(DS,y,"k--",lw=2.5,label="no mass")
        else: ax.plot(DS,y,"o-",color=cc[j],lw=3.5,ms=9,label="%g×"%m); j+=1
    ax.set_xlabel("observer distance from CoE D (Gpc)"); ax.set_ylabel(yl); ax.set_title("Expanding, H₀ = 70, R = 14.3 Gpc")
    ax.legend(frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(.5,-.18)); fig.tight_layout(); fig.savefig("integ_h70D_%s_vs_D.png"%tag,dpi=110); plt.close(fig)
